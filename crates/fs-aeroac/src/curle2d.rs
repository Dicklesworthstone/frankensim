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

/// Locked Strouhal number from the pinned 2-D CentralMoment ladder
/// (`slot_half = 6`, the executed Re sweep).
///
/// Log-linear interpolation in Reynolds number between neighbouring
/// rows. Outside the table the nearest endpoint is returned — that is
/// an observation bound, not an extrapolation claim. The jump near
/// Re ~ 10³ is in the recorded data (stage change), not smoothed away.
/// This does **not** mint 3-D broadband: every row is tonal.
#[must_use]
pub fn strouhal_at_reynolds(reynolds: f64) -> Option<f64> {
    if !(reynolds > 0.0 && reynolds.is_finite()) {
        return None;
    }
    let mut rows: Vec<(f64, f64)> = crate::regime::PINNED_2D_CENTRAL_MOMENT_TONAL
        .iter()
        .filter(|row| (row.slot_half - 6.0).abs() < 1.0e-12 && row.ran_in_regime)
        .map(|row| (row.reynolds, row.strouhal))
        .collect();
    if rows.is_empty() {
        return None;
    }
    rows.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(core::cmp::Ordering::Equal));
    if reynolds <= rows[0].0 {
        return Some(rows[0].1);
    }
    if reynolds >= rows[rows.len() - 1].0 {
        return Some(rows[rows.len() - 1].1);
    }
    for pair in rows.windows(2) {
        let (re0, st0) = pair[0];
        let (re1, st1) = pair[1];
        if reynolds >= re0 && reynolds <= re1 {
            let ln0 = det::ln(re0);
            let t = (det::ln(reynolds) - ln0) / (det::ln(re1) - ln0);
            return Some(st0 + t * (st1 - st0));
        }
    }
    None
}

/// Tonal lock frequency [Hz]: `f = St(Re) · U / δ`.
///
/// `slot_width` is the same `delta` the pinned table uses (slot half
/// in lattice units, or the physical slot scale the caller binds).
#[must_use]
pub fn tonal_lock_frequency(reynolds: f64, u_jet: f64, slot_width: f64) -> Option<f64> {
    if !(u_jet > 0.0 && slot_width > 0.0 && u_jet.is_finite() && slot_width.is_finite()) {
        return None;
    }
    strouhal_at_reynolds(reynolds).map(|st| st * u_jet / slot_width)
}

/// Compact 2-D Curle dipole of a locked slot-jet tone at an observer.
///
/// Force per unit span is the caller-authored lift coefficient times
/// `½ ρ U² δ`. Shape and scaling only — see [`crate::SCOPE_STATEMENT`].
/// This is observer-side radiation of the existing 2-D lock, not a
/// 3-D jet-noise spectrum.
///
/// # Errors
/// [`AeroacError`] from the dipole kernel, or missing Strouhal lock.
#[allow(clippy::too_many_arguments)] // one coherent observer record
pub fn tonal_dipole_observer(
    density: f64,
    u_jet: f64,
    slot_width: f64,
    reynolds: f64,
    lift_coeff: f64,
    sound_speed: f64,
    observer: [f64; 2],
    source: [f64; 2],
) -> Result<DipoleField, AeroacError> {
    if ![density, u_jet, slot_width, lift_coeff, sound_speed]
        .iter()
        .all(|v| v.is_finite() && *v > 0.0)
    {
        return Err(AeroacError::NonFinite {
            what: "tonal dipole inputs",
        });
    }
    let freq =
        tonal_lock_frequency(reynolds, u_jet, slot_width).ok_or(AeroacError::InvalidParameter {
            what: "no pinned 2D Strouhal lock at this Reynolds number",
        })?;
    let omega = 2.0 * core::f64::consts::PI * freq;
    let k = omega / sound_speed;
    let force_amp = 0.5 * density * u_jet * u_jet * slot_width * lift_coeff;
    dipole_pressure(
        [C64::new(force_amp, 0.0), C64::new(0.0, 0.0)],
        k,
        observer,
        source,
    )
}

/// Observer-only amplitude modulation by a locked 2-D tone.
///
/// Multiplies an existing pressure history by
/// `1 + depth · sin(2π f t)`. It does **not** enter a reed/bore lock
/// loop. `depth` is clamped to `[0, 1]`.
pub fn modulate_observer_by_tone(pressure: &mut [f64], dt: f64, frequency_hz: f64, depth: f64) {
    if !(dt > 0.0 && frequency_hz > 0.0 && dt.is_finite() && frequency_hz.is_finite()) {
        return;
    }
    let depth = depth.clamp(0.0, 1.0);
    let omega = 2.0 * core::f64::consts::PI * frequency_hz;
    for (i, sample) in pressure.iter_mut().enumerate() {
        let phase = omega * dt * i as f64;
        *sample *= 1.0 + depth * det::sin(phase);
    }
}
