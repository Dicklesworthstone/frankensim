//! Compose a [`BernoulliAperture`] onto a characteristic line.
//!
//! The time port is the TMM reflectance FIR
//! ([`crate::driving_point`]). Reed lay is an [`fs_dcontact`] obstacle
//! whose Hunt–Crossley `χ` is matched to the existing viscous damper
//! (not a private impact law). Observer absorption is applied once
//! by realize (ISO 9613-1). Vocal folds, lip reeds, and relief
//! valves are the same objects.

use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::BernoulliAperture;
use crate::driving_point::characteristic_line;
use crate::thin_plate::PlateBank;
use crate::unilateral_contact::{slit_contact_force, slit_lay};
use fs_dcontact::Obstacle;
use fs_duct::{Duct, Termination};
use fs_material::gas::GasState;
use fs_math::det;
use fs_scenario::BeatingReed;

/// Realize mouthpiece pressure from a beating reed on a TMM bore.
///
/// # Errors
/// Domain, TMM, or solver refusals.
pub fn realize_reed_bore(
    physics: &Duct,
    gas: &GasState,
    reed: BeatingReed,
    termination: Termination,
    plates: &mut PlateBank,
    listener_m: f64,
    sample_rate_hz: u32,
    n: usize,
) -> Result<Vec<f64>, AcousticRealizeError> {
    if !(reed.rest_opening_m > 0.0
        && reed.width_m > 0.0
        && reed.closing_pressure_pa > 0.0
        && reed.blowing_pressure_pa >= 0.0
        && reed.attack_s >= 0.0
        && reed.mass_kg >= 0.0
        && reed.stiffness_n_m >= 0.0
        && reed.rest_opening_m.is_finite())
    {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "reed parameters must be physical and finite",
        });
    }
    let inlet_r = physics
        .segments
        .first()
        .ok_or(AcousticRealizeError::InvalidDescription {
            what: "duct has no segments",
        })?
        .outlet_radius();
    let area_bore = core::f64::consts::PI * inlet_r * inlet_r;
    let zc = gas.density * gas.sound_speed / area_bore;
    let mut line =
        characteristic_line(physics, gas, termination, sample_rate_hz, n, zc).map_err(map_drive)?;
    let dt = 1.0 / f64::from(sample_rate_hz);
    let mut reed_y = reed.rest_opening_m;
    let mut reed_v = 0.0;
    let lay = if reed.mass_kg > 0.0 {
        let k_lay = 1.0e7 * reed.width_m;
        let (_k, r_damp) = reed_structural(reed);
        let chi = r_damp / (k_lay * reed.rest_opening_m * reed.rest_opening_m).max(1.0e-18);
        Some(
            slit_lay(k_lay, 2.0)
                .and_then(|o| o.with_internal_loss(chi))
                .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?,
        )
    } else {
        None
    };
    let mut p_bore_hist = vec![0.0; n];
    let mut p_plus_prev = 5.0;
    let mut f_jet_prev = 0.0;
    let _ = line
        .push(p_plus_prev)
        .map_err(|_| AcousticRealizeError::Reed {
            what: "characteristic line left the finite set",
        })?;
    for i in 0..n {
        let p_m = blowing_envelope(reed, i as f64 * dt);
        let p_minus = line.incoming();
        let u_body = plates.volume_velocity();
        let p_plus = if reed.mass_kg > 0.0 {
            let (pp, y, v) = step_massive_reed(
                reed,
                gas.density,
                zc,
                p_minus,
                p_m,
                reed_y,
                reed_v,
                dt,
                u_body,
                lay.as_ref(),
            )?;
            reed_y = y;
            reed_v = v;
            pp
        } else {
            solve_reed_wave(
                reed,
                gas.density,
                zc,
                0.0,
                p_minus,
                p_m,
                p_plus_prev,
                u_body,
            )?
        };
        p_plus_prev = p_plus;
        let p_minus_now = line.push(p_plus).map_err(|_| AcousticRealizeError::Reed {
            what: "characteristic line left the finite set",
        })?;
        let p_bore = p_plus + p_minus_now;
        let mut p_obs = p_bore;
        p_obs += plates.drive_and_radiate(p_bore * area_bore, dt, gas.density, listener_m)?;
        // Compact far-field dipole of the slit force. This is the
        // same Green's family as the string dipole, not a 3-D jet
        // and not fs-aeroac 2-D Curle (cylindrical Hankel).
        let opening = if reed.mass_kg > 0.0 {
            reed_y.max(0.0)
        } else {
            aperture_of(reed).opening_m(p_m - p_bore).max(0.0)
        };
        let f_jet = (p_m - p_bore) * reed.width_m * opening;
        let dfdt = (f_jet - f_jet_prev) / dt;
        f_jet_prev = f_jet;
        p_obs += gas.density * dfdt / (4.0 * core::f64::consts::PI * listener_m * gas.sound_speed);
        p_bore_hist[i] = p_obs;
    }
    Ok(p_bore_hist)
}

fn map_drive(err: crate::driving_point::DrivingPointError) -> AcousticRealizeError {
    match err {
        crate::driving_point::DrivingPointError::Invalid { what } => {
            AcousticRealizeError::Reed { what }
        }
        crate::driving_point::DrivingPointError::Duct(d) => AcousticRealizeError::Duct(d),
        crate::driving_point::DrivingPointError::Realize(_) => AcousticRealizeError::Reed {
            what: "characteristic realization refused",
        },
        crate::driving_point::DrivingPointError::Discrete(_) => AcousticRealizeError::Reed {
            what: "characteristic line left the finite set",
        },
    }
}

fn blowing_envelope(reed: BeatingReed, t: f64) -> f64 {
    if reed.attack_s <= 0.0 {
        return reed.blowing_pressure_pa;
    }
    if t >= reed.attack_s {
        return reed.blowing_pressure_pa;
    }
    let x = t / reed.attack_s;
    reed.blowing_pressure_pa * 0.5 * (1.0 - det::cos(core::f64::consts::PI * x))
}

fn reed_structural(reed: BeatingReed) -> (f64, f64) {
    let face = reed.width_m * 0.025;
    let k = if reed.stiffness_n_m > 0.0 {
        reed.stiffness_n_m
    } else {
        reed.closing_pressure_pa * face / reed.rest_opening_m
    };
    let r_damp = 2.0 * 0.35 * det::sqrt((k * reed.mass_kg).max(0.0));
    (k, r_damp)
}

fn aperture_of(reed: BeatingReed) -> BernoulliAperture {
    BernoulliAperture {
        rest_opening_m: reed.rest_opening_m,
        width_m: reed.width_m,
        closing_pressure_pa: reed.closing_pressure_pa,
    }
}

fn step_massive_reed(
    reed: BeatingReed,
    rho: f64,
    zc: f64,
    p_minus: f64,
    p_m: f64,
    y: f64,
    v: f64,
    dt: f64,
    u_body: f64,
    lay: Option<&Obstacle>,
) -> Result<(f64, f64, f64), AcousticRealizeError> {
    let face = reed.width_m * 0.025;
    let (k, r_damp) = reed_structural(reed);
    let p_plus = solve_reed_wave(reed, rho, zc, 0.0, p_minus, p_m, 2.0 * p_minus, u_body)?;
    let p_bore = p_plus + p_minus;
    let dp = p_m - p_bore;
    let mut acc = (-k * (y - reed.rest_opening_m) - r_damp * v - face * dp) / reed.mass_kg;
    let mut y1 = y + dt * v;
    let mut v1 = v + dt * acc;
    if let Some(obstacle) = lay {
        let contact = slit_contact_force(obstacle, y1)
            .map_err(|e| AcousticRealizeError::Nonlinear(e.to_string()))?;
        let hc = obstacle.dissipative_modal_forces(1, &[y1, v], &[v])[0];
        acc += (contact + hc) / reed.mass_kg;
        v1 = v + dt * acc;
        y1 = y + dt * v1;
    }
    if !y1.is_finite() || !v1.is_finite() {
        return Err(AcousticRealizeError::Reed {
            what: "massive reed left the finite set",
        });
    }
    Ok((p_plus, y1, v1))
}

pub(crate) fn solve_reed_wave(
    reed: BeatingReed,
    rho: f64,
    zc: f64,
    r0: f64,
    p_minus_hist: f64,
    p_m: f64,
    guess: f64,
    u_body: f64,
) -> Result<f64, AcousticRealizeError> {
    let denom = (1.0 - r0).clamp(-0.999, 0.999);
    let closed_plus = p_minus_hist / denom;
    let p_bore_closed = closed_plus + (p_minus_hist + r0 * closed_plus);
    if p_m - p_bore_closed >= reed.closing_pressure_pa {
        return Ok(closed_plus);
    }
    let span = (2.0 * reed.closing_pressure_pa)
        .max(2.0 * p_m.abs())
        .max(1.0);
    let mut lo = guess - span;
    let mut hi = guess + span;
    let mut f_lo = reed_flow_mismatch(reed, rho, zc, r0, p_minus_hist, p_m, lo, u_body);
    let mut f_hi = reed_flow_mismatch(reed, rho, zc, r0, p_minus_hist, p_m, hi, u_body);
    if f_lo * f_hi > 0.0 {
        lo = -span;
        hi = span;
        f_lo = reed_flow_mismatch(reed, rho, zc, r0, p_minus_hist, p_m, lo, u_body);
        f_hi = reed_flow_mismatch(reed, rho, zc, r0, p_minus_hist, p_m, hi, u_body);
    }
    if f_lo * f_hi > 0.0 {
        let mut best = guess;
        let mut best_a =
            reed_flow_mismatch(reed, rho, zc, r0, p_minus_hist, p_m, guess, u_body).abs();
        for k in 0..21 {
            let x = -span + (2.0 * span) * k as f64 / 20.0;
            let a = reed_flow_mismatch(reed, rho, zc, r0, p_minus_hist, p_m, x, u_body).abs();
            if a < best_a {
                best_a = a;
                best = x;
            }
        }
        return Ok(best);
    }
    let mut mid = 0.5 * (lo + hi);
    for _ in 0..48 {
        mid = 0.5 * (lo + hi);
        let f_mid = reed_flow_mismatch(reed, rho, zc, r0, p_minus_hist, p_m, mid, u_body);
        if f_mid.abs() < 1.0e-8 * (1.0 + p_m.abs()) {
            return Ok(mid);
        }
        if f_lo * f_mid <= 0.0 {
            hi = mid;
            f_hi = f_mid;
        } else {
            lo = mid;
            f_lo = f_mid;
        }
    }
    let _ = f_hi;
    Ok(mid)
}

fn reed_flow_mismatch(
    reed: BeatingReed,
    rho: f64,
    zc: f64,
    r0: f64,
    p_minus_hist: f64,
    p_m: f64,
    p_plus: f64,
    u_body: f64,
) -> f64 {
    let p_minus = p_minus_hist + r0 * p_plus;
    let p_bore = p_plus + p_minus;
    let dp = p_m - p_bore;
    let flow = aperture_of(reed).volume_flow(dp, rho);
    let u_wave = (p_plus - p_minus) / zc;
    flow + u_body - u_wave
}
