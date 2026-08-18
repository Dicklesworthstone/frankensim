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
use crate::thin_plate::PlateBank;
use crate::unilateral_contact::slit_contact_force;
use fs_dcontact::Obstacle;
use fs_duct::{Duct, Termination};
use fs_material::gas::GasState;
use fs_math::det;
use fs_scenario::BeatingReed;

/// Realize mouthpiece pressure from a beating reed on a TMM bore.
///
/// Since the block render API landed (music bead 3ez8g.2.1) this is a
/// thin wrapper over [`crate::render::ReedBoreVoice`]: the voice owns the
/// per-sample loop; this one-shot path builds it, renders `n` samples in
/// one block, and swaps the caller's plate bank in and out so refusal
/// semantics (bank untouched on error) and results stay byte-identical to
/// the pre-API code. Two paths, one loop body — they cannot drift.
///
/// # Errors
/// Domain, TMM, or solver refusals.
#[allow(clippy::too_many_arguments)] // one coherent realization record
pub fn realize_reed_bore(
    physics: &Duct,
    gas: &GasState,
    reed: BeatingReed,
    termination: Termination,
    plates: &mut PlateBank,
    listener_m: f64,
    sample_rate_hz: u32,
    n: usize,
    wall: Option<&fs_phs::WallPin>,
) -> Result<Vec<f64>, AcousticRealizeError> {
    let mut voice = crate::render::ReedBoreVoice::new(
        physics,
        gas,
        reed,
        termination,
        PlateBank::default(),
        listener_m,
        sample_rate_hz,
        n,
        wall,
    )?;
    // Move the caller's bank in only after every fallible admission step
    // has succeeded, and back out unconditionally, so a refusal leaves it
    // untouched and success returns the post-render plate state exactly
    // as the pre-API loop did.
    core::mem::swap(voice.plate_bank_mut(), plates);
    let mut p_bore_hist = vec![0.0; n];
    let result = voice.step_block(&mut p_bore_hist);
    core::mem::swap(voice.plate_bank_mut(), plates);
    result?;
    Ok(p_bore_hist)
}

pub(crate) fn map_drive(err: crate::driving_point::DrivingPointError) -> AcousticRealizeError {
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

pub(crate) fn blowing_envelope(reed: BeatingReed, t: f64) -> f64 {
    if reed.attack_s <= 0.0 {
        return reed.blowing_pressure_pa;
    }
    if t >= reed.attack_s {
        return reed.blowing_pressure_pa;
    }
    let x = t / reed.attack_s;
    reed.blowing_pressure_pa * 0.5 * (1.0 - det::cos(core::f64::consts::PI * x))
}

pub(crate) fn reed_structural(reed: BeatingReed) -> (f64, f64) {
    let face = reed.width_m * 0.025;
    let k = if reed.stiffness_n_m > 0.0 {
        reed.stiffness_n_m
    } else {
        reed.closing_pressure_pa * face / reed.rest_opening_m
    };
    let r_damp = 2.0 * 0.35 * det::sqrt((k * reed.mass_kg).max(0.0));
    (k, r_damp)
}

pub(crate) fn aperture_of(reed: BeatingReed) -> BernoulliAperture {
    BernoulliAperture {
        rest_opening_m: reed.rest_opening_m,
        width_m: reed.width_m,
        closing_pressure_pa: reed.closing_pressure_pa,
    }
}

#[allow(clippy::too_many_arguments)] // one coherent reed-step record
pub(crate) fn step_massive_reed(
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

#[allow(clippy::too_many_arguments)] // one coherent junction record
#[allow(clippy::unnecessary_wraps)] // uniform Result surface across the reed solvers
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
            let x = -span + (2.0 * span) * f64::from(k) / 20.0;
            let a = reed_flow_mismatch(reed, rho, zc, r0, p_minus_hist, p_m, x, u_body).abs();
            if a < best_a {
                best_a = a;
                best = x;
            }
        }
        return Ok(best);
    }
    let mut mid = f64::midpoint(lo, hi);
    for _ in 0..48 {
        mid = f64::midpoint(lo, hi);
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

#[allow(clippy::too_many_arguments)] // one coherent junction record
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
