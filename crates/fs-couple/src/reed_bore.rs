//! Compose a [`BernoulliAperture`] onto a [`TravelingWaveLine`].
//!
//! A clarinet reed on a bore is one filling of that composition.
//! Vocal folds, lip reeds, and relief valves are the same objects.

use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::BernoulliAperture;
use crate::traveling_wave_line::{TravelingWaveError, TravelingWaveLine};
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
    let mut line = TravelingWaveLine::from_duct(physics, gas, termination, sample_rate_hz, n, zc)
        .map_err(map_line)?;
    let dt = 1.0 / f64::from(sample_rate_hz);
    let mut reed_y = reed.rest_opening_m;
    let mut reed_v = 0.0;
    let mut p_bore_hist = vec![0.0; n];
    let mut p_plus_prev = 5.0;
    let _ = line.push(p_plus_prev);
    for i in 0..n {
        let p_m = blowing_envelope(reed, i as f64 * dt);
        let p_minus = line.incoming();
        let p_plus = if reed.mass_kg > 0.0 {
            let (pp, y, v) =
                step_massive_reed(reed, gas.density, zc, p_minus, p_m, reed_y, reed_v, dt)?;
            reed_y = y;
            reed_v = v;
            pp
        } else {
            solve_reed_wave(reed, gas.density, zc, 0.0, p_minus, p_m, p_plus_prev)?
        };
        p_plus_prev = p_plus;
        let p_minus_now = line.push(p_plus);
        p_bore_hist[i] = p_plus + p_minus_now;
    }
    Ok(p_bore_hist)
}

fn map_line(err: TravelingWaveError) -> AcousticRealizeError {
    match err {
        TravelingWaveError::Invalid { what } => AcousticRealizeError::Reed { what },
        TravelingWaveError::Duct(e) => AcousticRealizeError::Duct(e),
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
) -> Result<(f64, f64, f64), AcousticRealizeError> {
    let face = reed.width_m * 0.025;
    let k = if reed.stiffness_n_m > 0.0 {
        reed.stiffness_n_m
    } else {
        reed.closing_pressure_pa * face / reed.rest_opening_m
    };
    let r_damp = 2.0 * 0.35 * det::sqrt(k * reed.mass_kg);
    let p_plus = solve_reed_wave(reed, rho, zc, 0.0, p_minus, p_m, 2.0 * p_minus)?;
    let p_bore = p_plus + p_minus;
    let dp = p_m - p_bore;
    let mut acc = (-k * (y - reed.rest_opening_m) - r_damp * v - face * dp) / reed.mass_kg;
    let mut y1 = y + dt * v;
    let mut v1 = v + dt * acc;
    if y1 < 0.0 {
        // Hunt–Crossley beating: k_c |x|^1.5 (1 + 1.5 χ ẋ)
        let depth = -y1;
        let kc = 1.0e7 * reed.width_m;
        let contact = kc * depth.powf(1.5) * (1.0 + 1.5 * 0.6 * v1.max(0.0));
        acc -= contact / reed.mass_kg;
        v1 = v + dt * acc;
        y1 = (y + dt * v1).min(0.0);
        v1 = v1.min(0.0);
    }
    if !y1.is_finite() || !v1.is_finite() {
        return Err(AcousticRealizeError::Reed {
            what: "massive reed left the finite set",
        });
    }
    Ok((p_plus, y1, v1))
}

fn solve_reed_wave(
    reed: BeatingReed,
    rho: f64,
    zc: f64,
    r0: f64,
    p_minus_hist: f64,
    p_m: f64,
    guess: f64,
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
    let mut f_lo = reed_flow_mismatch(reed, rho, zc, r0, p_minus_hist, p_m, lo);
    let mut f_hi = reed_flow_mismatch(reed, rho, zc, r0, p_minus_hist, p_m, hi);
    if f_lo * f_hi > 0.0 {
        lo = -span;
        hi = span;
        f_lo = reed_flow_mismatch(reed, rho, zc, r0, p_minus_hist, p_m, lo);
        f_hi = reed_flow_mismatch(reed, rho, zc, r0, p_minus_hist, p_m, hi);
    }
    if f_lo * f_hi > 0.0 {
        let mut best = guess;
        let mut best_a = reed_flow_mismatch(reed, rho, zc, r0, p_minus_hist, p_m, guess).abs();
        for k in 0..21 {
            let x = -span + (2.0 * span) * k as f64 / 20.0;
            let a = reed_flow_mismatch(reed, rho, zc, r0, p_minus_hist, p_m, x).abs();
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
        let f_mid = reed_flow_mismatch(reed, rho, zc, r0, p_minus_hist, p_m, mid);
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
) -> f64 {
    let p_minus = p_minus_hist + r0 * p_plus;
    let p_bore = p_plus + p_minus;
    let dp = p_m - p_bore;
    let flow = aperture_of(reed).volume_flow(dp, rho);
    let u_wave = (p_plus - p_minus) / zc;
    flow - u_wave
}
