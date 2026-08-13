//! Time-domain beating-reed × viscothermal bore (wave-variable loop).
//!
//! The bore is the reflection function of the TMM input impedance. The
//! reed is the quasistatic Bernoulli valve of Fletcher / McIntyre–
//! Schumacher–Woodhouse: `y = H max(0, 1 − Δp/P_c)`,
//! `U = w y sgn(Δp) √(2|Δp|/ρ)`. There is no instrument crate.

use crate::acoustic_realize::AcousticRealizeError;
use fs_duct::{Duct, LossModel, Termination, input_impedance};
use fs_material::gas::GasState;
use fs_math::det;
use fs_scenario::BeatingReed;

/// Realize mouthpiece pressure from a beating reed on a TMM bore.
///
/// # Errors
/// Domain, TMM, or Newton refusals.
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
        && reed.rest_opening_m.is_finite()
        && reed.width_m.is_finite()
        && reed.closing_pressure_pa.is_finite()
        && reed.blowing_pressure_pa.is_finite())
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
    let area = core::f64::consts::PI * inlet_r * inlet_r;
    let zc = gas.density * gas.sound_speed / area;
    let r = reflection_function(physics, gas, termination, sample_rate_hz, n, zc)?;
    let dt = 1.0 / f64::from(sample_rate_hz);
    let mut p_plus = vec![0.0; n];
    let mut p_bore = vec![0.0; n];
    let mut prev = 5.0;
    p_plus[0] = 5.0;
    for i in 0..n {
        let t = i as f64 * dt;
        let p_m = blowing_envelope(reed, t);
        let p_minus_hist = convolve_tail(&r, &p_plus, i);
        let p_plus_i = solve_reed_wave(reed, gas.density, zc, r[0], p_minus_hist, p_m, prev)?;
        p_plus[i] = p_plus_i;
        prev = p_plus_i;
        let p_minus = p_minus_hist + r[0] * p_plus_i;
        p_bore[i] = p_plus_i + p_minus;
    }
    Ok(p_bore)
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

fn reflection_function(
    physics: &Duct,
    gas: &GasState,
    termination: Termination,
    sample_rate_hz: u32,
    n: usize,
    zc: f64,
) -> Result<Vec<f64>, AcousticRealizeError> {
    // Round-trip delay 2L/c of the cylinder chain, polarity from the
    // termination, magnitude from the viscothermal TMM at the quarter-
    // wave. This is the exact lossless kernel −δ(t − 2L/c) (open) or
    // +δ(t − 2L/c) (closed), with |R| taken from AllRegime Z_in.
    let length: f64 = physics
        .segments
        .iter()
        .map(|s| match *s {
            fs_duct::Segment::Cylinder { length, .. } | fs_duct::Segment::Cone { length, .. } => {
                length
            }
            fs_duct::Segment::ToneHole { .. } => 0.0,
        })
        .sum();
    if !(length > 0.0) {
        return Err(AcousticRealizeError::InvalidDescription {
            what: "reed-bore needs positive cylinder length",
        });
    }
    let dt = 1.0 / f64::from(sample_rate_hz);
    let delay = (2.0 * length / gas.sound_speed / dt).round();
    if !(delay >= 2.0 && delay < n as f64 - 1.0) {
        return Err(AcousticRealizeError::Reed {
            what: "round-trip delay does not fit the realized history",
        });
    }
    let delay = delay as usize;
    let omega0 = core::f64::consts::PI * gas.sound_speed / (2.0 * length);
    let z = input_impedance(physics, gas, omega0, LossModel::AllRegime, termination)
        .map_err(AcousticRealizeError::Duct)?
        .impedance;
    let (rr, ri) = reflect(z.re, z.im, zc);
    let mag = det::sqrt(rr * rr + ri * ri).clamp(0.2, 0.99);
    let sign = match termination {
        Termination::Closed => 1.0,
        _ => -1.0,
    };
    let n_fft = n.next_power_of_two();
    let mut r = vec![0.0; n_fft];
    r[delay] = sign * mag;
    Ok(r)
}

fn reflect(zr: f64, zi: f64, zc: f64) -> (f64, f64) {
    // R = (Z − Zc)/(Z + Zc)
    let dr = zr + zc;
    let di = zi;
    let den = dr * dr + di * di;
    if den < 1.0e-30 {
        return (0.0, 0.0);
    }
    let nr = zr - zc;
    ((nr * dr + zi * di) / den, (zi * dr - nr * di) / den)
}

fn convolve_tail(r: &[f64], p_plus: &[f64], n: usize) -> f64 {
    let mut acc = 0.0;
    let last = r.len().min(n + 1);
    for k in 1..last {
        acc += r[k] * p_plus[n - k];
    }
    acc
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
    let open = (1.0 - dp / reed.closing_pressure_pa).clamp(0.0, 1.0);
    let y = reed.rest_opening_m * open;
    let flow = if y <= 0.0 || dp.abs() < 1.0e-12 {
        0.0
    } else {
        reed.width_m * y * dp.signum() * det::sqrt(2.0 * dp.abs() / rho)
    };
    let u_wave = (p_plus - p_minus) / zc;
    flow - u_wave
}
