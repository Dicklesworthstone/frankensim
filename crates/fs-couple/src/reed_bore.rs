//! Time-domain beating-reed × viscothermal bore (wave-variable loop).
//!
//! The bore is the reflection function of the TMM input impedance. The
//! reed is the quasistatic Bernoulli valve of Fletcher / McIntyre–
//! Schumacher–Woodhouse: `y = H max(0, 1 − Δp/P_c)`,
//! `U = w y sgn(Δp) √(2|Δp|/ρ)`. There is no instrument crate.

use crate::acoustic_realize::AcousticRealizeError;
use fs_duct::{Duct, LossModel, Termination, input_impedance};
use fs_fft::{C64 as FftC64, Fft};
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
    for i in 0..n {
        let t = i as f64 * dt;
        let p_m = blowing_envelope(reed, t);
        let p_minus_hist = convolve_tail(&r, &p_plus, i);
        let p_plus_i = solve_reed_wave(
            reed,
            gas.density,
            zc,
            r[0],
            p_minus_hist,
            p_m,
        )?;
        p_plus[i] = p_plus_i;
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
    let n_fft = n.next_power_of_two();
    let dt = 1.0 / f64::from(sample_rate_hz);
    let mut buf = vec![FftC64::new(0.0, 0.0); n_fft];
    for (k, bin) in buf.iter_mut().enumerate().take(n_fft / 2 + 1) {
        if k == 0 {
            // DC: a closed reed-facing cylinder reflects as +1; an
            // ideal-open DC is −1. Use the TMM at a tiny ω for the rest.
            let omega = 2.0 * core::f64::consts::PI / (n_fft as f64 * dt);
            let z = input_impedance(physics, gas, omega, LossModel::AllRegime, termination)
                .map_err(AcousticRealizeError::Duct)?
                .impedance;
            let r = reflect(z.re, z.im, zc);
            *bin = FftC64::new(r.0, 0.0);
            continue;
        }
        let omega = core::f64::consts::TAU * k as f64 / (n_fft as f64 * dt);
        let z = input_impedance(physics, gas, omega, LossModel::AllRegime, termination)
            .map_err(AcousticRealizeError::Duct)?
            .impedance;
        let (re, im) = reflect(z.re, z.im, zc);
        *bin = FftC64::new(re, im);
    }
    for k in 1..n_fft / 2 {
        let c = buf[k];
        buf[n_fft - k] = FftC64::new(c.re, -c.im);
    }
    let fft = Fft::new(n_fft);
    let mut scratch = vec![FftC64::new(0.0, 0.0); n_fft];
    fft.inverse(&mut buf, &mut scratch);
    Ok(buf.iter().map(|c| c.re).collect())
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
) -> Result<f64, AcousticRealizeError> {
    // Scalar Newton on p⁺. Residual: U_reed − (p⁺ − p⁻)/Zc.
    let mut p_plus = 0.0;
    for _ in 0..16 {
        let (res, deriv) = reed_residual(reed, rho, zc, r0, p_minus_hist, p_m, p_plus);
        if res.abs() < 1.0e-10 * (1.0 + p_m.abs()) {
            return Ok(p_plus);
        }
        let step = if deriv.abs() < 1.0e-18 {
            -res * 1.0e-3
        } else {
            -res / deriv
        };
        p_plus += step.clamp(-0.5 * reed.closing_pressure_pa, 0.5 * reed.closing_pressure_pa);
        if !p_plus.is_finite() {
            return Err(AcousticRealizeError::Reed {
                what: "reed-bore Newton left the finite set",
            });
        }
    }
    Err(AcousticRealizeError::Reed {
        what: "reed-bore Newton did not converge",
    })
}

fn reed_residual(
    reed: BeatingReed,
    rho: f64,
    zc: f64,
    r0: f64,
    p_minus_hist: f64,
    p_m: f64,
    p_plus: f64,
) -> (f64, f64) {
    let eps = 1.0e-2;
    let f = reed_flow_mismatch(reed, rho, zc, r0, p_minus_hist, p_m, p_plus);
    let fp = reed_flow_mismatch(reed, rho, zc, r0, p_minus_hist, p_m, p_plus + eps);
    (f, (fp - f) / eps)
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


