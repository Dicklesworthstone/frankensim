//! Coupled Fourier heat diffusion and local molecular energy reservoirs.
//!
//! For each original thermal mode, cp(s) = cp_inf + sum Cj/(1+s*tauj)
//! gives beta(s) = [1-R sum w*s/(s*cp(s)+cp_eq*lambda)]/p.
//! Keep cp_eq*lambda: lambda already contains the material conductivity.
//!
//! Stack capacity C=diag(cp_inf,Cj) and the positive thermal/reservoir
//! conductance L; q selects each temperature with sqrt(w). Woodbury gives
//! K(s)=p+p*R*s*q' [s*(C-R*q*q')+L]^-1 q. Since C-R*q*q' is SPD,
//! ONE symmetric eigensolve yields positive Maxwell compression arms.
//! The existing mechanical/thermal network then owns storage and heat.
//!
//! This retains the original frequency-independent Fourier conductivity and
//! local Debye reservoirs. It supplies no separate diffusion of molecular
//! internal energies, finite wall conductivity, or nonuniform acoustic pressure.
//! The pressure-driven energy equation and equation of state are the standard
//! linear thermoviscous equations; the positive realization above is their
//! algebraic reduction with this repository's original molecular law:
//! https://www.comsol.com/blogs/theory-of-thermoviscous-acoustics-thermal-and-viscous-losses/

use fs_material::relaxation::RelaxationState;
use fs_math::{c64::C64, det};

use super::spectrum::{self, MAX_COORDINATES, Spectrum};

#[derive(Clone, Copy, Debug)]
pub(super) struct ThermalMode {
    pub weight: f64,
    pub rate_per_s: f64,
}

fn positive(x: f64) -> bool {
    x.is_finite() && x > 0.0
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), &'static str> {
    if cancelled() {
        Err("molecular cavity thermal preparation cancelled")
    } else {
        Ok(())
    }
}

/// Resolve the original thermal/molecular law, including its dissipative part.
/// The 33-point band check verifies this algebraic realization, separately
/// from the thermal owner's equilibrium Bessel approximation admission.
pub(super) fn realize(
    law: &RelaxationState,
    modes: &[ThermalMode],
    band_rad_s: [f64; 2],
    maximum_coordinates: usize,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Spectrum, &'static str> {
    poll(cancelled)?;
    if maximum_coordinates == 0
        || maximum_coordinates > MAX_COORDINATES
        || modes.is_empty()
        || modes
            .iter()
            .any(|m| !positive(m.weight) || !positive(m.rate_per_s))
        || !positive(band_rad_s[0])
        || !positive(band_rad_s[1])
        || band_rad_s[0] >= band_rad_s[1]
    {
        return Err("molecular cavity needs positive thermal modes, band and bounded coordinates");
    }
    let reservoirs: Vec<_> = law
        .processes()
        .iter()
        .filter(|p| p.heat_capacity_j_kg_k > 0.0)
        .collect();
    let width = reservoirs
        .len()
        .checked_add(1)
        .ok_or("molecular cavity coordinate count overflow")?;
    let n = modes
        .len()
        .checked_mul(width)
        .filter(|&n| n <= maximum_coordinates)
        .ok_or("molecular cavity thermal coordinate budget exceeded")?;
    let gas = law.gas();
    let r = gas.specific_gas_constant;
    let cp_inf = law.frozen_cv_j_kg_k() + r;
    let total_weight: f64 = modes.iter().map(|m| m.weight).sum();
    if (total_weight - 1.0).abs() > 2e-12 || !positive(cp_inf) {
        return Err("molecular cavity thermal measure lost its storage moment");
    }
    let mut a = vec![0.0; n * n];
    let mut b = vec![0.0; n];
    for (block, mode) in modes.iter().enumerate() {
        poll(cancelled)?;
        let i = block * width;
        b[i] = det::sqrt(mode.weight / cp_inf);
        a[i * n + i] = mode.rate_per_s * (gas.specific_heat_cp / cp_inf);
        for (j, process) in reservoirs.iter().enumerate() {
            let k = i + j + 1;
            let rate = process.relaxation_time_s.recip();
            let capacity_ratio = process.heat_capacity_j_kg_k / cp_inf;
            if !positive(rate) || !positive(capacity_ratio) {
                return Err("molecular cavity reservoir scales are unrepresentable");
            }
            let coupling = -rate * det::sqrt(capacity_ratio);
            a[i * n + i] += rate * capacity_ratio;
            a[i * n + k] = coupling;
            a[k * n + i] = coupling;
            a[k * n + k] = rate;
        }
    }
    let result =
        spectrum::realize(gas.pressure, r, a, &b, cancelled).map_err(|reason| match reason {
            "cavity caloric preparation cancelled" => {
                "molecular cavity thermal preparation cancelled"
            }
            other => other,
        })?;
    for sample in 0..=32 {
        poll(cancelled)?;
        let omega = if sample == 32 {
            band_rad_s[1]
        } else {
            band_rad_s[0] * det::pow(band_rad_s[1] / band_rad_s[0], sample as f64 / 32.0)
        };
        let cp = law
            .caloric_response(omega)
            .map_err(|_| "molecular cavity caloric response is unrepresentable")?
            .cp_j_kg_k
            .conj();
        let temperature = modes.iter().fold(C64::ZERO, |sum, mode| {
            let denominator = cp + C64::new(0.0, -gas.specific_heat_cp * (mode.rate_per_s / omega));
            sum + denominator.recip().scale(mode.weight)
        });
        let exact = (C64::ONE - temperature.scale(r))
            .recip()
            .scale(gas.pressure);
        let actual = result.bulk_at(gas.pressure, omega);
        if !positive(exact.re)
            || !positive(exact.im)
            || !positive(actual.re)
            || !positive(actual.im)
            || (actual - exact).abs() > 1e-8 * exact.abs()
            || (actual.im - exact.im).abs() > 1e-8 * exact.im
        {
            return Err(
                "molecular cavity thermal realization exceeds its sampled algebraic tolerance",
            );
        }
    }
    Ok(result)
}

#[cfg(test)]
mod tests;
