//! Exact damaged-area compression measure on the ORIGINAL damage path.
//!
//! A rigid-friction mixture may use `d p_c` without changing the bond's free
//! energy or introducing an elastic friction spring. This module supplies only
//! that path measure; it neither chooses a coefficient nor advances a second
//! damage history. In particular total normal traction (which may include
//! tensile bonds and viscosity) is not a substitute for compression pressure.

use super::*;

impl<'a> CohesivePoint<'a> {
    /// Original timed material trial, its mean traction [Pa], and
    /// `integral_0^1 d(s) K_compression max(-g_n(s),0) ds` [Pa].
    ///
    /// All three outputs share the same original irreversible damage path,
    /// including first frozen-BK initiation, unloading and complete failure.
    /// The extra scalar is nonnegative and adds no energy or material state.
    /// It supports a declared damaged-area rigid-Coulomb mixture; this is not
    /// an implementation of a separate elastic friction internal variable.
    /// Analytic, allocation-free, bounded work, with no endpoint-damage shortcut.
    pub fn trial_over_time_with_damage_weighted_compression(
        &self,
        jump_m: [f64; 3],
        dt_s: f64,
        energy_tolerance_j_m2: f64,
    ) -> Result<(CohesiveRateStep<'a>, [f64; 3], f64), CohesiveError> {
        let (step, mean) =
            self.trial_over_time_with_mean_traction(jump_m, dt_s, energy_tolerance_j_m2)?;
        let pressure = damaged_compression(self, &step.equilibrium.point)?;
        Ok((step, mean, pressure))
    }
}

fn damaged_compression(
    old: &CohesivePoint<'_>,
    end: &CohesivePoint<'_>,
) -> Result<f64, CohesiveError> {
    let compression = positive_mean(-old.jump[0], -end.jump[0]);
    if compression == 0.0 {
        return Ok(0.0);
    }
    let damage = old.damage();
    if damage == 1.0 {
        return finite(old.law.compression_stiffness * compression);
    }
    let (weights, onset, failure, threshold, lower) = match old.law.mixed_mode {
        Some(p) => {
            let Some(h) = end.mixed_history else {
                return Ok(0.0);
            };
            let lower = if old.mixed_history.is_none() {
                super::super::mixed_mode::initiate(p, old.jump, end.jump)?.1
            } else {
                0.0
            };
            (
                [
                    det::sqrt(p.normal_stiffness_pa_m),
                    det::sqrt(p.tangential_stiffness_pa_m),
                    det::sqrt(p.tangential_stiffness_pa_m),
                ],
                h.onset_sqrt_energy,
                h.failure_sqrt_energy,
                old.maximum.max(h.onset_sqrt_energy),
                lower,
            )
        }
        None => (
            [1.0; 3],
            old.law.onset,
            old.law.failure,
            old.maximum.max(old.law.onset),
            0.0,
        ),
    };
    let weighted_old = core::array::from_fn(|i| weights[i] * old.jump[i]);
    let weighted_end = core::array::from_fn(|i| weights[i] * end.jump[i]);
    let end_radius = norm(active(weighted_end))?;
    if end_radius <= threshold {
        return finite(old.law.compression_stiffness * (damage * compression));
    }
    let start = outward_crossing(weighted_old, weighted_end, threshold, lower)?;
    let finish = if end_radius >= failure {
        outward_crossing(weighted_old, weighted_end, failure, start)?
    } else {
        1.0
    };
    let at = |s: f64| core::array::from_fn(|i| (1.0 - s) * old.jump[i] + s * end.jump[i]);
    let mut mean = start * damage * positive_mean(-old.jump[0], -at(start)[0]);
    let switch = normal_switch(old.jump[0], end.jump[0]);
    let cut = if switch > start && switch < finish {
        switch
    } else {
        finish
    };
    for (a, b) in [(start, cut), (cut, finish)] {
        if b > a {
            let x = at(a);
            let y = at(b);
            let p = [(-x[0]).max(0.0), (-y[0]).max(0.0)];
            if p != [0.0; 2] {
                let u = core::array::from_fn(|i| weights[i] * active(x)[i]);
                let v = core::array::from_fn(|i| weights[i] * active(y)[i]);
                mean += (b - a) * softening_pressure(u, v, p, onset, failure)?;
            }
        }
    }
    mean += (1.0 - finish) * positive_mean(-at(finish)[0], -end.jump[0]);
    // No negative loss or damaged area above the complete reference area is
    // repaired by clipping. A numerically inadmissible measure is refused.
    if !mean.is_finite() || mean < 0.0 || mean > compression * (1.0 + 64.0 * f64::EPSILON) {
        return Err(CohesiveError::Invalid("damaged compression path measure"));
    }
    finite(old.law.compression_stiffness * mean)
}

// Mean d(r(t))*p(t), with affine nonnegative compression p and weighted active
// jump q on one outward softening segment. Damage is f/(f-o)*(1-o/|q|).
fn softening_pressure(
    a: [f64; 3],
    b: [f64; 3],
    p: [f64; 2],
    onset: f64,
    failure: f64,
) -> Result<f64, CohesiveError> {
    let r0 = positive(norm(a)?)?;
    let u = a.map(|x| x / r0);
    let d = core::array::from_fn(|i| (b[i] - a[i]) / r0);
    let ud = dot(u, d);
    let dd = dot(d, d);
    let factor = failure / (failure - onset);
    let average = 0.5 * p[0] + 0.5 * p[1];
    let change = p[1] - p[0];
    let result = if 2.0 * ud.abs() + dd <= 0.001 {
        // Original inverse-radius moments, kept as small differences so a
        // newly damaged area is not obtained by subtracting two pressures.
        let mut moment = [0.0; 2];
        let mut binomial = 1.0;
        for k in 1..=8 {
            binomial *= (-0.5 - f64::from(k - 1)) / f64::from(k);
            let mut choose = 1.0;
            for j in 0..=k {
                if j > 0 {
                    choose *= f64::from(k - j + 1) / f64::from(j);
                }
                let coefficient = choose * powi(2.0 * ud, k - j) * powi(dd, j);
                moment[0] += binomial * coefficient / f64::from(k + j + 1);
                moment[1] += binomial * coefficient / f64::from(k + j + 2);
            }
        }
        (factor * ((r0 - onset) / r0)) * average
            - (factor * (onset / r0)) * (p[0] * moment[0] + change * moment[1])
    } else {
        let length = norm(d)?;
        let direction = d.map(|x| x / length);
        let projection = dot(u, direction);
        let end_radius = norm(core::array::from_fn(|i| u[i] + d[i]))?;
        let radial_rate = (length + 2.0 * projection) / (end_radius + 1.0);
        let logarithm = if length >= 1.0 {
            det::ln(length) + det::ln(1.0 / length + (1.0 + radial_rate) / (1.0 + projection))
        } else {
            det::ln(1.0 + length * (1.0 + radial_rate) / (1.0 + projection))
        };
        let j0 = logarithm / length;
        let j1 = (radial_rate - projection * j0) / length;
        factor * (average - (onset / r0) * (p[0] * j0 + change * j1))
    };
    finite(result)
}
