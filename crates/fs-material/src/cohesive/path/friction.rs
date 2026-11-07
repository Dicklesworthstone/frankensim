//! Exact damaged-area compression measure on the ORIGINAL damage path.
//!
//! A rigid-friction mixture may use `d p_c` without changing the bond's free
//! energy or introducing an elastic friction spring. This module supplies only
//! that path measure; it neither chooses a coefficient nor advances a second
//! damage history. In particular total normal traction (which may include
//! tensile bonds and viscosity) is not a substitute for compression pressure.

use super::*;
use fs_ad::Dual64;

type D = Dual64<3>;

impl<'a> CohesivePoint<'a> {
    /// Original timed material trial, its mean traction [Pa], and
    /// `integral_0^1 d(s) K_compression max(-g_n(s),0) ds` [Pa].
    ///
    /// All three outputs share the same original irreversible damage path,
    /// including first selected-envelope initiation, unloading and complete failure.
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

    /// Original timed trial, mean traction and damaged-area compression, plus
    /// their optional endpoint-jump derivatives at fixed accepted history and dt.
    /// The fourth output is d(mean traction)/d(jump) [Pa/m]; the fifth is
    /// d(integral d p_compression)/d(jump) [Pa/m]. These are separate derivatives:
    /// bond viscosity contributes to mean traction but not compression pressure.
    ///
    /// The pressure gradient differentiates the SAME analytic pressure integral,
    /// including moving initiation, rupture and normal-sign crossings. A first
    /// mixed-mode envelope retains its original initiation sensitivities; an accepted
    /// envelope stays frozen. No quadrature, finite difference or new history
    /// is introduced. Ordinary pressure callers execute no derivative work.
    ///
    /// Each None independently declines a classical derivative at a held normal
    /// or loading/unloading kink, or an unrepresentable/degenerate crossing. A
    /// declined derivative never changes the original physical trial or pressure.
    pub fn trial_over_time_with_damage_weighted_compression_tangent(
        &self,
        jump_m: [f64; 3],
        dt_s: f64,
        energy_tolerance_j_m2: f64,
    ) -> Result<
        (
            CohesiveRateStep<'a>,
            [f64; 3],
            f64,
            Option<[[f64; 3]; 3]>,
            Option<[f64; 3]>,
        ),
        CohesiveError,
    > {
        let (step, mean, tangent) =
            self.trial_over_time_with_mean_traction_tangent(jump_m, dt_s, energy_tolerance_j_m2)?;
        let pressure = damaged_compression(self, &step.equilibrium.point)?;
        let gradient = damaged_compression_gradient(self, &step.equilibrium.point);
        Ok((step, mean, pressure, tangent, gradient))
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

// Only derivative arithmetic is optional. The authoritative trial and pressure
// have passed their original domain, dissipation and nonnegative-measure gates.
fn damaged_compression_gradient(
    old: &CohesivePoint<'_>,
    end: &CohesivePoint<'_>,
) -> Option<[f64; 3]> {
    if old.jump[0] == 0.0 && end.jump[0] == 0.0 {
        // An entire path on the contact plane has distinct compression/opening
        // first variations. Do not advertise either as a classical derivative.
        return None;
    }
    let a = old.jump.map(D::constant);
    let b = core::array::from_fn::<_, 3, _>(|i| D::variable(end.jump[i], i));
    let compression = positive_mean(-a[0], -b[0]);
    if compression.re == 0.0 {
        return Some([0.0; 3]);
    }
    let stiffness = D::constant(old.law.compression_stiffness);
    let damage = old.damage();
    if damage == 1.0 {
        return finite_gradient(stiffness * compression);
    }
    let (weights, onset, failure, lower) = if let Some(p) = old.law.mixed_mode {
        let Some(history) = end.mixed_history else {
            // Below first initiation no damaged area exists anywhere on the
            // original straight path, even when the bond supports compression.
            return Some([0.0; 3]);
        };
        let weights = [
            det::sqrt(p.normal_stiffness_pa_m),
            det::sqrt(p.tangential_stiffness_pa_m),
            det::sqrt(p.tangential_stiffness_pa_m),
        ];
        if old.mixed_history.is_some() {
            (
                weights,
                D::constant(history.onset_sqrt_energy),
                D::constant(history.failure_sqrt_energy),
                0.0,
            )
        } else {
            if !p.smooth_initiation(history.onset_direction) { return None; }
            let (_, lower) = super::super::mixed_mode::initiate(p, old.jump, end.jump).ok()?;
            let (da, dgc) = super::super::mixed_mode::initiation_gradient(
                p, old.jump, end.jump, lower, history,
            )
            .ok()?;
            let onset = D {
                re: history.onset_sqrt_energy,
                eps: da,
            };
            let fracture = D {
                re: history.fracture_energy_j_m2,
                eps: dgc,
            };
            let mut failure = D::constant(2.0) * (fracture / onset);
            failure.re = history.failure_sqrt_energy;
            (weights, onset, failure, lower)
        }
    } else {
        (
            [1.0; 3],
            D::constant(old.law.onset),
            D::constant(old.law.failure),
            0.0,
        )
    };
    let threshold = if old.maximum > onset.re {
        D::constant(old.maximum)
    } else {
        onset
    };
    let weighted_old = core::array::from_fn(|i| weights[i] * old.jump[i]);
    let weighted_end = core::array::from_fn(|i| weights[i] * end.jump[i]);
    let end_radius = norm(active(weighted_end)).ok()?;
    if old.jump == end.jump && end_radius >= threshold.re {
        return None;
    }
    if end_radius <= threshold.re {
        return finite_gradient(stiffness * (D::constant(damage) * compression));
    }
    let start = super::tangent::crossing(weighted_old, weighted_end, weights, threshold, lower)?;
    let finish = if end_radius >= failure.re {
        super::tangent::crossing(weighted_old, weighted_end, weights, failure, start.re)?
    } else {
        D::one()
    };
    let at = |s: D| core::array::from_fn::<_, 3, _>(|i| (D::one() - s) * a[i] + s * b[i]);
    let mut mean = start * D::constant(damage) * positive_mean(-a[0], -at(start)[0]);
    let switch = normal_switch(a[0], b[0]);
    let cut = if switch.re > start.re && switch.re < finish.re {
        switch
    } else {
        finish
    };
    for (left, right) in [(start, cut), (cut, finish)] {
        if right.re > left.re {
            let mut u = at(left);
            let mut v = at(right);
            if 0.5 * u[0].re + 0.5 * v[0].re < 0.0 {
                // Use the INTERIOR compression branch at both ends. A moving
                // damage/failure crossing can coincide with g_n=0; clamping its
                // full dual value would incorrectly erase the endpoint motion.
                let mut p = [-u[0], -v[0]];
                p[0].re = p[0].re.max(0.0);
                p[1].re = p[1].re.max(0.0);
                u[0] = D::zero();
                v[0] = D::zero();
                let u = core::array::from_fn(|i| D::constant(weights[i]) * u[i]);
                let v = core::array::from_fn(|i| D::constant(weights[i]) * v[i]);
                mean = mean + (right - left) * softening_pressure(u, v, p, onset, failure).ok()?;
            }
        }
    }
    mean = mean + (D::one() - finish) * positive_mean(-at(finish)[0], -b[0]);
    finite_gradient(stiffness * mean)
}

fn finite_gradient(value: D) -> Option<[f64; 3]> {
    (value.re.is_finite() && value.eps.iter().all(|entry| entry.is_finite())).then_some(value.eps)
}

// Mean d(r(t))*p(t), with affine nonnegative compression p and weighted active
// jump q on one outward softening segment. Damage is f/(f-o)*(1-o/|q|).
// Its f64 and forward-dual instances retain identical primal arithmetic.
fn softening_pressure<S: Real>(
    a: [S; 3],
    b: [S; 3],
    p: [S; 2],
    onset: S,
    failure: S,
) -> Result<S, CohesiveError> {
    let r0 = scalar_norm(a)?;
    positive(r0.value())?;
    let u = a.map(|x| x / r0);
    let d = core::array::from_fn(|i| (b[i] - a[i]) / r0);
    let ud = scalar_dot(u, d);
    let dd = scalar_dot(d, d);
    let factor = failure / (failure - onset);
    let average = S::from_f64(0.5) * p[0] + S::from_f64(0.5) * p[1];
    let change = p[1] - p[0];
    let result = if 2.0 * ud.value().abs() + dd.value() <= 0.001 {
        // Original inverse-radius moments, kept as small differences so a
        // newly damaged area is not obtained by subtracting two pressures.
        let mut moment = [S::zero(); 2];
        let mut binomial = 1.0;
        for k in 1..=8 {
            binomial *= (-0.5 - f64::from(k - 1)) / f64::from(k);
            let mut choose = 1.0;
            for j in 0..=k {
                if j > 0 {
                    choose *= f64::from(k - j + 1) / f64::from(j);
                }
                let coefficient =
                    S::from_f64(choose) * powi(S::from_f64(2.0) * ud, k - j) * powi(dd, j);
                moment[0] = moment[0]
                    + S::from_f64(binomial) * coefficient / S::from_f64(f64::from(k + j + 1));
                moment[1] = moment[1]
                    + S::from_f64(binomial) * coefficient / S::from_f64(f64::from(k + j + 2));
            }
        }
        (factor * ((r0 - onset) / r0)) * average
            - (factor * (onset / r0)) * (p[0] * moment[0] + change * moment[1])
    } else {
        let length = scalar_norm(d)?;
        let direction = d.map(|x| x / length);
        let projection = scalar_dot(u, direction);
        let end_radius = scalar_norm(core::array::from_fn(|i| u[i] + d[i]))?;
        let radial_rate = (length + S::from_f64(2.0) * projection) / (end_radius + S::one());
        let logarithm = if length.value() >= 1.0 {
            length.ln()
                + (S::one() / length + (S::one() + radial_rate) / (S::one() + projection)).ln()
        } else {
            (S::one() + length * (S::one() + radial_rate) / (S::one() + projection)).ln()
        };
        let j0 = logarithm / length;
        let j1 = (radial_rate - projection * j0) / length;
        factor * (average - (onset / r0) * (p[0] * j0 + change * j1))
    };
    finite(result.value())?;
    Ok(result)
}
