//! Analytic vector traction integrals on the constitutive owner's straight path.
//!
//! Endpoint traction and its tangent remain endpoint quantities. These averages
//! supply an actual work-conjugate force to a dynamic owner, including shear and
//! changing direction; a scalar work divided by a normal displacement cannot.

use super::{
    CohesiveError, CohesivePoint, CohesiveRateStep, CohesiveStep, active, det, dot, finite, norm,
    outward_crossing, positive,
};

mod friction;

impl<'a> CohesivePoint<'a> {
    /// Return the original endpoint trial and the mean equilibrium traction
    /// `integral_0^1 t(g_old + s * (g_new-g_old)) ds`, in Pa. Its dot product
    /// with the complete three-component jump change equals the original path
    /// work, within the supplied absolute J/m2 gate. No damage is committed.
    ///
    /// The mean is an analytic path integral, not endpoint or midpoint traction.
    /// The returned step's tangents still differentiate ENDPOINT traction, not
    /// this mean. Compression, retained damage, first frozen-BK initiation, and
    /// complete failure use the original constitutive history. Allocation-free,
    /// bounded work; physical viscosity requires the explicit timed analogue.
    pub fn trial_with_mean_traction(
        &self,
        jump_m: [f64; 3],
        energy_tolerance_j_m2: f64,
    ) -> Result<(CohesiveStep<'a>, [f64; 3]), CohesiveError> {
        let step = self.trial(jump_m, energy_tolerance_j_m2)?;
        let mean = mean_equilibrium(self, &step.point)?;
        check_work(
            self.jump,
            jump_m,
            mean,
            step.energy.work,
            energy_tolerance_j_m2,
        )?;
        Ok((step, mean))
    }

    /// Timed counterpart of `trial_with_mean_traction`. The mean includes the
    /// surviving bonds' physical dashpot traction at the declared constant jump
    /// rate. Its work equals intrinsic work plus the separately reported viscous
    /// dissipation. Endpoint traction/tangents retain their original meaning.
    pub fn trial_over_time_with_mean_traction(
        &self,
        jump_m: [f64; 3],
        dt_s: f64,
        energy_tolerance_j_m2: f64,
    ) -> Result<(CohesiveRateStep<'a>, [f64; 3]), CohesiveError> {
        let (step, mean_intact) =
            self.trial_over_time_integrated(jump_m, dt_s, energy_tolerance_j_m2)?;
        let mut mean = mean_equilibrium(self, &step.equilibrium.point)?;
        if let Some(viscosity) = self.law.viscosity {
            let eta = [
                viscosity.normal_pa_s_m(),
                viscosity.tangential_pa_s_m(),
                viscosity.tangential_pa_s_m(),
            ];
            for i in 0..3 {
                let rate = (jump_m[i] - self.jump[i]) / dt_s;
                mean[i] = finite(mean[i] + mean_intact * (eta[i] * rate))?;
            }
        }
        check_work(
            self.jump,
            jump_m,
            mean,
            step.total_work_j_m2(),
            energy_tolerance_j_m2,
        )?;
        Ok((step, mean))
    }
}

fn check_work(
    old: [f64; 3],
    end: [f64; 3],
    mean: [f64; 3],
    work: f64,
    tolerance: f64,
) -> Result<(), CohesiveError> {
    let delta = core::array::from_fn(|i| end[i] - old[i]);
    let defect = finite(dot(mean, delta) - work)?;
    if defect.abs() > tolerance {
        return Err(CohesiveError::EnergyDefect(defect));
    }
    Ok(())
}

fn mean_equilibrium(
    old: &CohesivePoint<'_>,
    end: &CohesivePoint<'_>,
) -> Result<[f64; 3], CohesiveError> {
    if old.jump == end.jump {
        return Ok(end.traction_pa());
    }
    let stiffness = old.law.mixed_mode.map_or([old.law.stiffness; 3], |p| {
        [
            p.normal_stiffness_pa_m,
            p.tangential_stiffness_pa_m,
            p.tangential_stiffness_pa_m,
        ]
    });
    let intact = if old.law.mixed_mode.is_some() {
        old.mixed_integrity()
    } else {
        old.law.intact(old.maximum)
    };
    let mut mean = [
        old.law.compression_stiffness * negative_mean(old.jump[0], end.jump[0]),
        0.0,
        0.0,
    ];
    if intact == 0.0 {
        return checked(mean);
    }
    let (weights, onset, failure, threshold, lower) = match old.law.mixed_mode {
        Some(p) => {
            let Some(h) = end.mixed_history else {
                add_elastic(&mut mean, stiffness, intact, 1.0, old.jump, end.jump);
                return checked(mean);
            };
            let lower = if old.mixed_history.is_none() {
                super::mixed_mode::initiate(p, old.jump, end.jump)?.1
            } else {
                0.0
            };
            (
                stiffness.map(det::sqrt),
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
        add_elastic(&mut mean, stiffness, intact, 1.0, old.jump, end.jump);
        return checked(mean);
    }
    let start = outward_crossing(weighted_old, weighted_end, threshold, lower)?;
    let finish = if end_radius >= failure {
        outward_crossing(weighted_old, weighted_end, failure, start)?
    } else {
        1.0
    };
    let at = |s: f64| core::array::from_fn(|i| (1.0 - s) * old.jump[i] + s * end.jump[i]);
    add_elastic(&mut mean, stiffness, intact, start, old.jump, at(start));
    let switch = normal_switch(old.jump[0], end.jump[0]);
    let cut = if switch > start && switch < finish {
        switch
    } else {
        finish
    };
    for (a, b) in [(start, cut), (cut, finish)] {
        if b > a {
            let u = core::array::from_fn(|i| weights[i] * active(at(a))[i]);
            let v = core::array::from_fn(|i| weights[i] * active(at(b))[i]);
            let integral = softening_mean(u, v, onset, failure)?;
            for i in 0..3 {
                mean[i] += (b - a) * ((stiffness[i] / weights[i]) * integral[i]);
            }
        }
    }
    checked(mean)
}

fn checked(values: [f64; 3]) -> Result<[f64; 3], CohesiveError> {
    for value in values {
        finite(value)?;
    }
    Ok(values)
}

// Fraction of a sign-changing affine normal jump. Scaling avoids overflow in
// |a|+|b|. A same-sign segment has no internal tension/compression switch.
fn normal_switch(a: f64, b: f64) -> f64 {
    if (a < 0.0 && b > 0.0) || (a > 0.0 && b < 0.0) {
        let scale = a.abs().max(b.abs());
        let x = a.abs() / scale;
        x / (x + b.abs() / scale)
    } else {
        1.0
    }
}

fn positive_mean(a: f64, b: f64) -> f64 {
    if a >= 0.0 && b >= 0.0 {
        0.5 * a + 0.5 * b
    } else if a <= 0.0 && b <= 0.0 {
        0.0
    } else if a > 0.0 {
        0.5 * a * normal_switch(a, b)
    } else {
        0.5 * b * normal_switch(b, a)
    }
}

fn negative_mean(a: f64, b: f64) -> f64 {
    -positive_mean(-a, -b)
}

fn add_elastic(
    mean: &mut [f64; 3],
    stiffness: [f64; 3],
    intact: f64,
    duration: f64,
    a: [f64; 3],
    b: [f64; 3],
) {
    let jump = [
        positive_mean(a[0], b[0]),
        0.5 * a[1] + 0.5 * b[1],
        0.5 * a[2] + 0.5 * b[2],
    ];
    for i in 0..3 {
        mean[i] += duration * ((intact * stiffness[i]) * jump[i]);
    }
}

// Average (1-d) q for affine weighted jump q, on a monotonically outward,
// entirely softening interval with no normal sign switch. The mean traction is
// K_i/weight_i times this result. No quadrature or vector work projection is used.
fn softening_mean(
    a: [f64; 3],
    b: [f64; 3],
    onset: f64,
    failure: f64,
) -> Result<[f64; 3], CohesiveError> {
    let r0 = positive(norm(a)?)?;
    let u = a.map(|x| x / r0);
    let delta = core::array::from_fn::<_, 3, _>(|i| b[i] - a[i]);
    let d = delta.map(|x| x / r0);
    let length = norm(d)?;
    let ud = dot(u, d);
    let dd = dot(d, d);
    let average = core::array::from_fn::<_, 3, _>(|i| 0.5 * a[i] + 0.5 * b[i]);
    let result = if 2.0 * ud.abs() + dd <= 0.001 {
        // Moments of (1 + 2 u.d t + d.d t^2)^(-1/2) - 1. Eight
        // binomial terms have < 2e-28 absolute truncation in this range.
        // Integrating both 1 and t retains the transverse traction. Keeping
        // the small difference explicitly avoids loss of traction near failure.
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
        let base = (onset / r0) * ((failure - r0) / (failure - onset));
        let correction = (onset / r0) * (failure / (failure - onset));
        core::array::from_fn(|i| {
            base * average[i] + correction * (a[i] * moment[0] + delta[i] * moment[1])
        })
    } else {
        let direction = d.map(|x| x / length);
        let projection = dot(u, direction);
        let end_radius = norm(core::array::from_fn(|i| u[i] + d[i]))?;
        let radial_rate = (length + 2.0 * projection) / (end_radius + 1.0);
        let logarithm = if length >= 1.0 {
            det::ln(length) + det::ln(1.0 / length + (1.0 + radial_rate) / (1.0 + projection))
        } else {
            det::ln(1.0 + length * (1.0 + radial_rate) / (1.0 + projection))
        };
        let inverse_radius_mean = logarithm / length;
        let factor = onset / (failure - onset);
        core::array::from_fn(|i| {
            let unit_mean = direction[i] * radial_rate
                + (u[i] - projection * direction[i]) * inverse_radius_mean;
            factor * (failure * unit_mean - average[i])
        })
    };
    checked(result)
}

fn powi(x: f64, n: i32) -> f64 {
    let mut result = 1.0;
    for _ in 0..n {
        result *= x;
    }
    result
}

#[cfg(test)]
mod tests;
