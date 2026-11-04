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
mod tangent;

use fs_ad::Real;

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
        self.trial_over_time_mean(jump_m, dt_s, energy_tolerance_j_m2)
            .map(|(step, mean, _)| (step, mean))
    }

    // Keep the original rate owner's mean intact fraction for the optional
    // path tangent, without replaying its physical trial or work calculation.
    fn trial_over_time_mean(
        &self,
        jump_m: [f64; 3],
        dt_s: f64,
        energy_tolerance_j_m2: f64,
    ) -> Result<(CohesiveRateStep<'a>, [f64; 3], f64), CohesiveError> {
        let (step, mean_intact) =
            self.trial_over_time_integrated(jump_m, dt_s, energy_tolerance_j_m2)?;
        let mut mean = mean_equilibrium(self, &step.equilibrium.point)?;
        if let Some(viscosity) = self.law.viscosity.filter(|_| mean_intact > 0.0) {
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
        Ok((step, mean, mean_intact))
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

fn checked<S: Real>(values: [S; 3]) -> Result<[S; 3], CohesiveError> {
    for value in values {
        finite(value.value())?;
    }
    Ok(values)
}

// These fixed-size primitives are shared by f64 path forces and their optional
// forward derivatives. Primal arithmetic keeps its original ordering; no
// derivative work is done by ordinary mean-traction callers.
fn scalar_dot<S: Real>(a: [S; 3], b: [S; 3]) -> S {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn scalar_norm<S: Real>(value: [S; 3]) -> Result<S, CohesiveError> {
    checked(value)?;
    let scale = value.iter().fold(0.0_f64, |a, b| a.max(b.value().abs()));
    if scale == 0.0 {
        // The zero norm is used only for branch selection in the short-segment
        // series. Its undefined derivative must never enter a nonzero formula.
        return Ok(S::zero());
    }
    let scale = S::from_f64(scale);
    let v = value.map(|x| x / scale);
    let result = scale * scalar_dot(v, v).sqrt();
    finite(result.value())?;
    Ok(result)
}

// Fraction of a sign-changing affine normal jump. Scaling avoids overflow in
// |a|+|b|. A same-sign segment has no internal tension/compression switch.
fn normal_switch<S: Real>(a: S, b: S) -> S {
    if (a.value() < 0.0 && b.value() > 0.0) || (a.value() > 0.0 && b.value() < 0.0) {
        let scale = S::from_f64(a.value().abs().max(b.value().abs()));
        let x = a.abs() / scale;
        x / (x + b.abs() / scale)
    } else {
        S::one()
    }
}

fn positive_mean<S: Real>(a: S, b: S) -> S {
    let half = S::from_f64(0.5);
    if a.value() >= 0.0 && b.value() >= 0.0 {
        half * a + half * b
    } else if a.value() <= 0.0 && b.value() <= 0.0 {
        S::zero()
    } else if a.value() > 0.0 {
        half * a * normal_switch(a, b)
    } else {
        half * b * normal_switch(b, a)
    }
}

fn negative_mean<S: Real>(a: S, b: S) -> S {
    -positive_mean(-a, -b)
}

fn add_elastic<S: Real>(
    mean: &mut [S; 3],
    stiffness: [f64; 3],
    intact: f64,
    duration: S,
    a: [S; 3],
    b: [S; 3],
) {
    let half = S::from_f64(0.5);
    let jump = [
        positive_mean(a[0], b[0]),
        half * a[1] + half * b[1],
        half * a[2] + half * b[2],
    ];
    for i in 0..3 {
        mean[i] = mean[i] + duration * (S::from_f64(intact * stiffness[i]) * jump[i]);
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
    softening_integrals::<f64, false>(a, b, onset, failure).map(|(mean, _)| mean)
}

// The derivative also needs the mean intact fraction for physical viscosity.
// Computing it uses the SAME inverse-radius moments as vector traction. The
// const parameter omits this output in the ordinary f64 vector-force path.
fn softening_integrals<S: Real, const WITH_INTEGRITY: bool>(
    a: [S; 3],
    b: [S; 3],
    onset: f64,
    failure: f64,
) -> Result<([S; 3], S), CohesiveError> {
    let r0 = scalar_norm(a)?;
    positive(r0.value())?;
    let u = a.map(|x| x / r0);
    let delta = core::array::from_fn::<_, 3, _>(|i| b[i] - a[i]);
    let d = delta.map(|x| x / r0);
    let length = scalar_norm(d)?;
    let ud = scalar_dot(u, d);
    let dd = scalar_dot(d, d);
    let half = S::from_f64(0.5);
    let average = core::array::from_fn::<_, 3, _>(|i| half * a[i] + half * b[i]);
    let (result, integrity) = if 2.0 * ud.value().abs() + dd.value() <= 0.001 {
        // Moments of (1 + 2 u.d t + d.d t^2)^(-1/2) - 1. Eight
        // binomial terms have < 2e-28 absolute truncation in this range.
        // Integrating both 1 and t retains the transverse traction. Keeping
        // the small difference explicitly avoids loss of traction near failure.
        let mut moment = [S::zero(); 2];
        let mut binomial = 1.0;
        for k in 1..=8 {
            binomial *= (-0.5 - f64::from(k - 1)) / f64::from(k);
            let mut choose = 1.0;
            for j in 0..=k {
                if j > 0 {
                    choose *= f64::from(k - j + 1) / f64::from(j);
                }
                let coefficient = S::from_f64(choose)
                    * powi(S::from_f64(2.0) * ud, k - j)
                    * powi(dd, j);
                moment[0] = moment[0]
                    + S::from_f64(binomial) * coefficient / S::from_f64(f64::from(k + j + 1));
                moment[1] = moment[1]
                    + S::from_f64(binomial) * coefficient / S::from_f64(f64::from(k + j + 2));
            }
        }
        let base = (S::from_f64(onset) / r0)
            * ((S::from_f64(failure) - r0) / S::from_f64(failure - onset));
        let correction = (S::from_f64(onset) / r0) * S::from_f64(failure / (failure - onset));
        (
            core::array::from_fn(|i| {
                base * average[i] + correction * (a[i] * moment[0] + delta[i] * moment[1])
            }),
            if WITH_INTEGRITY {
                base + correction * moment[0]
            } else {
                S::zero()
            },
        )
    } else {
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
        let inverse_radius_mean = logarithm / length;
        let factor = S::from_f64(onset / (failure - onset));
        (
            core::array::from_fn(|i| {
                let unit_mean = direction[i] * radial_rate
                    + (u[i] - projection * direction[i]) * inverse_radius_mean;
                factor * (S::from_f64(failure) * unit_mean - average[i])
            }),
            if WITH_INTEGRITY {
                factor * ((S::from_f64(failure) / r0) * inverse_radius_mean - S::one())
            } else {
                S::zero()
            },
        )
    };
    finite(integrity.value())?;
    Ok((checked(result)?, integrity))
}

fn powi<S: Real>(x: S, n: i32) -> S {
    let mut result = S::one();
    for _ in 0..n {
        result = result * x;
    }
    result
}

#[cfg(test)]
mod tests;
