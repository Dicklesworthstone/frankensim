//! Derivative of the original straight-path mean traction at fixed old history.
//!
//! Differentiate the shared analytic primitives, including motion of the damage,
//! failure and normal-contact crossings. The endpoint tangent is a different
//! quantity. The old maximum, frozen BK envelope and all dissipation references
//! are constants in this derivative; no trial can advance accepted history.

use super::*;
use fs_ad::Dual64;

type D = Dual64<3>;

impl<'a> CohesivePoint<'a> {
    /// Original equilibrium trial and path mean, plus
    /// d(mean traction)/d(endpoint jump) in Pa/m at fixed accepted history.
    ///
    /// Common-mode damage and an already selected frozen-BK envelope are
    /// differentiated through unloading, softening, rupture and normal sign
    /// changes. This is the derivative of the vector path integral, NOT the
    /// returned step's endpoint tangent. No quadrature or difference scale is
    /// introduced; arithmetic and scratch storage are bounded and allocation-free.
    ///
    /// None explicitly requests a consumer's existing numerical/semismooth
    /// fallback: first mixed-mode initiation is not covered, nor a path held
    /// identically at a normal or damage loading/unloading kink, or an
    /// unrepresentable/degenerate crossing derivative. The original trial/mean
    /// still obey every material-domain and energy gate. No selected one-sided
    /// branch is mislabeled as a classical derivative.
    pub fn trial_with_mean_traction_tangent(
        &self,
        jump_m: [f64; 3],
        energy_tolerance_j_m2: f64,
    ) -> Result<(CohesiveStep<'a>, [f64; 3], Option<[[f64; 3]; 3]>), CohesiveError> {
        let (step, mean) = self.trial_with_mean_traction(jump_m, energy_tolerance_j_m2)?;
        let tangent = mean_tangent(self, &step.point, None);
        Ok((step, mean, tangent))
    }

    /// Timed path-mean derivative at fixed accepted history and fixed positive
    /// dt. Includes BOTH eta_i * mean(1-d) / dt and the sensitivity of the
    /// evolving intact-fraction integral multiplying the actual jump rate.
    /// Physical viscous loss and intrinsic fracture remain separate, unchanged.
    ///
    /// The same explicit None cases as trial_with_mean_traction_tangent apply.
    /// A newly rupturing path still has its pre-rupture dashpot derivative; a
    /// path starting failed has no surviving bond dashpots.
    pub fn trial_over_time_with_mean_traction_tangent(
        &self,
        jump_m: [f64; 3],
        dt_s: f64,
        energy_tolerance_j_m2: f64,
    ) -> Result<(CohesiveRateStep<'a>, [f64; 3], Option<[[f64; 3]; 3]>), CohesiveError> {
        let (step, mean, intact) =
            self.trial_over_time_mean(jump_m, dt_s, energy_tolerance_j_m2)?;
        let tangent = mean_tangent(self, &step.equilibrium.point, Some((dt_s, intact)));
        Ok((step, mean, tangent))
    }
}

// The authoritative trial and its path work have already succeeded. Optional
// derivative arithmetic may decline support; it cannot alter the physical
// response or turn a valid force into a fabricated zero Jacobian.
fn mean_tangent(
    old: &CohesivePoint<'_>,
    end: &CohesivePoint<'_>,
    timed: Option<(f64, f64)>,
) -> Option<[[f64; 3]; 3]> {
    if old.jump[0] == 0.0 && end.jump[0] == 0.0 {
        return None;
    }
    if old.law.mixed_mode.is_some()
        && old.mixed_history.is_none()
        && end.mixed_history.is_some()
    {
        // The first onset mode itself depends on the endpoint. Freezing it here
        // would silently omit initiation/toughness sensitivity.
        return None;
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
    let a = old.jump.map(D::constant);
    let b = core::array::from_fn(|i| D::variable(end.jump[i], i));
    let mut mean = [
        D::constant(old.law.compression_stiffness) * negative_mean(a[0], b[0]),
        D::zero(),
        D::zero(),
    ];
    let mut mean_intact = D::constant(intact);
    if intact > 0.0 {
        let envelope = if old.law.mixed_mode.is_some() {
            old.mixed_history.map(|h| {
                (
                    stiffness.map(det::sqrt),
                    h.onset_sqrt_energy,
                    h.failure_sqrt_energy,
                )
            })
        } else {
            Some(([1.0; 3], old.law.onset, old.law.failure))
        };
        if let Some((weights, onset, failure)) = envelope {
            let threshold = old.maximum.max(onset);
            let weighted_old = core::array::from_fn(|i| weights[i] * old.jump[i]);
            let weighted_end = core::array::from_fn(|i| weights[i] * end.jump[i]);
            let end_radius = norm(active(weighted_end)).ok()?;
            if old.jump == end.jump && end_radius >= threshold {
                // At a held maximum, loading and unloading have different
                // first variations even though their force values coincide.
                return None;
            }
            if end_radius <= threshold {
                add_elastic(&mut mean, stiffness, intact, D::one(), a, b);
            } else {
                let start = crossing(weighted_old, weighted_end, weights, threshold, 0.0)?;
                let finish = if end_radius >= failure {
                    crossing(weighted_old, weighted_end, weights, failure, start.re)?
                } else {
                    D::one()
                };
                let at = |s: D| {
                    core::array::from_fn(|i| (D::one() - s) * a[i] + s * b[i])
                };
                add_elastic(&mut mean, stiffness, intact, start, a, at(start));
                mean_intact = D::constant(intact) * start;
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
                        // A segment lies on one side of the normal switch.
                        // Select its INTERIOR branch for both endpoints: when
                        // a damage/failure crossing coincides with g_n=0, the
                        // endpoint's own derivative need not stay on that plane.
                        if 0.5 * u[0].re + 0.5 * v[0].re <= 0.0 {
                            u[0] = D::zero();
                            v[0] = D::zero();
                        } else {
                            u[0].re = u[0].re.max(0.0);
                            v[0].re = v[0].re.max(0.0);
                        }
                        let u = core::array::from_fn(|i| D::constant(weights[i]) * u[i]);
                        let v = core::array::from_fn(|i| D::constant(weights[i]) * v[i]);
                        let (traction, integrity) =
                            softening_integrals::<D, true>(u, v, onset, failure).ok()?;
                        for i in 0..3 {
                            mean[i] = mean[i]
                                + (right - left)
                                    * (D::constant(stiffness[i] / weights[i]) * traction[i]);
                        }
                        mean_intact = mean_intact + (right - left) * integrity;
                    }
                }
            }
        } else {
            // A mixed-mode point that remains below first initiation is elastic.
            add_elastic(&mut mean, stiffness, intact, D::one(), a, b);
        }
    }
    if let (Some(viscosity), Some((dt_s, original_mean_intact))) = (old.law.viscosity, timed) {
        if original_mean_intact > 0.0 {
            // Preserve the original rate integral's primal rounding. Its
            // derivative comes from the shared inverse-radius path primitive.
            mean_intact.re = original_mean_intact;
            let eta = [
                viscosity.normal_pa_s_m(),
                viscosity.tangential_pa_s_m(),
                viscosity.tangential_pa_s_m(),
            ];
            for i in 0..3 {
                let rate = (b[i] - a[i]) / D::constant(dt_s);
                mean[i] = mean[i] + mean_intact * (D::constant(eta[i]) * rate);
            }
        } else if intact > 0.0 {
            // An underflow/clamped integral has no resolved derivative here.
            return None;
        }
    }
    let tangent = mean.map(|entry| entry.eps);
    tangent.iter().flatten().all(|value| value.is_finite()).then_some(tangent)
}

// Reuse the original accurately scaled outward root for its value. Implicit
// differentiation of |q_old+s(q_end-q_old)|=radius gives
// ds/dg_end,j = -s * weight_j * unit_q_j / (unit_q dot delta_q).
// Radius is fixed by accepted history, never by a new endpoint-selected mode.
// Normal compression has unit_q_n=0. Scaling avoids squaring a metre/energy
// radius or differentiating through the root's numerical bracketing/clamps.
fn crossing(
    old: [f64; 3],
    end: [f64; 3],
    weights: [f64; 3],
    radius: f64,
    lower: f64,
) -> Option<D> {
    let fraction = outward_crossing(old, end, radius, lower).ok()?;
    let unit = active(core::array::from_fn(|i| {
        ((1.0 - fraction) * old[i] + fraction * end[i]) / radius
    }));
    let delta = core::array::from_fn(|i| end[i] - old[i]);
    let length = norm(delta).ok()?;
    if length == 0.0 {
        return None;
    }
    let direction = delta.map(|value| value / length);
    let projection = dot(unit, direction);
    if !projection.is_finite() || projection <= 64.0 * f64::EPSILON {
        return None;
    }
    let eps = core::array::from_fn(|i| {
        -(fraction * unit[i] / projection) * (weights[i] / length)
    });
    eps.iter().all(|value| value.is_finite()).then_some(D { re: fraction, eps })
}
