//! First physical failure on a straight shear-stage segment at held D_n.
//!
//! Radius probes eliminate the normal spring analytically; they never rebuild
//! fixtures or requery material cards. The limiting radius is bracketed on the
//! original full-stick branch. The published response is the canonical last
//! admitted stage just inside that boundary, retaining all contact/apparatus
//! storage. No slip, snap, release, heat, or post-event motion is invented.

use super::{
    JkrPullOffBranch, JkrPullOffError, JkrPullOffFixture, JkrPullOffResponse, JkrPullOffState,
    JkrShearStageTransition, Work, finite,
};
use crate::normal_patch::adhesive::mixed_mode::loading::JkrShearBoundaryKind;
use core::f64::consts::PI;
use fs_exec::CancelGate;

/// A located first physical boundary; source/applicability failures remain
/// errors instead of being relabelled as a slip or peel event.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrShearStageEvent {
    pub kind: JkrShearBoundaryKind,
    /// Numerical bracket for the critical radius, not an interval certificate.
    /// The returned contacting response stays on its admitted side; near a
    /// fold its radius can differ by more than the stage rounding guard.
    pub critical_radius_bracket_m: [f64; 2],
    /// Parameter on D_t(s)=(1-s)D_t(0)+sD_t(requested), in [0,1].
    pub accepted_fraction: f64,
}

impl JkrPullOffFixture {
    /// Follow a shear-stage segment until its first sourced slip or radius fold.
    ///
    /// Below that boundary this is the same reversible motion as
    /// `advance_tangential_stage`. When the request crosses it, return the last
    /// admitted stage, the limiting radius bracket, and unused vector travel.
    /// Reversals first unload along the same branch, then reach the outward
    /// crossing of the limiting-magnitude circle; zero net tangential direction
    /// is not assumed. Interior expansion and every original source screen
    /// still apply to the accepted segment.
    ///
    /// A 128-epsilon relative magnitude guard keeps the canonical endpoint on
    /// the admissible side of the numerical boundary. It is not a friction
    /// regularization or a finite physical slip distance. All radius searches,
    /// the single successor construction and the actual equilibrium share one
    /// original evaluation budget. Boundary storage and previous released
    /// energy remain intact for a separately selected event owner.
    pub fn advance_tangential_stage_until_event(
        &self,
        previous: JkrPullOffState,
        requested_stage_m: [f64; 2],
        gate: &CancelGate,
    ) -> Result<JkrShearStageTransition, JkrPullOffError> {
        if previous.identity != self.identity {
            return Err(JkrPullOffError::StateMismatch);
        }
        for value in requested_stage_m {
            finite(value)?;
        }
        let requested = requested_stage_m[0].hypot(requested_stage_m[1]);
        finite(requested)?;
        let current = self.tangential_stage[0].hypot(self.tangential_stage[1]);
        let mut budget = Work::new(self.limits.maximum_evaluations, gate);
        let before = self.shear_stage_before(previous, &mut budget)?;
        if previous.branch == JkrPullOffBranch::Detached || requested <= current {
            return self.move_shear_from_before(before, requested_stage_m, &mut budget);
        }

        let (magnitude, kind, bracket) = self.shear_stage_limit(before, &mut budget)?;
        if requested <= magnitude {
            return self.move_shear_from_before(before, requested_stage_m, &mut budget);
        }
        if magnitude < current * (1.0 - 256.0 * f64::EPSILON) {
            return Err(JkrPullOffError::NoConvergence);
        }
        let retained = (magnitude * (1.0 - 128.0 * f64::EPSILON)).max(current);
        let (fraction, accepted) =
            outward_crossing(self.tangential_stage, requested_stage_m, retained)?;
        if fraction == 0.0 && retained > current {
            return Err(JkrPullOffError::NoConvergence);
        }
        let mut result = self.move_shear_from_before(before, accepted, &mut budget)?;
        result.event = Some(JkrShearStageEvent {
            kind,
            critical_radius_bracket_m: bracket,
            accepted_fraction: fraction,
        });
        result.requested_stage_m = requested_stage_m;
        result.untraversed_stage_m =
            core::array::from_fn(|i| requested_stage_m[i] - result.source.tangential_stage[i]);
        for value in result.untraversed_stage_m {
            finite(value)?;
        }
        budget.poll()?;
        Ok(result)
    }

    fn shear_stage_limit(
        &self,
        before: JkrPullOffResponse,
        budget: &mut Work<'_>,
    ) -> Result<(f64, JkrShearBoundaryKind, [f64; 2]), JkrPullOffError> {
        let a = before
            .contact
            .ok_or(JkrPullOffError::StateMismatch)?
            .contact_radius_m;
        let current = self.tangential_stage[0].hypot(self.tangential_stage[1]);
        if before.static_capacity_n == 0.0 {
            return Ok((
                current,
                if before.stage_stiffness_n_per_m.is_none() {
                    JkrShearBoundaryKind::Coincident
                } else {
                    JkrShearBoundaryKind::GrossSlip
                },
                [a; 2],
            ));
        }
        if before.state.normal_stage_m == self.boundary.response.state.normal_stage_m {
            return Ok((current, self.boundary.kind, self.boundary.radius_bracket_m));
        }
        let normal = self.contact.normal_source();
        let kernel = HeldNormal {
            e: normal.normal_state().elastic().reduced_modulus_pa(),
            w: normal.normal_state().adhesion_energy_j_per_m2(),
            r: normal.effective_radius_m(),
            k: self.contact.tangential_stiffness_per_radius,
            kn: self.normal_stiffness,
            kt: self.tangential_stiffness,
            dn: before.state.normal_stage_m,
            strength: self.shear.law().static_shear_strength_pa(),
        };
        let start = kernel.sample(a, budget)?;
        if start.fold_derivative >= 0.0 {
            return Err(JkrPullOffError::NoConvergence);
        }
        let (fold, fold_bracket) =
            kernel.last_admitted(0.0, start, self.limits.maximum_iterations, false, budget)?;
        // f(a)-k*u_capacity(a)^2/2 is strictly concave. On the retained
        // branch its crossing is unique because |u|/a grows as a decreases.
        let uncertainty = 64.0 * f64::EPSILON * fold.f.abs().max(fold.capacity_energy);
        if fold.slip_excess > uncertainty {
            if start.slip_excess > 0.0 {
                return Err(JkrPullOffError::NoConvergence);
            }
            let (slip, bracket) = kernel.last_admitted(
                fold.radius,
                start,
                self.limits.maximum_iterations,
                true,
                budget,
            )?;
            Ok((
                kernel.magnitude(slip)?,
                if bracket[0] > fold_bracket[1] {
                    JkrShearBoundaryKind::GrossSlip
                } else {
                    JkrShearBoundaryKind::Coincident
                },
                bracket,
            ))
        } else {
            Ok((
                kernel.magnitude(fold)?,
                if fold.slip_excess.abs() <= uncertainty {
                    JkrShearBoundaryKind::Coincident
                } else {
                    JkrShearBoundaryKind::RadiusInstability
                },
                fold_bracket,
            ))
        }
    }
}

struct HeldNormal {
    e: f64,
    w: f64,
    r: f64,
    k: f64,
    kn: f64,
    kt: f64,
    dn: f64,
    strength: f64,
}
#[derive(Clone, Copy)]
struct RadiusTrial {
    radius: f64,
    f: f64,
    fold_derivative: f64,
    capacity_energy: f64,
    slip_excess: f64,
}
impl HeldNormal {
    fn sample(&self, a: f64, budget: &mut Work<'_>) -> Result<RadiusTrial, JkrPullOffError> {
        budget.probe()?;
        let c = 2.0 * self.e / self.kn;
        let q = 1.0 + c * a;
        let hertz = a * (a / self.r);
        let delta = self.dn / q + (c * a / q) * (hertz / 3.0);
        let g = hertz - delta;
        let f = 2.0 * PI * self.w * a - (self.e * g) * g;
        let fp = 2.0 * PI * self.w - 2.0 * self.e * g * (2.0 * a / self.r - c * g / q);
        let raw_derivative = fp + 2.0 * f / (a + self.kt / self.k);
        // Below the positive-opening/positive-f interval, extend only the
        // bracket sign. No contacting response is published from that region.
        let fold_derivative = if g <= 0.0 || (f <= 0.0 && fp > 0.0) {
            1.0
        } else {
            raw_derivative
        };
        let capacity_shift = PI * self.strength * a / self.k;
        let capacity_energy = (0.5 * self.k * capacity_shift) * capacity_shift;
        let slip_excess = f - capacity_energy;
        for value in [
            q,
            delta,
            g,
            f,
            fp,
            raw_derivative,
            capacity_energy,
            slip_excess,
        ] {
            finite(value)?;
        }
        Ok(RadiusTrial {
            radius: a,
            f,
            fold_derivative,
            capacity_energy,
            slip_excess,
        })
    }

    fn magnitude(&self, trial: RadiusTrial) -> Result<f64, JkrPullOffError> {
        if trial.f <= 0.0 {
            return Err(JkrPullOffError::NoConvergence);
        }
        let magnitude = (2.0 * trial.f / self.k).sqrt() * (1.0 + self.k * trial.radius / self.kt);
        finite(magnitude)?;
        Ok(magnitude)
    }

    fn last_admitted(
        &self,
        mut lo: f64,
        mut hi: RadiusTrial,
        maximum_iterations: usize,
        slip: bool,
        budget: &mut Work<'_>,
    ) -> Result<(RadiusTrial, [f64; 2]), JkrPullOffError> {
        for _ in 0..maximum_iterations {
            let mid = f64::midpoint(lo, hi.radius);
            if mid == lo || mid == hi.radius || hi.radius - lo <= 16.0 * f64::EPSILON * hi.radius {
                return Ok((hi, [lo, hi.radius]));
            }
            let trial = self.sample(mid, budget)?;
            let sign = if slip {
                trial.slip_excess
            } else {
                trial.fold_derivative
            };
            if sign > 0.0 {
                lo = mid;
            } else {
                hi = trial;
            }
        }
        Err(JkrPullOffError::NoConvergence)
    }
}

/// Outward intersection with the limiting-magnitude circle. A stable quadratic
/// root handles loading, unloading followed by reversal, and immediate failure.
/// Coordinates are scaled before subtraction/squaring to avoid overflow.
fn outward_crossing(
    from: [f64; 2],
    to: [f64; 2],
    magnitude: f64,
) -> Result<(f64, [f64; 2]), JkrPullOffError> {
    let scale = from.into_iter().chain(to).map(f64::abs).fold(0.0, f64::max);
    let p = from.map(|x| x / scale);
    let q = to.map(|x| x / scale);
    let d = [q[0] - p[0], q[1] - p[1]];
    let a = d[0] * d[0] + d[1] * d[1];
    let b = p[0] * d[0] + p[1] * d[1];
    let initial = p[0].hypot(p[1]);
    let limit = magnitude / scale;
    let c = ((initial - limit) * (initial + limit)).min(0.0);
    let root = (b * b - a * c).sqrt();
    let t = if b >= 0.0 {
        if b + root == 0.0 {
            0.0
        } else {
            -c / (b + root)
        }
    } else {
        (-b + root) / a
    };
    finite(t)?;
    if !(0.0..=1.0 + 16.0 * f64::EPSILON).contains(&t) {
        return Err(JkrPullOffError::NoConvergence);
    }
    let t = t.min(1.0);
    let stage = core::array::from_fn(|i| ((1.0 - t) * p[i] + t * q[i]) * scale);
    for value in stage {
        finite(value)?;
    }
    Ok((t, stage))
}
