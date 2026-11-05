//! Combined peel/shear loading of the original finite apparatus.
//!
//! On one straight stage path, normal withdrawal and nondecreasing shear-stage
//! magnitude make the retained contact radius decrease. The tips are eliminated
//! analytically; a quadratic supplies path fraction at each radius. One scalar
//! fold and, when earlier, one static-capacity root locate the first event.
//! No fixture reconstruction or source query occurs inside those searches.

use super::events::JkrShearStageEvent;
use super::{
    JkrPullOffBranch, JkrPullOffError, JkrPullOffFixture, JkrPullOffResponse, JkrPullOffState,
    Work, finite,
};
use crate::normal_patch::adhesive::mixed_mode::loading::JkrShearBoundaryKind;
use core::f64::consts::PI;
use fs_blake3::ContentHash;
use fs_exec::CancelGate;

/// One accepted prefix of a straight [normal, shear-1, shear-2] stage path.
/// Accept the source and response together. Original sources, spring constants,
/// total storage and previously released energy remain with this same owner.
#[derive(Clone, Debug, PartialEq)]
pub struct JkrCombinedStageTransition {
    original_fixture_identity: ContentHash,
    source: JkrPullOffFixture,
    before: JkrPullOffResponse,
    response: JkrPullOffResponse,
    event: Option<JkrShearStageEvent>,
    requested_stages_m: [f64; 3],
    untraversed_stages_m: [f64; 3],
    mechanical_work_j: f64,
    evaluations: usize,
}
impl JkrCombinedStageTransition {
    pub const fn original_fixture_identity(&self) -> ContentHash {
        self.original_fixture_identity
    }
    pub const fn source(&self) -> &JkrPullOffFixture {
        &self.source
    }
    pub fn into_source(self) -> JkrPullOffFixture {
        self.source
    }
    pub const fn before(&self) -> &JkrPullOffResponse {
        &self.before
    }
    pub const fn response(&self) -> &JkrPullOffResponse {
        &self.response
    }
    pub const fn event(&self) -> Option<JkrShearStageEvent> {
        self.event
    }
    pub const fn requested_stages_m(&self) -> [f64; 3] {
        self.requested_stages_m
    }
    pub const fn untraversed_stages_m(&self) -> [f64; 3] {
        self.untraversed_stages_m
    }
    /// Integral N dD_n + T_1 dD_1 + T_2 dD_2, equal to the complete
    /// reversible potential change. It is not a shear-only work receipt.
    pub const fn mechanical_work_j(&self) -> f64 {
        self.mechanical_work_j
    }
    pub const fn evaluations(&self) -> usize {
        self.evaluations
    }
}

impl JkrPullOffFixture {
    /// Move all three apparatus stages together until their first slip/fold.
    ///
    /// An attached path requires D_n(1)<=D_n(0) and
    /// D_t(0).(D_t(1)-D_t(0))>=0: withdrawal with a shear magnitude that never
    /// decreases on this straight segment. Nonmonotone combined motion refuses;
    /// held-normal reversals have their separate existing continuation owner.
    /// Detached motion retains its history and accepts any nonoverlapping path.
    ///
    /// Write c=2E*/K_n, g=(a²/R+c*2a³/(3R)-D_n)/(1+c*a),
    /// u=K_t D_t/(K_t+k*a), k=8G*. Equilibrium is
    /// F=E*g²+k|u|²/2-2*pi*W*a=0. On g>0, g_aa>0 and F_aa>0.
    /// Along the admitted path F_s>=0, so the stable radius decreases until its
    /// unique fold. The static-capacity ratio and both shear screens increase;
    /// indentation decreases. Radius screens are monotone and the original
    /// nominal-strain screen is convex in radius, so endpoint maxima suffice.
    /// This proves the complete accepted prefix, without endpoint-only admission
    /// of an arbitrary multiaxis loading path.
    ///
    /// The event retains the canonical equilibrium just inside the boundary,
    /// its numerical critical-radius bracket and all three unused travels.
    /// A 128-epsilon stage-coordinate guard is a numerical allowance, not a
    /// physical slip length or an interval certificate. One original work budget
    /// covers all probes and the single successor construction. No snap, slip,
    /// new release, heat partition or postevent motion is inferred.
    pub fn withdraw_and_shear_until_event(
        &self,
        previous: JkrPullOffState,
        requested_stages_m: [f64; 3],
        gate: &CancelGate,
    ) -> Result<JkrCombinedStageTransition, JkrPullOffError> {
        if previous.identity != self.identity {
            return Err(JkrPullOffError::StateMismatch);
        }
        self.check_stage(requested_stages_m[0])?;
        for value in requested_stages_m {
            finite(value)?;
        }
        let from = [
            previous.normal_stage_m,
            self.tangential_stage[0],
            self.tangential_stage[1],
        ];
        let mut budget = Work::new(self.limits.maximum_evaluations, gate);
        let before = self.shear_stage_before(previous, &mut budget)?;
        let mut accepted = requested_stages_m;
        let mut event = None;
        if previous.branch == JkrPullOffBranch::Attached && from != requested_stages_m {
            let path = StagePath::new(self, from, requested_stages_m)?;
            let (limit, kind, bracket) = self.combined_stage_limit(before, &path, &mut budget)?;
            if limit.fraction < 1.0 {
                let boundary = interpolate(from, requested_stages_m, limit.fraction);
                let scale = norm(from).max(norm(boundary));
                let guard = 128.0 * f64::EPSILON * (scale / norm(path.increment));
                let fraction = (limit.fraction - guard).max(0.0);
                accepted = interpolate(from, requested_stages_m, fraction);
                event = Some(JkrShearStageEvent {
                    kind,
                    critical_radius_bracket_m: bracket,
                    accepted_fraction: fraction,
                });
            }
        }
        let (source, response) =
            self.stage_successor(before, accepted[0], [accepted[1], accepted[2]], &mut budget)?;
        let actual = [
            response.state.normal_stage_m,
            source.tangential_stage[0],
            source.tangential_stage[1],
        ];
        let untraversed_stages_m = core::array::from_fn(|i| requested_stages_m[i] - actual[i]);
        let mechanical_work_j = response.total_stored_energy_j - before.total_stored_energy_j;
        for value in untraversed_stages_m.into_iter().chain([mechanical_work_j]) {
            finite(value)?;
        }
        budget.poll()?;
        Ok(JkrCombinedStageTransition {
            original_fixture_identity: self.identity,
            source,
            before,
            response,
            event,
            requested_stages_m,
            untraversed_stages_m,
            mechanical_work_j,
            evaluations: budget.evaluations,
        })
    }

    fn combined_stage_limit(
        &self,
        before: JkrPullOffResponse,
        path: &StagePath,
        budget: &mut Work<'_>,
    ) -> Result<(RadiusTrial, JkrShearBoundaryKind, [f64; 2]), JkrPullOffError> {
        let radius = before
            .contact
            .ok_or(JkrPullOffError::StateMismatch)?
            .contact_radius_m;
        let start = RadiusTrial {
            radius,
            fraction: 0.0,
            curvature: before.apparatus_radius_curvature_n_per_m,
            capacity: before.static_capacity_n,
            margin: before.static_margin_n,
        };
        if before.static_capacity_n == 0.0 && (path.increment[1] != 0.0 || path.increment[2] != 0.0)
        {
            return Ok((
                start,
                if before.stage_stiffness_n_per_m.is_none() {
                    JkrShearBoundaryKind::Coincident
                } else {
                    JkrShearBoundaryKind::GrossSlip
                },
                [radius; 2],
            ));
        }
        if before.state.normal_stage_m == self.boundary.response.state.normal_stage_m {
            return Ok((start, self.boundary.kind, self.boundary.radius_bracket_m));
        }
        if start.curvature <= 0.0 {
            return Err(JkrPullOffError::NoConvergence);
        }
        let (fold, fold_bracket) =
            path.last_admitted(0.0, start, false, self.limits.maximum_iterations, budget)?;
        // The zero-shear path has no active gross-slip constraint, even if the
        // original source declares zero junction strength.
        let shear_active = path.from[1] != 0.0
            || path.from[2] != 0.0
            || path.increment[1] != 0.0
            || path.increment[2] != 0.0;
        let uncertainty = 64.0 * f64::EPSILON * fold.capacity;
        if shear_active && fold.margin < -uncertainty {
            let (slip, bracket) = path.last_admitted(
                fold.radius,
                start,
                true,
                self.limits.maximum_iterations,
                budget,
            )?;
            Ok((
                slip,
                if bracket[0] > fold_bracket[1] {
                    JkrShearBoundaryKind::GrossSlip
                } else {
                    JkrShearBoundaryKind::Coincident
                },
                bracket,
            ))
        } else {
            Ok((
                fold,
                if shear_active && fold.margin.abs() <= uncertainty {
                    JkrShearBoundaryKind::Coincident
                } else {
                    JkrShearBoundaryKind::RadiusInstability
                },
                fold_bracket,
            ))
        }
    }
}

struct StagePath {
    e: f64,
    w: f64,
    r: f64,
    k: f64,
    kn: f64,
    kt: f64,
    strength: f64,
    from: [f64; 3],
    increment: [f64; 3],
}
#[derive(Clone, Copy)]
struct RadiusTrial {
    radius: f64,
    fraction: f64,
    curvature: f64,
    capacity: f64,
    margin: f64,
}
impl StagePath {
    fn new(
        source: &JkrPullOffFixture,
        from: [f64; 3],
        to: [f64; 3],
    ) -> Result<Self, JkrPullOffError> {
        let increment = core::array::from_fn(|i| to[i] - from[i]);
        for value in increment {
            finite(value)?;
        }
        if increment[0] > 0.0 {
            return Err(JkrPullOffError::Invalid(
                "combined attached loading requires normal withdrawal",
            ));
        }
        let scale = from[1]
            .abs()
            .max(from[2].abs())
            .max(to[1].abs())
            .max(to[2].abs());
        if scale > 0.0 {
            let radial = (from[1] / scale) * (increment[1] / scale)
                + (from[2] / scale) * (increment[2] / scale);
            if radial < 0.0 {
                return Err(JkrPullOffError::Invalid(
                    "combined attached loading requires nondecreasing shear magnitude",
                ));
            }
        }
        let normal = source.contact.normal_source();
        Ok(Self {
            e: normal.normal_state().elastic().reduced_modulus_pa(),
            w: normal.normal_state().adhesion_energy_j_per_m2(),
            r: normal.effective_radius_m(),
            k: source.contact.tangential_stiffness_per_radius,
            kn: source.normal_stiffness,
            kt: source.tangential_stiffness,
            strength: source.shear.law().static_shear_strength_pa(),
            from,
            increment,
        })
    }

    /// Upper quadratic root is the only possible solution with g>0 along this
    /// outward path. None marks the region below its tensile/real-root domain;
    /// such a probe is never published or labelled as a physical event.
    fn sample(
        &self,
        radius: f64,
        budget: &mut Work<'_>,
    ) -> Result<Option<RadiusTrial>, JkrPullOffError> {
        budget.probe()?;
        let c = 2.0 * self.e / self.kn;
        let q = 1.0 + c * radius;
        let hertz = radius * (radius / self.r);
        let g0 = hertz - self.from[0] / q - (c * radius / q) * hertz / 3.0;
        let dg = -self.increment[0] / q;
        let shear_scale = self.kt / (self.kt + self.k * radius);
        let u0 = [self.from[1] * shear_scale, self.from[2] * shear_scale];
        let du = [
            self.increment[1] * shear_scale,
            self.increment[2] * shear_scale,
        ];
        let aa = (self.e * dg) * dg + 0.5 * self.k * (du[0] * du[0] + du[1] * du[1]);
        let bb = (self.e * g0) * dg + 0.5 * self.k * (u0[0] * du[0] + u0[1] * du[1]);
        let cc = (self.e * g0) * g0 + 0.5 * self.k * (u0[0] * u0[0] + u0[1] * u0[1])
            - 2.0 * PI * self.w * radius;
        for value in [q, g0, dg, aa, bb, cc] {
            finite(value)?;
        }
        let scale = aa.abs().max(bb.abs()).max(cc.abs());
        let (a, b, c0) = (aa / scale, bb / scale, cc / scale);
        if !a.is_finite() || a <= 0.0 {
            return Err(JkrPullOffError::NoConvergence);
        }
        let discriminant = b.mul_add(b, -a * c0);
        if discriminant < 0.0 {
            return Ok(None);
        }
        let root = discriminant.sqrt();
        let fraction = if b >= 0.0 {
            if b + root == 0.0 {
                0.0
            } else {
                -c0 / (b + root)
            }
        } else {
            (root - b) / a
        };
        finite(fraction)?;
        let g = dg.mul_add(fraction, g0);
        if fraction < 0.0 || g <= 0.0 {
            return Ok(None);
        }
        let u = [
            du[0].mul_add(fraction, u0[0]),
            du[1].mul_add(fraction, u0[1]),
        ];
        let u_norm = u[0].hypot(u[1]);
        let curvature = 2.0 * self.e * g * (2.0 * radius / self.r - c * g / q)
            - (self.k * u_norm) * (self.k * u_norm / (self.kt + self.k * radius))
            - 2.0 * PI * self.w;
        let capacity = (PI * self.strength * radius) * radius;
        let margin = capacity - self.k * radius * u_norm;
        for value in [g, curvature, capacity, margin] {
            finite(value)?;
        }
        Ok(Some(RadiusTrial {
            radius,
            fraction,
            curvature,
            capacity,
            margin,
        }))
    }

    fn last_admitted(
        &self,
        mut lo: f64,
        mut hi: RadiusTrial,
        slip: bool,
        maximum_iterations: usize,
        budget: &mut Work<'_>,
    ) -> Result<(RadiusTrial, [f64; 2]), JkrPullOffError> {
        for _ in 0..maximum_iterations {
            let mid = f64::midpoint(lo, hi.radius);
            if mid == lo || mid == hi.radius || hi.radius - lo <= 16.0 * f64::EPSILON * hi.radius {
                return Ok((hi, [lo, hi.radius]));
            }
            let Some(trial) = self.sample(mid, budget)? else {
                if slip {
                    return Err(JkrPullOffError::NoConvergence);
                }
                lo = mid;
                continue;
            };
            if (if slip { trial.margin } else { trial.curvature }) < 0.0 {
                lo = mid;
            } else {
                hi = trial;
            }
        }
        Err(JkrPullOffError::NoConvergence)
    }
}

fn norm(value: [f64; 3]) -> f64 {
    value[0].hypot(value[1]).hypot(value[2])
}
fn interpolate(from: [f64; 3], to: [f64; 3], fraction: f64) -> [f64; 3] {
    if fraction == 0.0 {
        from
    } else if fraction == 1.0 {
        to
    } else {
        core::array::from_fn(|i| {
            if from[i] == to[i] {
                from[i]
            } else {
                fraction.mul_add(to[i], (1.0 - fraction) * from[i])
            }
        })
    }
}
