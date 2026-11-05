//! Symmetric discrete gradient of the original full-stick Griffith potential.
//!
//! The uncondensed potential is polynomial in (a,delta,u). Its symmetric
//! product differences give exact work without subtracting endpoint energies.
//! The exact secant of radius stationarity eliminates the radial work. This
//! avoids both quadrature and division by a vanishing coordinate increment.
//! Original finite equilibrium residuals remain visible in the work residual.
//!
//! A separate bounded quadratic-path admission encloses every point using
//! Bezier convex hulls and monotonic radius corners. It never accepts a path
//! by sampling it. Subdivision tightens an inconclusive enclosure; exhaustion
//! refuses. These are conservative floating-point screens, not interval
//! certificates or a post-slip, post-fold, or bond-nucleation law.

use super::{JkrMixedModeControl, JkrMixedModeError, JkrMixedModeResponse, SourceJkrMixedMode};
use crate::interface_binding::adhesive_shear::{AdhesiveShearError, SourceAdhesiveShear};
use core::f64::consts::PI;
use fs_exec::CancelGate;

#[derive(Clone, Debug, PartialEq)]
pub enum JkrMixedModePathError {
    Source(JkrMixedModeError),
    Shear(AdhesiveShearError),
    Invalid(&'static str),
    StateMismatch,
    Cancelled,
    Budget,
    UnresolvedEnvelope,
    UnresolvedSecant,
    StaticCapacity { force_n: f64, capacity_n: f64 },
}
impl From<JkrMixedModeError> for JkrMixedModePathError {
    fn from(value: JkrMixedModeError) -> Self {
        Self::Source(value)
    }
}
impl From<AdhesiveShearError> for JkrMixedModePathError {
    fn from(value: AdhesiveShearError) -> Self {
        Self::Shear(value)
    }
}
impl core::fmt::Display for JkrMixedModePathError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "original mixed-mode JKR path refused: {self:?}")
    }
}
impl std::error::Error for JkrMixedModePathError {}

/// Endpoint work exchange; whole physical-path admission is a separate call.
/// Every effort is positive into the original contact potential. Native body
/// reactions have the opposite sign. The endpoint derivative is generally
/// nonsymmetric even though each original constitutive stiffness is symmetric.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrMixedModePath {
    pub before: JkrMixedModeResponse,
    pub contact: JkrMixedModeResponse,
    pub mean_force_n: [f64; 3],
    pub end_tangent_n_per_m: [[f64; 3]; 3],
    /// Factored difference of the ORIGINAL normal/shear/surface potential.
    pub reversible_work_j: f64,
    /// Mean effort dot coordinate increment minus original potential work.
    pub work_residual_j: f64,
    pub evaluations: usize,
}

/// Immutable original endpoint for repeated native geometry/Newton proposals.
/// Only this source can create or advance the cached response. Caller-editable
/// diagnostic fields never become constitutive state. Preparation costs one
/// original equilibrium query; each changed trial endpoint costs one more.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PreparedJkrMixedModePath<'a> {
    source: &'a SourceJkrMixedMode,
    contact: JkrMixedModeResponse,
}

/// A sealed successor accompanies the ordinary work/tangent diagnostics.
/// Whole-path and static-capacity admission remain the physical caller's job.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PreparedJkrMixedModeTrial<'a> {
    endpoint: PreparedJkrMixedModePath<'a>,
    path: JkrMixedModePath,
}

impl<'a> PreparedJkrMixedModePath<'a> {
    pub const fn source(&self) -> &'a SourceJkrMixedMode {
        self.source
    }
    pub const fn contact(&self) -> &JkrMixedModeResponse {
        &self.contact
    }
    pub fn position_m(&self) -> [f64; 3] {
        [
            self.contact.indentation_m,
            self.contact.tangential_displacement_m[0],
            self.contact.tangential_displacement_m[1],
        ]
    }

    /// Propose from the same original accepted radius without solving it again.
    /// Exact coincidence uses no source query. Every failure leaves this cache
    /// unchanged; accepting a trial is explicit through its sealed endpoint.
    pub fn trial(
        &self,
        to: [f64; 3],
        gate: &CancelGate,
    ) -> Result<PreparedJkrMixedModeTrial<'a>, JkrMixedModePathError> {
        let path = self
            .source
            .discrete_path_from_response(self.contact, to, gate)?;
        Ok(PreparedJkrMixedModeTrial {
            endpoint: Self {
                source: self.source,
                contact: path.contact,
            },
            path,
        })
    }
}

impl<'a> PreparedJkrMixedModeTrial<'a> {
    pub const fn endpoint(&self) -> PreparedJkrMixedModePath<'a> {
        self.endpoint
    }
    pub const fn path(&self) -> &JkrMixedModePath {
        &self.path
    }
    pub fn into_parts(self) -> (PreparedJkrMixedModePath<'a>, JkrMixedModePath) {
        (self.endpoint, self.path)
    }
}

/// Conservative bounds over every point of an admitted quadratic trajectory.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrMixedModePathBounds {
    pub evaluations: usize,
    /// Bound on the Euclidean norm of the original three-component effort.
    pub maximum_force_n: f64,
    pub minimum_radius_m: f64,
    pub maximum_radius_m: f64,
    pub maximum_shear_displacement_to_radius: f64,
}

impl SourceJkrMixedMode {
    /// Prepare an original, strictly stable endpoint for allocation-free trials.
    /// This performs one source query and the coincident radial-tangent check;
    /// it neither establishes a physical trajectory nor admits static capacity.
    pub fn prepare_discrete_path(
        &self,
        from: [f64; 3],
        gate: &CancelGate,
    ) -> Result<PreparedJkrMixedModePath<'_>, JkrMixedModePathError> {
        let path = self.discrete_path(from, from, gate)?;
        Ok(PreparedJkrMixedModePath {
            source: self,
            contact: path.contact,
        })
    }

    /// Symmetric, second-order, work-consistent endpoint force and its analytic
    /// derivative with respect to `to`, in [approach,shear1,shear2] order.
    /// Source and radial secant singularities refuse; no fold continuation is
    /// inferred. Call `admit_quadratic_path` to admit an inertial trajectory and
    /// its actual sourced static capacity before publishing a physical step.
    pub fn discrete_path(
        &self,
        from: [f64; 3],
        to: [f64; 3],
        gate: &CancelGate,
    ) -> Result<JkrMixedModePath, JkrMixedModePathError> {
        poll(gate)?;
        let before = self.evaluate_displacement_controlled(from[0], [from[1], from[2]])?;
        let mut path = self.discrete_path_from_response(before, to, gate)?;
        path.evaluations += 1;
        Ok(path)
    }

    /// Reuse an immutable ORIGINAL response already evaluated by this module's
    /// motion or prepared-path owner. This is not a public response-taking API:
    /// caller-editable diagnostic fields cannot become constitutive authority.
    /// Counts only new source evaluations (zero at coincidence, one otherwise).
    pub(super) fn discrete_path_from_response(
        &self,
        before: JkrMixedModeResponse,
        to: [f64; 3],
        gate: &CancelGate,
    ) -> Result<JkrMixedModePath, JkrMixedModePathError> {
        poll(gate)?;
        if before.model_identity() != self.identity()
            || before.normal_state_identity() != self.normal.normal_state().identity()
            || before.control != JkrMixedModeControl::Displacement
        {
            return Err(JkrMixedModePathError::StateMismatch);
        }
        let from = [
            before.indentation_m,
            before.tangential_displacement_m[0],
            before.tangential_displacement_m[1],
        ];
        let first_tangent = before
            .coupled_stiffness_n_per_m
            .ok_or(JkrMixedModePathError::UnresolvedSecant)?;
        if from == to {
            return Ok(JkrMixedModePath {
                before,
                contact: before,
                mean_force_n: forces(before),
                end_tangent_n_per_m: first_tangent.map(|row| row.map(|x| 0.5 * x)),
                reversible_work_j: 0.0,
                work_residual_j: 0.0,
                evaluations: 0,
            });
        }
        poll(gate)?;
        let contact = self.evaluate_displacement_controlled(to[0], [to[1], to[2]])?;
        if contact.coupled_stiffness_n_per_m.is_none() {
            return Err(JkrMixedModePathError::UnresolvedSecant);
        }
        let e = self.normal.normal_state().elastic().reduced_modulus_pa();
        let k = self.tangential_stiffness_per_radius;
        let r = self.normal.radius_m;
        let w = self.normal.normal_state().adhesion_energy_j_per_m2();
        let a0 = before.contact_radius_m;
        let a1 = contact.contact_radius_m;
        let a = a0.midpoint(a1);
        let da = a1 - a0;
        let q: [f64; 3] = core::array::from_fn(|i| from[i].midpoint(to[i]));
        let dq: [f64; 3] = core::array::from_fn(|i| to[i] - from[i]);
        // A = dU/da. Retain the original source's actual residual instead of
        // silently declaring approximate equilibrium to be exact arithmetic.
        let radial0 = 2.0 * PI * a0 * before.griffith_residual_j_m2;
        let radial1 = 2.0 * PI * a1 * contact.griffith_residual_j_m2;
        let opening = (a0 * (a0 / r)).midpoint(a1 * (a1 / r)) - q[0];
        let curvature = 4.0 * e * (a / r) * opening - 2.0 * PI * w;
        if !curvature.is_finite() || curvature <= 0.0 {
            return Err(JkrMixedModePathError::UnresolvedSecant);
        }
        // Exact stationarity secant: A1-A0 = C*da + c dot dq.
        let cross = [-2.0 * e * opening, k * q[1], k * q[2]];
        let depth = a * (a / r);
        let radius_square_step = da * (da / r);
        let normal_factor = q[0] - depth / 3.0 - radius_square_step / 4.0;
        let uncondensed = [2.0 * e * a * normal_factor, k * a * q[1], k * a * q[2]];
        let radial_factor = q[0] / (3.0 * r) - (a / r) * (a / r) - (da / r) * (da / r) / 20.0;
        let radial_inner = (a / r) * dq[0] + da * radial_factor;
        let radial = radial0.midpoint(radial1) + e * da * radial_inner;
        let ratio = radial / curvature;
        let mean_force_n = core::array::from_fn(|i| uncondensed[i] - ratio * cross[i]);
        let reversible_work_j = dot(uncondensed, dq) + radial * da;
        let work_residual_j = dot(mean_force_n, dq) - reversible_work_j;

        // Differentiate only the admitted equilibrium, never its root solver.
        // da1/dq = -U_aq/U_aa. Total derivatives of A0 and A1 are zero on this
        // branch; their finite arithmetic residual remains in value/work audits.
        let endpoint_opening = a1 * (a1 / r) - to[0];
        let endpoint_cross = [-2.0 * e * endpoint_opening, k * to[1], k * to[2]];
        let radius_gradient =
            endpoint_cross.map(|x| -x / contact.displacement_radius_curvature_n_per_m);
        let mut end_tangent_n_per_m = [[0.0; 3]; 3];
        for j in 0..3 {
            let normal_coordinate = if j == 0 { 1.0 } else { 0.0 };
            let dr = radius_gradient[j];
            let opening_derivative = (a1 / r) * dr - 0.5 * normal_coordinate;
            let curvature_derivative =
                (4.0 * e / r) * (0.5 * dr * opening + a * opening_derivative);
            let radial_factor_derivative =
                normal_coordinate / (6.0 * r) - (a / r) * (dr / r) - (da / r) * (dr / r) / 10.0;
            let inner_derivative = 0.5 * (dr / r) * dq[0]
                + (a / r) * normal_coordinate
                + dr * radial_factor
                + da * radial_factor_derivative;
            let radial_derivative = e * (dr * radial_inner + da * inner_derivative);
            let ratio_derivative = (radial_derivative - ratio * curvature_derivative) / curvature;
            let base_derivative = [
                2.0 * e
                    * (0.5 * dr * normal_factor
                        + a * (0.5 * normal_coordinate - (a / (3.0 * r) + da / (2.0 * r)) * dr)),
                0.5 * k * (dr * q[1] + if j == 1 { a } else { 0.0 }),
                0.5 * k * (dr * q[2] + if j == 2 { a } else { 0.0 }),
            ];
            let cross_derivative = [
                -2.0 * e * opening_derivative,
                if j == 1 { 0.5 * k } else { 0.0 },
                if j == 2 { 0.5 * k } else { 0.0 },
            ];
            for i in 0..3 {
                end_tangent_n_per_m[i][j] =
                    base_derivative[i] - ratio_derivative * cross[i] - ratio * cross_derivative[i];
            }
        }
        if !reversible_work_j.is_finite()
            || !work_residual_j.is_finite()
            || mean_force_n.iter().any(|x| !x.is_finite())
            || end_tangent_n_per_m.iter().flatten().any(|x| !x.is_finite())
        {
            return Err(JkrMixedModePathError::UnresolvedSecant);
        }
        Ok(JkrMixedModePath {
            before,
            contact,
            mean_force_n,
            end_tangent_n_per_m,
            reversible_work_j,
            work_residual_j,
            evaluations: 1,
        })
    }

    /// Admit an entire quadratic Bezier path, including its hidden interior
    /// excursions. For a midpoint inertial step the control points are
    /// `[q0, q0 + dt*v0/2, q1]`; a straight path uses the endpoint midpoint.
    /// All original normal/shear screens and the independently sourced static
    /// strength apply to the evolving REAL contact area. A failed enclosure
    /// subdivides deterministically under one shared evaluation budget.
    pub fn admit_quadratic_path(
        &self,
        shear: &SourceAdhesiveShear,
        bezier: [[f64; 3]; 3],
        maximum_evaluations: usize,
        gate: &CancelGate,
    ) -> Result<JkrMixedModePathBounds, JkrMixedModePathError> {
        poll(gate)?;
        if shear.normal_state_identity() != self.normal.normal_state().identity() {
            return Err(JkrMixedModePathError::StateMismatch);
        }
        if maximum_evaluations == 0 || bezier.iter().flatten().any(|x| !x.is_finite()) {
            return Err(JkrMixedModePathError::Invalid(
                "finite path and positive evaluation budget",
            ));
        }
        let mut work = Work {
            gate,
            maximum: maximum_evaluations,
            evaluations: 0,
        };
        let strength = shear.law().static_shear_strength_pa();
        for point in [bezier[0], bezier[2]] {
            self.admit_point(point, strength, &mut work)?;
        }
        let mut bounds = self.admit_hull(bezier, strength, 0, &mut work)?;
        bounds.evaluations = work.evaluations;
        Ok(bounds)
    }

    fn admit_point(
        &self,
        point: [f64; 3],
        strength: f64,
        work: &mut Work<'_>,
    ) -> Result<(), JkrMixedModePathError> {
        let response = work.evaluate(self, point)?;
        if response.marginal {
            return Err(JkrMixedModePathError::UnresolvedEnvelope);
        }
        let force_n = response.tangential_force_n[0].hypot(response.tangential_force_n[1]);
        let capacity_n = strength * response.contact_area_m2;
        if !capacity_n.is_finite() || force_n > capacity_n {
            return Err(JkrMixedModePathError::StaticCapacity {
                force_n,
                capacity_n,
            });
        }
        Ok(())
    }

    fn admit_hull(
        &self,
        points: [[f64; 3]; 3],
        strength: f64,
        depth: usize,
        work: &mut Work<'_>,
    ) -> Result<JkrMixedModePathBounds, JkrMixedModePathError> {
        poll(work.gate)?;
        if let Some(bounds) = self.hull_bounds(points, strength, work)? {
            return Ok(bounds);
        }
        if depth >= 32 {
            return Err(JkrMixedModePathError::UnresolvedEnvelope);
        }
        let left = midpoint(points[0], points[1]);
        let right = midpoint(points[1], points[2]);
        let middle = midpoint(left, right);
        // This actual point can disprove admission, but cannot prove a segment.
        self.admit_point(middle, strength, work)?;
        let lo = self.admit_hull([points[0], left, middle], strength, depth + 1, work)?;
        let hi = self.admit_hull([middle, right, points[2]], strength, depth + 1, work)?;
        Ok(JkrMixedModePathBounds {
            evaluations: 0,
            maximum_force_n: lo.maximum_force_n.max(hi.maximum_force_n),
            minimum_radius_m: lo.minimum_radius_m.min(hi.minimum_radius_m),
            maximum_radius_m: lo.maximum_radius_m.max(hi.maximum_radius_m),
            maximum_shear_displacement_to_radius: lo
                .maximum_shear_displacement_to_radius
                .max(hi.maximum_shear_displacement_to_radius),
        })
    }

    fn hull_bounds(
        &self,
        points: [[f64; 3]; 3],
        strength: f64,
        work: &mut Work<'_>,
    ) -> Result<Option<JkrMixedModePathBounds>, JkrMixedModePathError> {
        let low = points.iter().map(|p| p[0]).fold(f64::INFINITY, f64::min);
        let high = points
            .iter()
            .map(|p| p[0])
            .fold(f64::NEG_INFINITY, f64::max);
        let radii = points.map(|p| p[1].hypot(p[2]));
        let maximum_shear = radii.into_iter().fold(0.0_f64, f64::max);
        if !maximum_shear.is_finite() {
            return Ok(None);
        }
        // Every convex hull point lies in the ball centered at each vertex
        // with radius its furthest vertex. Reverse triangle inequality gives
        // a rotation-invariant lower norm bound which tightens on subdivision.
        let mut minimum_shear = 0.0_f64;
        for i in 0..3 {
            let span = points
                .iter()
                .map(|p| (p[1] - points[i][1]).hypot(p[2] - points[i][2]))
                .fold(0.0_f64, f64::max);
            minimum_shear = minimum_shear.max((radii[i] - span).next_down().max(0.0));
        }
        let maximum_shear = if maximum_shear == 0.0 {
            0.0
        } else {
            maximum_shear.next_up()
        };
        // On the stable tensile branch a_delta>0 and a_|u|²<0. These two
        // original corner solves enclose ALL radii in the entire hull. A bad
        // corner is inconclusive (it need not lie on the physical path).
        let lower = match work.evaluate(self, [low, maximum_shear, 0.0]) {
            Ok(response) if !response.marginal => response,
            Err(JkrMixedModePathError::Source(_)) | Ok(_) => return Ok(None),
            Err(error) => return Err(error),
        };
        let upper = match work.evaluate(self, [high, minimum_shear, 0.0]) {
            Ok(response) if !response.marginal => response,
            Err(JkrMixedModePathError::Source(_)) | Ok(_) => return Ok(None),
            Err(error) => return Err(error),
        };
        let a_min = lower.contact_radius_m.next_down();
        let a_max = upper.contact_radius_m.next_up();
        if a_min <= 0.0 || a_max < a_min {
            return Ok(None);
        }
        // All original radius screens are monotone or convex in a; their
        // maxima therefore occur at these two bounds (including normal strain).
        if self.normal.check_applicability(a_min).is_err()
            || self.normal.check_applicability(a_max).is_err()
        {
            return Ok(None);
        }
        let ratio = if maximum_shear == 0.0 {
            0.0
        } else {
            (maximum_shear / a_min).next_up()
        };
        let shear_stress = self.tangential_stiffness_per_radius * ratio / PI;
        let yield_stress = self
            .normal
            .normal_state()
            .elastic()
            .limiting_yield_stress_pa();
        if !ratio.is_finite()
            || !shear_stress.is_finite()
            || ratio > self.envelope.maximum_tangential_displacement_to_radius
            || shear_stress / yield_stress > self.envelope.maximum_mean_shear_stress_to_yield
            || shear_stress > strength
        {
            return Ok(None);
        }
        let e = self.normal.normal_state().elastic().reduced_modulus_pa();
        let normal_force = 2.0
            * e
            * a_max
            * (low.abs().max(high.abs()) + a_max * (a_max / self.normal.radius_m) / 3.0);
        let shear_force = self.tangential_stiffness_per_radius * a_max * maximum_shear;
        let maximum_force_n = normal_force.hypot(shear_force).next_up();
        if !maximum_force_n.is_finite() {
            return Ok(None);
        }
        Ok(Some(JkrMixedModePathBounds {
            evaluations: 0,
            maximum_force_n,
            minimum_radius_m: a_min,
            maximum_radius_m: a_max,
            maximum_shear_displacement_to_radius: ratio,
        }))
    }
}

struct Work<'a> {
    gate: &'a CancelGate,
    maximum: usize,
    evaluations: usize,
}
impl Work<'_> {
    fn evaluate(
        &mut self,
        source: &SourceJkrMixedMode,
        point: [f64; 3],
    ) -> Result<JkrMixedModeResponse, JkrMixedModePathError> {
        poll(self.gate)?;
        if self.evaluations >= self.maximum {
            return Err(JkrMixedModePathError::Budget);
        }
        self.evaluations += 1;
        Ok(source.evaluate_displacement_controlled(point[0], [point[1], point[2]])?)
    }
}
fn poll(gate: &CancelGate) -> Result<(), JkrMixedModePathError> {
    if gate.is_requested() {
        Err(JkrMixedModePathError::Cancelled)
    } else {
        Ok(())
    }
}
fn forces(response: JkrMixedModeResponse) -> [f64; 3] {
    [
        response.normal_force_n,
        response.tangential_force_n[0],
        response.tangential_force_n[1],
    ]
}
fn midpoint(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|i| a[i].midpoint(b[i]))
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}
