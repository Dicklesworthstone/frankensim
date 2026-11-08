//! Work-conjugate normal, two-axis shear and twist on one relaxed JKR area.
//!
//! This is the discrete gradient of the original combined potential, including
//! its surface term. Polynomial differences and the exact radial-stationarity
//! secant eliminate radius work without quadrature, finite differences or
//! subtracting nearly equal endpoint energies. The endpoint tangent follows
//! implicit equilibrium, not iterations of the radius solver.
//!
//! A coordinate box encloses every admitted straight path. The motion owner
//! supplies the extrema of its actual quadratic trajectory to the same box
//! admission. A box touching a fold, an opening edge or a source envelope
//! refuses; callers can subdivide. These are conservative floating-point
//! screens, not certified interval bounds. There is no sliding, irreversible
//! loss, heat conversion or post-detachment state in this full-stick owner.

use super::{JkrCombinedError, JkrCombinedResponse, SourceJkrMixedMode};
use core::f64::consts::PI;
use fs_exec::CancelGate;

#[derive(Clone, Debug, PartialEq)]
pub enum JkrCombinedPathError {
    Source(JkrCombinedError),
    Cancelled,
    /// No further original-source query may start.
    Budget,
    /// Marginal endpoint or unresolved radial secant.
    UnresolvedSecant,
    /// The complete coordinate box cannot establish attached contact.
    UnresolvedEnvelope,
}
impl From<JkrCombinedError> for JkrCombinedPathError {
    fn from(value: JkrCombinedError) -> Self {
        Self::Source(value)
    }
}
impl core::fmt::Display for JkrCombinedPathError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "original combined JKR path refused: {self:?}")
    }
}
impl std::error::Error for JkrCombinedPathError {}

/// Coordinate order [approach (m), shear1 (m), shear2 (m), twist (rad)].
/// Efforts enter the original potential; the body reactions are opposite.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrCombinedPath {
    pub before: JkrCombinedResponse,
    pub contact: JkrCombinedResponse,
    /// Work-conjugate [normal force, shear1 force, shear2 force, twist moment].
    pub mean_effort: [f64; 4],
    /// Derivative of mean efforts with respect to the endpoint coordinates.
    /// Generally nonsymmetric, unlike the pointwise condensed stiffness.
    pub end_tangent: [[f64; 4]; 4],
    /// Factored change of the complete original contact potential [J].
    pub reversible_work_j: f64,
    /// Mean effort dot increment minus the original potential work [J].
    pub work_residual_j: f64,
    pub minimum_contact_radius_m: f64,
    pub maximum_contact_radius_m: f64,
    /// Actual original-source queries, at most five for an uncached path.
    pub evaluations: usize,
}

/// Sealed original endpoint, retained without allocations for repeated trials.
/// Mutable diagnostic responses cannot be supplied as constitutive history.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PreparedJkrCombinedPath<'a> {
    source: &'a SourceJkrMixedMode,
    contact: JkrCombinedResponse,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PreparedJkrCombinedTrial<'a> {
    endpoint: PreparedJkrCombinedPath<'a>,
    path: JkrCombinedPath,
}
impl<'a> PreparedJkrCombinedPath<'a> {
    pub const fn coordinates(&self) -> [f64; 4] {
        self.contact.coordinates
    }
    pub const fn source(&self) -> &'a SourceJkrMixedMode {
        self.source
    }
    pub fn trial(
        &self,
        to: [f64; 4],
        gate: &CancelGate,
    ) -> Result<PreparedJkrCombinedTrial<'a>, JkrCombinedPathError> {
        let path = self
            .source
            .combined_path_from_response(self.contact, to, gate)?;
        Ok(PreparedJkrCombinedTrial {
            endpoint: Self {
                source: self.source,
                contact: path.contact,
            },
            path,
        })
    }
}
impl<'a> PreparedJkrCombinedTrial<'a> {
    pub const fn path(&self) -> &JkrCombinedPath {
        &self.path
    }
    pub const fn endpoint(&self) -> PreparedJkrCombinedPath<'a> {
        self.endpoint
    }
    pub fn into_parts(self) -> (PreparedJkrCombinedPath<'a>, JkrCombinedPath) {
        (self.endpoint, self.path)
    }
}

impl SourceJkrMixedMode {
    pub fn prepare_combined_discrete_path(
        &self,
        from: [f64; 4],
        gate: &CancelGate,
    ) -> Result<PreparedJkrCombinedPath<'_>, JkrCombinedPathError> {
        let path = self.combined_discrete_path(from, from, gate)?;
        Ok(PreparedJkrCombinedPath {
            source: self,
            contact: path.contact,
        })
    }

    /// Symmetric second-order discrete gradient on one relaxed contact area.
    /// Fixed four-component arithmetic and at most five source queries; no
    /// allocations, numerical quadrature or independent shear/twist caps.
    pub fn combined_discrete_path(
        &self,
        from: [f64; 4],
        to: [f64; 4],
        gate: &CancelGate,
    ) -> Result<JkrCombinedPath, JkrCombinedPathError> {
        poll(gate)?;
        let before = self.evaluate_combined_displacement_controlled(from)?;
        let mut path = self.combined_path_from_response(before, to, gate)?;
        path.evaluations += 1;
        Ok(path)
    }

    fn combined_path_from_response(
        &self,
        before: JkrCombinedResponse,
        to: [f64; 4],
        gate: &CancelGate,
    ) -> Result<JkrCombinedPath, JkrCombinedPathError> {
        let mut evaluations = 0;
        self.combined_path_counted(before, to, &mut evaluations, usize::MAX, gate)
    }

    /// Only the sibling motion owner can supply a sealed cached endpoint.
    /// Rejected proposals retain the actual number of source queries consumed.
    pub(super) fn combined_path_counted(
        &self,
        before: JkrCombinedResponse,
        to: [f64; 4],
        evaluations: &mut usize,
        maximum_evaluations: usize,
        gate: &CancelGate,
    ) -> Result<JkrCombinedPath, JkrCombinedPathError> {
        poll(gate)?;
        let started = *evaluations;
        let from = before.coordinates;
        let tangent = strictly_attached(before)?;
        if from == to {
            return Ok(JkrCombinedPath {
                before,
                contact: before,
                mean_effort: before.force,
                end_tangent: tangent.map(|row| row.map(|x| 0.5 * x)),
                reversible_work_j: 0.0,
                work_residual_j: 0.0,
                minimum_contact_radius_m: before.contact_radius_m,
                maximum_contact_radius_m: before.contact_radius_m,
                evaluations: 0,
            });
        }
        charge(evaluations, maximum_evaluations)?;
        let contact = self.evaluate_combined_displacement_controlled(to)?;
        strictly_attached(contact)?;
        let bounds = [
            core::array::from_fn(|i| from[i].min(to[i])),
            core::array::from_fn(|i| from[i].max(to[i])),
        ];
        let [minimum_contact_radius_m, maximum_contact_radius_m] = self.admit_combined_box(
            before,
            contact,
            bounds,
            evaluations,
            maximum_evaluations,
            gate,
        )?;
        let e = self.normal.normal_state().elastic().reduced_modulus_pa();
        let k = 8.0 * self.reduced_shear_modulus_pa();
        let b = (8.0 / 3.0) * self.torsional_reduced_shear_modulus_pa();
        let r = self.normal.effective_radius_m();
        let w = self.normal.normal_state().adhesion_energy_j_per_m2();
        let a0 = before.contact_radius_m;
        let a1 = contact.contact_radius_m;
        let a = a0.midpoint(a1);
        let da = a1 - a0;
        let q: [f64; 4] = core::array::from_fn(|i| from[i].midpoint(to[i]));
        let dq: [f64; 4] = core::array::from_fn(|i| to[i] - from[i]);
        let radius_square = a * a + 0.25 * da * da;
        let radius_cube = a * (a * a + 0.75 * da * da);
        let twist_square = q[3] * q[3] + 0.25 * dq[3] * dq[3];
        let opening = (a0 * (a0 / r)).midpoint(a1 * (a1 / r)) - q[0];
        // A=U_a. Its exact secant is A1-A0 = C*da + cross dot dq.
        let curvature = 4.0 * e * (a / r) * opening - 2.0 * PI * w + 6.0 * b * a * twist_square;
        if !curvature.is_finite() || curvature <= 0.0 {
            return Err(JkrCombinedPathError::UnresolvedSecant);
        }
        let cross = [
            -2.0 * e * opening,
            k * q[1],
            k * q[2],
            6.0 * b * radius_square * q[3],
        ];
        let normal_factor = q[0] - a * (a / r) / 3.0 - da * (da / r) / 4.0;
        let uncondensed = [
            2.0 * e * a * normal_factor,
            k * a * q[1],
            k * a * q[2],
            2.0 * b * radius_cube * q[3],
        ];
        let radial_factor = q[0] / (3.0 * r) - (a / r).powi(2) - (da / r).powi(2) / 20.0;
        let normal_inner = (a / r) * dq[0] + da * radial_factor;
        let torsion_inner = 3.0 * a * q[3] * dq[3] + 0.5 * da * twist_square;
        // For k*a*|u|^2/2 the radial discrete gradient is exactly the mean
        // endpoint U_a contribution. The original full residual retains it.
        let radial = before.radial_residual_n.midpoint(contact.radial_residual_n)
            + da * (e * normal_inner - b * torsion_inner);
        let ratio = radial / curvature;
        let mean_effort = core::array::from_fn(|i| uncondensed[i] - ratio * cross[i]);
        let reversible_work_j = dot(uncondensed, dq) + radial * da;
        let work_residual_j = dot(mean_effort, dq) - reversible_work_j;

        // Implicit derivative of ORIGINAL equilibrium; no root differentiation.
        let endpoint_cross = [
            -2.0 * e * (a1 * (a1 / r) - to[0]),
            k * to[1],
            k * to[2],
            6.0 * b * a1 * a1 * to[3],
        ];
        let radius_gradient =
            endpoint_cross.map(|x| -x / contact.displacement_radius_curvature_n_per_m);
        let mut end_tangent = [[0.0; 4]; 4];
        for j in 0..4 {
            let normal_coordinate = if j == 0 { 1.0 } else { 0.0 };
            let twist_coordinate = if j == 3 { 1.0 } else { 0.0 };
            let dr = radius_gradient[j];
            let opening_derivative = (a1 / r) * dr - 0.5 * normal_coordinate;
            let curvature_derivative = (4.0 * e / r)
                * (0.5 * dr * opening + a * opening_derivative)
                + 6.0 * b * (0.5 * dr * twist_square + a * to[3] * twist_coordinate);
            let factor_derivative =
                normal_coordinate / (6.0 * r) - (a / r) * (dr / r) - (da / r) * (dr / r) / 10.0;
            let normal_derivative = 0.5 * (dr / r) * dq[0]
                + (a / r) * normal_coordinate
                + dr * radial_factor
                + da * factor_derivative;
            let torsion_derivative = 1.5 * dr * q[3] * dq[3]
                + 3.0 * a * (0.5 * dq[3] + q[3]) * twist_coordinate
                + 0.5 * dr * twist_square
                + 0.5 * da * to[3] * twist_coordinate;
            let radial_derivative = e * (dr * normal_inner + da * normal_derivative)
                - b * (dr * torsion_inner + da * torsion_derivative);
            let ratio_derivative = (radial_derivative - ratio * curvature_derivative) / curvature;
            let base_derivative = [
                2.0 * e
                    * (0.5 * dr * normal_factor
                        + a * (0.5 * normal_coordinate - (a / (3.0 * r) + da / (2.0 * r)) * dr)),
                0.5 * k * (dr * q[1] + if j == 1 { a } else { 0.0 }),
                0.5 * k * (dr * q[2] + if j == 2 { a } else { 0.0 }),
                b * (3.0 * a1 * a1 * dr * q[3] + radius_cube * twist_coordinate),
            ];
            let cross_derivative = [
                -2.0 * e * opening_derivative,
                if j == 1 { 0.5 * k } else { 0.0 },
                if j == 2 { 0.5 * k } else { 0.0 },
                6.0 * b * (a1 * dr * q[3] + 0.5 * radius_square * twist_coordinate),
            ];
            for i in 0..4 {
                end_tangent[i][j] =
                    base_derivative[i] - ratio_derivative * cross[i] - ratio * cross_derivative[i];
            }
        }
        if !reversible_work_j.is_finite()
            || !work_residual_j.is_finite()
            || mean_effort.iter().any(|x| !x.is_finite())
            || end_tangent.iter().flatten().any(|x| !x.is_finite())
        {
            return Err(JkrCombinedPathError::UnresolvedSecant);
        }
        Ok(JkrCombinedPath {
            before,
            contact,
            mean_effort,
            end_tangent,
            reversible_work_j,
            work_residual_j,
            minimum_contact_radius_m,
            maximum_contact_radius_m,
            evaluations: *evaluations - started,
        })
    }

    /// Shared source admission for straight paths and inertial trajectories.
    /// Corners are componentwise [minimum,maximum] coordinates.
    pub(super) fn admit_combined_box(
        &self,
        before: JkrCombinedResponse,
        contact: JkrCombinedResponse,
        coordinates: [[f64; 4]; 2],
        evaluations: &mut usize,
        maximum_evaluations: usize,
        gate: &CancelGate,
    ) -> Result<[f64; 2], JkrCombinedPathError> {
        if coordinates.iter().flatten().any(|x| !x.is_finite())
            || (0..4).any(|i| coordinates[0][i] > coordinates[1][i])
        {
            return Err(JkrCombinedPathError::UnresolvedEnvelope);
        }
        let low = coordinates[0][0];
        let high = coordinates[1][0];
        let maximum: [f64; 3] =
            core::array::from_fn(|i| coordinates[0][i + 1].abs().max(coordinates[1][i + 1].abs()));
        let minimum: [f64; 3] = core::array::from_fn(|i| {
            let lo = coordinates[0][i + 1];
            let hi = coordinates[1][i + 1];
            if lo <= 0.0 && hi >= 0.0 {
                0.0
            } else {
                lo.abs().min(hi.abs())
            }
        });
        let mut evaluate = |q: [f64; 4]| -> Result<JkrCombinedResponse, JkrCombinedPathError> {
            poll(gate)?;
            for cached in [before, contact] {
                if cached.coordinates[0] == q[0]
                    && (1..4).all(|i| cached.coordinates[i].abs() == q[i])
                {
                    return Ok(cached);
                }
            }
            charge(evaluations, maximum_evaluations)?;
            let point = self
                .evaluate_combined_displacement_controlled(q)
                .map_err(|_| JkrCombinedPathError::UnresolvedEnvelope)?;
            strictly_attached(point).map_err(|_| JkrCombinedPathError::UnresolvedEnvelope)?;
            Ok(point)
        };
        // At fixed shear/twist, delta(a)=a^2/R-sqrt(Q(a)/E), where Q is
        // a concave quadratic. Hence delta''(a)>0 on the opening branch:
        // its stable interval has no hidden gap between these approaches.
        // At fixed delta, radius decreases with either loading magnitude.
        let lower = evaluate([low, maximum[0], maximum[1], maximum[2]])?;
        evaluate([high, maximum[0], maximum[1], maximum[2]])?;
        let upper = evaluate([high, minimum[0], minimum[1], minimum[2]])?;
        let bounds = [
            lower.contact_radius_m.next_down(),
            upper.contact_radius_m.next_up(),
        ];
        if bounds[0] <= 0.0
            || bounds[1] < bounds[0]
            || bounds
                .iter()
                .any(|&a| self.normal.check_applicability(a).is_err())
        {
            return Err(JkrCombinedPathError::UnresolvedEnvelope);
        }
        // Triangle bounds apply to the combined rim displacement and mean
        // traction, not independently admitted shear and twist quantities.
        let shear = maximum[0].hypot(maximum[1]);
        let ratio = (shear / bounds[0] + maximum[2]).next_up();
        let yield_stress = self
            .normal
            .normal_state()
            .elastic()
            .limiting_yield_stress_pa();
        let mean_stress = (8.0 * self.reduced_shear_modulus_pa() * (shear / bounds[0]) / PI
            + 2.0 * self.torsional_reduced_shear_modulus_pa() * maximum[2])
            .next_up();
        if !ratio.is_finite()
            || !mean_stress.is_finite()
            || ratio > self.envelope.maximum_tangential_displacement_to_radius
            || mean_stress / yield_stress > self.envelope.maximum_mean_shear_stress_to_yield
        {
            return Err(JkrCombinedPathError::UnresolvedEnvelope);
        }
        Ok(bounds)
    }
}

fn strictly_attached(point: JkrCombinedResponse) -> Result<[[f64; 4]; 4], JkrCombinedPathError> {
    if point.boundary.is_some()
        || !point.displacement_radius_curvature_n_per_m.is_finite()
        || point.displacement_radius_curvature_n_per_m <= 0.0
    {
        return Err(JkrCombinedPathError::UnresolvedSecant);
    }
    point
        .coupled_stiffness
        .ok_or(JkrCombinedPathError::UnresolvedSecant)
}
fn poll(gate: &CancelGate) -> Result<(), JkrCombinedPathError> {
    if gate.is_requested() {
        Err(JkrCombinedPathError::Cancelled)
    } else {
        Ok(())
    }
}
fn charge(evaluations: &mut usize, maximum: usize) -> Result<(), JkrCombinedPathError> {
    if *evaluations >= maximum {
        return Err(JkrCombinedPathError::Budget);
    }
    *evaluations += 1;
    Ok(())
}
fn dot(a: [f64; 4], b: [f64; 4]) -> f64 {
    (a[0] * b[0] + a[1] * b[1]) + (a[2] * b[2] + a[3] * b[3])
}
