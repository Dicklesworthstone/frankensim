//! Symmetric discrete normal force and torque on the original relaxed area.
//!
//! This differentiates the SAME I/III potential as the pointwise owner. The
//! torsional term is `B a^3 phi^2`, `B=8 Gbar/3`, from Popov, Hess and Willert
//! (2019), DOI 10.1007/978-3-662-58709-6_5, Eq. 5.3. Symmetric polynomial
//! differences and the exact secant of radius stationarity eliminate radial
//! work without numerical quadrature or endpoint-energy subtraction. The
//! analytic endpoint tangent differentiates equilibrium, not its root solver.
//!
//! Admission encloses the straight coordinate path in its normal/twist box.
//! A box crossing a fold, opening edge or source envelope refuses, even if a
//! narrower path might fit; the caller can subdivide. This does not admit an
//! inertial trajectory with interior overshoot. There is no partial-slip,
//! post-fold, rate-dependent, irreversible-loss or heat model. Bounds are
//! conservative floating-point screens, not certified error enclosures.

use super::{JkrTorsionalError, JkrTorsionalResponse, SourceJkrMixedMode};
use core::f64::consts::PI;
use fs_exec::CancelGate;

#[derive(Clone, Debug, PartialEq)]
pub enum JkrTorsionalPathError {
    Source(JkrTorsionalError),
    Cancelled,
    /// Exhausted original-source queries before another evaluation started.
    Budget,
    /// Endpoint or radial secant is singular or numerically unresolved.
    UnresolvedSecant,
    /// The whole coordinate box cannot establish stable attached contact.
    UnresolvedEnvelope,
}
impl From<JkrTorsionalError> for JkrTorsionalPathError {
    fn from(value: JkrTorsionalError) -> Self {
        Self::Source(value)
    }
}
impl core::fmt::Display for JkrTorsionalPathError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "original torsional JKR path refused: {self:?}")
    }
}
impl std::error::Error for JkrTorsionalPathError {}

/// Efforts enter the original contact potential; body reactions are opposite.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct JkrTorsionalPath {
    pub before: JkrTorsionalResponse,
    pub contact: JkrTorsionalResponse,
    /// Work-conjugate [normal force (N), torsional moment (Nm)].
    pub mean_effort: [f64; 2],
    /// Derivative of mean efforts with respect to endpoint [approach (m),
    /// twist (rad)]. Generally nonsymmetric, unlike the pointwise stiffness.
    pub end_tangent: [[f64; 2]; 2],
    /// Factored change in normal, torsional and surface energy together [J].
    pub reversible_work_j: f64,
    /// Mean effort dot coordinate increment minus original potential work.
    pub work_residual_j: f64,
    pub minimum_contact_radius_m: f64,
    pub maximum_contact_radius_m: f64,
    /// Number of original-source queries, at most five for an uncached path.
    pub evaluations: usize,
}

/// Immutable original endpoint for repeated geometry/Newton proposals. Public
/// diagnostic fields cannot be fed back as constitutive state. Copying this
/// small value allocates nothing and never replaces its source or material.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PreparedJkrTorsionalPath<'a> {
    source: &'a SourceJkrMixedMode,
    contact: JkrTorsionalResponse,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PreparedJkrTorsionalTrial<'a> {
    endpoint: PreparedJkrTorsionalPath<'a>,
    path: JkrTorsionalPath,
}
impl<'a> PreparedJkrTorsionalPath<'a> {
    pub const fn coordinates(&self) -> [f64; 2] {
        [self.contact.indentation_m, self.contact.twist_angle_rad]
    }
    pub const fn source(&self) -> &'a SourceJkrMixedMode {
        self.source
    }
    pub fn trial(
        &self,
        to: [f64; 2],
        gate: &CancelGate,
    ) -> Result<PreparedJkrTorsionalTrial<'a>, JkrTorsionalPathError> {
        let path = self
            .source
            .torsional_path_from_response(self.contact, to, gate)?;
        Ok(PreparedJkrTorsionalTrial {
            endpoint: Self {
                source: self.source,
                contact: path.contact,
            },
            path,
        })
    }
}
impl<'a> PreparedJkrTorsionalTrial<'a> {
    pub const fn path(&self) -> &JkrTorsionalPath {
        &self.path
    }
    pub const fn endpoint(&self) -> PreparedJkrTorsionalPath<'a> {
        self.endpoint
    }
    pub fn into_parts(self) -> (PreparedJkrTorsionalPath<'a>, JkrTorsionalPath) {
        (self.endpoint, self.path)
    }
}

impl SourceJkrMixedMode {
    /// Admit and retain one original strictly attached normal/twist endpoint.
    pub fn prepare_torsional_discrete_path(
        &self,
        from: [f64; 2],
        gate: &CancelGate,
    ) -> Result<PreparedJkrTorsionalPath<'_>, JkrTorsionalPathError> {
        let path = self.torsional_discrete_path(from, from, gate)?;
        Ok(PreparedJkrTorsionalPath {
            source: self,
            contact: path.contact,
        })
    }

    /// Symmetric, second-order discrete gradient and its endpoint derivative.
    /// The whole straight [approach,twist] path must remain strictly attached
    /// inside the original material/geometry envelope. At most five original
    /// queries and fixed scalar arithmetic are required, with no allocations.
    pub fn torsional_discrete_path(
        &self,
        from: [f64; 2],
        to: [f64; 2],
        gate: &CancelGate,
    ) -> Result<JkrTorsionalPath, JkrTorsionalPathError> {
        poll(gate)?;
        let before = self.evaluate_torsion_displacement_controlled(from[0], from[1])?;
        let mut path = self.torsional_path_from_response(before, to, gate)?;
        path.evaluations += 1;
        Ok(path)
    }

    fn torsional_path_from_response(
        &self,
        before: JkrTorsionalResponse,
        to: [f64; 2],
        gate: &CancelGate,
    ) -> Result<JkrTorsionalPath, JkrTorsionalPathError> {
        let mut evaluations = 0;
        self.torsional_path_counted(before, to, &mut evaluations, usize::MAX, gate)
    }

    /// Only a sealed original endpoint may enter from the sibling motion
    /// owner. Rejected probes retain the actual number of queries consumed.
    pub(super) fn torsional_path_counted(
        &self,
        before: JkrTorsionalResponse,
        to: [f64; 2],
        evaluations: &mut usize,
        maximum_evaluations: usize,
        gate: &CancelGate,
    ) -> Result<JkrTorsionalPath, JkrTorsionalPathError> {
        poll(gate)?;
        let started = *evaluations;
        let from = [before.indentation_m, before.twist_angle_rad];
        let tangent = strictly_attached(before)?;
        if from == to {
            return Ok(JkrTorsionalPath {
                before,
                contact: before,
                mean_effort: [before.normal_force_n, before.torsional_moment_n_m],
                end_tangent: tangent.map(|row| row.map(|x| 0.5 * x)),
                reversible_work_j: 0.0,
                work_residual_j: 0.0,
                minimum_contact_radius_m: before.contact_radius_m,
                maximum_contact_radius_m: before.contact_radius_m,
                evaluations: 0,
            });
        }
        charge(evaluations, maximum_evaluations)?;
        let contact = self.evaluate_torsion_displacement_controlled(to[0], to[1])?;
        strictly_attached(contact)?;
        let bounds = [
            core::array::from_fn(|i| from[i].min(to[i])),
            core::array::from_fn(|i| from[i].max(to[i])),
        ];
        let [minimum_contact_radius_m, maximum_contact_radius_m] = self.admit_torsional_box(
            before,
            contact,
            bounds,
            evaluations,
            maximum_evaluations,
            gate,
        )?;
        let e = self.normal.normal_state().elastic().reduced_modulus_pa();
        let b = (8.0 / 3.0) * self.torsional_reduced_shear_modulus_pa();
        let r = self.normal.effective_radius_m();
        let w = self.normal.normal_state().adhesion_energy_j_per_m2();
        let a0 = before.contact_radius_m;
        let a1 = contact.contact_radius_m;
        let a = a0.midpoint(a1);
        let da = a1 - a0;
        let q: [f64; 2] = core::array::from_fn(|i| from[i].midpoint(to[i]));
        let dq: [f64; 2] = core::array::from_fn(|i| to[i] - from[i]);
        let radius_square = a * a + 0.25 * da * da;
        let radius_cube = a * (a * a + 0.75 * da * da);
        let twist_square = q[1] * q[1] + 0.25 * dq[1] * dq[1];
        let opening = (a0 * (a0 / r)).midpoint(a1 * (a1 / r)) - q[0];
        // A=U_a. Its exact secant is A1-A0 = C*da + cross dot dq.
        let curvature = 4.0 * e * (a / r) * opening - 2.0 * PI * w + 6.0 * b * a * twist_square;
        if !curvature.is_finite() || curvature <= 0.0 {
            return Err(JkrTorsionalPathError::UnresolvedSecant);
        }
        let cross = [-2.0 * e * opening, 6.0 * b * radius_square * q[1]];
        let normal_factor = q[0] - a * (a / r) / 3.0 - da * (da / r) / 4.0;
        let uncondensed = [2.0 * e * a * normal_factor, 2.0 * b * radius_cube * q[1]];
        let radial_factor = q[0] / (3.0 * r) - (a / r).powi(2) - (da / r).powi(2) / 20.0;
        let normal_inner = (a / r) * dq[0] + da * radial_factor;
        let torsion_inner = 3.0 * a * q[1] * dq[1] + 0.5 * da * twist_square;
        // Rewriting the radial discrete gradient around endpoint equilibrium
        // retains tiny work. The source's actual residual is not erased.
        let radial0 = 2.0 * PI * a0 * before.griffith_residual_j_m2;
        let radial1 = 2.0 * PI * a1 * contact.griffith_residual_j_m2;
        let radial = radial0.midpoint(radial1) + da * (e * normal_inner - b * torsion_inner);
        let ratio = radial / curvature;
        let mean_effort = core::array::from_fn(|i| uncondensed[i] - ratio * cross[i]);
        let reversible_work_j = self.torsional_endpoint_work(before, contact)?;
        let work_residual_j = dot(mean_effort, dq) - reversible_work_j;

        // Implicit differentiation of the endpoint equilibrium, da1/dq1.
        let endpoint_cross = [
            -2.0 * e * (a1 * (a1 / r) - to[0]),
            6.0 * b * a1 * a1 * to[1],
        ];
        let radius_gradient =
            endpoint_cross.map(|x| -x / contact.displacement_radius_curvature_n_per_m);
        let mut end_tangent = [[0.0; 2]; 2];
        for j in 0..2 {
            let normal_coordinate = if j == 0 { 1.0 } else { 0.0 };
            let twist_coordinate = if j == 1 { 1.0 } else { 0.0 };
            let dr = radius_gradient[j];
            let opening_derivative = (a1 / r) * dr - 0.5 * normal_coordinate;
            let curvature_derivative = (4.0 * e / r)
                * (0.5 * dr * opening + a * opening_derivative)
                + 6.0 * b * (0.5 * dr * twist_square + a * to[1] * twist_coordinate);
            let factor_derivative =
                normal_coordinate / (6.0 * r) - (a / r) * (dr / r) - (da / r) * (dr / r) / 10.0;
            let normal_derivative = 0.5 * (dr / r) * dq[0]
                + (a / r) * normal_coordinate
                + dr * radial_factor
                + da * factor_derivative;
            let torsion_derivative = 1.5 * dr * q[1] * dq[1]
                + 3.0 * a * (0.5 * dq[1] + q[1]) * twist_coordinate
                + 0.5 * dr * twist_square
                + 0.5 * da * to[1] * twist_coordinate;
            let radial_derivative = e * (dr * normal_inner + da * normal_derivative)
                - b * (dr * torsion_inner + da * torsion_derivative);
            let ratio_derivative = (radial_derivative - ratio * curvature_derivative) / curvature;
            let base_derivative = [
                2.0 * e
                    * (0.5 * dr * normal_factor
                        + a * (0.5 * normal_coordinate - (a / (3.0 * r) + da / (2.0 * r)) * dr)),
                b * (3.0 * a1 * a1 * dr * q[1] + radius_cube * twist_coordinate),
            ];
            let cross_derivative = [
                -2.0 * e * opening_derivative,
                6.0 * b * (a1 * dr * q[1] + 0.5 * radius_square * twist_coordinate),
            ];
            for i in 0..2 {
                end_tangent[i][j] =
                    base_derivative[i] - ratio_derivative * cross[i] - ratio * cross_derivative[i];
            }
        }
        if !reversible_work_j.is_finite()
            || !work_residual_j.is_finite()
            || mean_effort.iter().any(|x| !x.is_finite())
            || end_tangent.iter().flatten().any(|x| !x.is_finite())
        {
            return Err(JkrTorsionalPathError::UnresolvedSecant);
        }
        Ok(JkrTorsionalPath {
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

    /// Factored original potential change, including a marginal endpoint.
    /// This part of the discrete gradient needs no inverse radial curvature
    /// or endpoint tangent, so controlled loading can reach a genuine fold.
    /// The arguments must be original, source-authenticated responses.
    pub(super) fn torsional_endpoint_work(
        &self,
        before: JkrTorsionalResponse,
        contact: JkrTorsionalResponse,
    ) -> Result<f64, JkrTorsionalPathError> {
        let e = self.normal.normal_state().elastic().reduced_modulus_pa();
        let b = (8.0 / 3.0) * self.torsional_reduced_shear_modulus_pa();
        let r = self.normal.effective_radius_m();
        let a0 = before.contact_radius_m;
        let a1 = contact.contact_radius_m;
        let a = a0.midpoint(a1);
        let da = a1 - a0;
        let delta = before.indentation_m.midpoint(contact.indentation_m);
        let dd = contact.indentation_m - before.indentation_m;
        let phi = before.twist_angle_rad.midpoint(contact.twist_angle_rad);
        let dp = contact.twist_angle_rad - before.twist_angle_rad;
        let phi_square = phi * phi + 0.25 * dp * dp;
        let radius_cube = a * (a * a + 0.75 * da * da);
        let normal_factor = delta - a * (a / r) / 3.0 - da * (da / r) / 4.0;
        let radial_factor = delta / (3.0 * r) - (a / r).powi(2) - (da / r).powi(2) / 20.0;
        let normal_inner = (a / r) * dd + da * radial_factor;
        let torsion_inner = 3.0 * a * phi * dp + 0.5 * da * phi_square;
        let radial0 = 2.0 * PI * a0 * before.griffith_residual_j_m2;
        let radial1 = 2.0 * PI * a1 * contact.griffith_residual_j_m2;
        let radial = radial0.midpoint(radial1) + da * (e * normal_inner - b * torsion_inner);
        let work =
            (2.0 * e * a * normal_factor) * dd + (2.0 * b * radius_cube * phi) * dp + radial * da;
        if work.is_finite() {
            Ok(work)
        } else {
            Err(JkrTorsionalPathError::UnresolvedSecant)
        }
    }

    /// Original radius/applicability bounds shared by straight and inertial
    /// paths. `coordinates` holds componentwise [minimum, maximum] corners.
    pub(super) fn admit_torsional_box(
        &self,
        before: JkrTorsionalResponse,
        contact: JkrTorsionalResponse,
        coordinates: [[f64; 2]; 2],
        evaluations: &mut usize,
        maximum_evaluations: usize,
        gate: &CancelGate,
    ) -> Result<[f64; 2], JkrTorsionalPathError> {
        if coordinates.iter().flatten().any(|x| !x.is_finite())
            || (0..2).any(|i| coordinates[0][i] > coordinates[1][i])
        {
            return Err(JkrTorsionalPathError::UnresolvedEnvelope);
        }
        let low = coordinates[0][0];
        let high = coordinates[1][0];
        let p = [coordinates[0][1], coordinates[1][1]];
        let maximum_twist = p[0].abs().max(p[1].abs());
        let minimum_twist = if p[0].is_sign_negative() == p[1].is_sign_negative() {
            p[0].abs().min(p[1].abs())
        } else {
            0.0
        };
        let mut evaluate =
            |delta: f64, twist: f64| -> Result<JkrTorsionalResponse, JkrTorsionalPathError> {
                poll(gate)?;
                for cached in [before, contact] {
                    if cached.indentation_m == delta && cached.twist_angle_rad.abs() == twist {
                        return Ok(cached);
                    }
                }
                charge(evaluations, maximum_evaluations)?;
                let point = self
                    .evaluate_torsion_displacement_controlled(delta, twist)
                    .map_err(|_| JkrTorsionalPathError::UnresolvedEnvelope)?;
                strictly_attached(point).map_err(|_| JkrTorsionalPathError::UnresolvedEnvelope)?;
                Ok(point)
            };
        // The ideal twist limit increases on delta<0 and decreases on delta>0.
        // Both extreme approaches at maximum twist must therefore attach;
        // testing only the minimum-radius corner would miss the opening edge.
        let lower = evaluate(low, maximum_twist)?;
        evaluate(high, maximum_twist)?;
        let upper = evaluate(high, minimum_twist)?;
        // On this stable opening branch a_delta>=0 and a_|phi|<=0. These
        // corners enclose ALL radii, including expansion on twist reversal.
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
            return Err(JkrTorsionalPathError::UnresolvedEnvelope);
        }
        Ok(bounds)
    }
}

fn strictly_attached(point: JkrTorsionalResponse) -> Result<[[f64; 2]; 2], JkrTorsionalPathError> {
    if point.boundary.is_some() || point.displacement_radius_curvature_n_per_m <= 0.0 {
        return Err(JkrTorsionalPathError::UnresolvedSecant);
    }
    point
        .coupled_stiffness
        .ok_or(JkrTorsionalPathError::UnresolvedSecant)
}
fn poll(gate: &CancelGate) -> Result<(), JkrTorsionalPathError> {
    if gate.is_requested() {
        Err(JkrTorsionalPathError::Cancelled)
    } else {
        Ok(())
    }
}
fn charge(evaluations: &mut usize, maximum: usize) -> Result<(), JkrTorsionalPathError> {
    if *evaluations >= maximum {
        return Err(JkrTorsionalPathError::Budget);
    }
    *evaluations += 1;
    Ok(())
}
fn dot(a: [f64; 2], b: [f64; 2]) -> f64 {
    a[0] * b[0] + a[1] * b[1]
}
