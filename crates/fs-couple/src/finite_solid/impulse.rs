//! Instantaneous momentum transfer on the existing finite-strain solid.
//!
//! A P1 material point uses the SAME barycentric map for velocity sampling and
//! transpose impulse distribution. Thus an impulse J acts at the actual current
//! point x, supplies angular impulse (x-origin) cross J, and does work through
//! the average point velocity. The original consistent reference mass supplies
//! M_ff delta-v_f = J_f; a diagonal/lumped surrogate would change the response.
//!
//! This is a fixed-configuration impulse, not a collision detector, restitution
//! or friction law, contact complementarity solver, or time step. The caller
//! owns the impulse law and any equal/opposite partner. Stationary supports
//! retain their reaction impulse. Failed solves/balance gates/cancellation leave
//! the accepted body untouched. The conservative velocity-independent material
//! force is unchanged, so q, a, physical time and step history are retained.

use super::{
    Csr, Cx, FiniteSolidBody, FiniteSolidError, FiniteSolidInitialization, MassDiagonal,
    SecondOrderProblem, SolveReport, checkpoint, inflate, second_order_acceleration,
};

/// A material point inside one of the original P1 tetrahedra.
/// Zero weights select faces, edges or vertices without a separate load map.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteSolidPoint {
    tetrahedron: usize,
    weights: [f64; 4],
}

impl FiniteSolidPoint {
    /// Admit nonnegative barycentric weights summing to one within roundoff.
    /// The roundoff-sized sum defect is normalized once. Extrapolation refuses.
    /// The tetrahedron index is checked against the body at each use.
    pub fn new(tetrahedron: usize, weights: [f64; 4]) -> Result<Self, FiniteSolidError> {
        let sum: f64 = weights.iter().sum();
        if weights
            .iter()
            .any(|w| !w.is_finite() || !(0.0..=1.0).contains(w))
            || !sum.is_finite()
            || (sum - 1.0).abs() > 64.0 * f64::EPSILON
        {
            return Err(invalid(
                "point requires nonnegative unit-sum barycentric weights",
            ));
        }
        Ok(Self {
            tetrahedron,
            weights: weights.map(|w| w / sum),
        })
    }

    pub const fn tetrahedron(self) -> usize {
        self.tetrahedron
    }
    pub const fn weights(self) -> [f64; 4] {
        self.weights
    }
}

/// Continuum momentum integrated with the original consistent P1 mass.
/// Angular momentum uses CURRENT deformed coordinates, about the stated origin.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteSolidMomentum {
    pub origin_m: [f64; 3],
    pub linear_kg_m_s: [f64; 3],
    pub angular_kg_m2_s: [f64; 3],
}

/// Bounded native mass solve and absolute physical impulse/work tolerances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteSolidImpulseControls {
    pub mass_solve: FiniteSolidInitialization,
    pub linear_tolerance_kg_m_s: f64,
    pub angular_tolerance_kg_m2_s: f64,
    pub work_tolerance_j: f64,
}

impl Default for FiniteSolidImpulseControls {
    fn default() -> Self {
        Self {
            mass_solve: FiniteSolidInitialization::default(),
            linear_tolerance_kg_m_s: 1e-10,
            angular_tolerance_kg_m2_s: 1e-10,
            work_tolerance_j: 1e-10,
        }
    }
}

/// Accepted physical change; no numerical defect is relabeled as heat.
#[derive(Clone, Debug)]
pub struct FiniteSolidImpulse {
    pub solve: SolveReport,
    pub before: FiniteSolidMomentum,
    pub after: FiniteSolidMomentum,
    pub applied: FiniteSolidMomentum,
    pub support: FiniteSolidMomentum,
    /// Impulse ON the solid, in original prescribed-DOF order [kg m/s].
    pub support_impulses_kg_m_s: Vec<(usize, f64)>,
    pub linear_residual_kg_m_s: [f64; 3],
    pub angular_residual_kg_m2_s: [f64; 3],
    pub kinetic_before_j: f64,
    pub kinetic_after_j: f64,
    /// Computed as average(v) dot M delta-v, without subtracting large energies.
    pub kinetic_change_j: f64,
    /// Applied impulse dotted with average velocity; stationary supports do no work.
    pub impulse_work_j: f64,
    pub work_residual_j: f64,
}

impl<'a> FiniteSolidBody<'a> {
    /// Observe the actual linear and angular momentum at an arbitrary world origin.
    pub fn momentum_about(
        &self,
        origin_m: [f64; 3],
        cx: &Cx<'_>,
    ) -> Result<FiniteSolidMomentum, FiniteSolidError> {
        checkpoint(cx)?;
        let positions = self.current_positions(cx)?;
        let velocity: Vec<f64> = self.velocity_m_s().into_iter().flatten().collect();
        let mut momentum = vec![0.0; velocity.len()];
        self.mass.spmv(&velocity, &mut momentum);
        moments(&positions, &momentum, origin_m, cx)
    }

    pub fn point_position_m(
        &self,
        point: FiniteSolidPoint,
        cx: &Cx<'_>,
    ) -> Result<[f64; 3], FiniteSolidError> {
        checkpoint(cx)?;
        let nodes = self.point_nodes(point)?;
        let mut result = [0.0; 3];
        for (&node, weight) in nodes.iter().zip(point.weights) {
            for (c, value) in result.iter_mut().enumerate() {
                let dof = 3 * node + c;
                let q = self
                    .free
                    .binary_search(&dof)
                    .map_or(self.fixed[dof], |row| self.state.q[row]);
                *value = weight.mul_add(self.problem.nodes_m[node][c] + q, *value);
            }
        }
        finite(&result)?;
        Ok(result)
    }

    /// The velocity trace paired with `point_impulse_candidate` by virtual work.
    pub fn point_velocity_m_s(
        &self,
        point: FiniteSolidPoint,
        cx: &Cx<'_>,
    ) -> Result<[f64; 3], FiniteSolidError> {
        checkpoint(cx)?;
        let nodes = self.point_nodes(point)?;
        let mut result = [0.0; 3];
        for (&node, weight) in nodes.iter().zip(point.weights) {
            for (c, value) in result.iter_mut().enumerate() {
                let v = self
                    .free
                    .binary_search(&(3 * node + c))
                    .map_or(0.0, |row| self.state.v[row]);
                *value = weight.mul_add(v, *value);
            }
        }
        finite(&result)?;
        Ok(result)
    }

    /// Directional point compliance d^T B M_ff^-1 B^T d [1/kg].
    /// Use a unit direction for contact effective mass; rescaling d rescales
    /// this value quadratically. A fully fixed direction returns zero.
    /// The same native mass solve and point map drive the eventual impulse.
    pub fn point_directional_inverse_mass_per_kg(
        &self,
        point: FiniteSolidPoint,
        direction_world: [f64; 3],
        controls: FiniteSolidInitialization,
        cx: &Cx<'_>,
    ) -> Result<f64, FiniteSolidError> {
        checkpoint(cx)?;
        finite(&direction_world)?;
        if magnitude(direction_world) == 0.0 {
            return Err(invalid("effective mass requires a nonzero direction"));
        }
        let impulse = self.point_distribution(point, direction_world)?;
        let (change, _) = self.impulse_response(&impulse, controls, cx)?;
        let value = self
            .free
            .iter()
            .zip(change)
            .fold(0.0, |s, (&dof, dv)| impulse[dof].mul_add(dv, s));
        if !value.is_finite() || value < 0.0 {
            return Err(invalid("point inverse mass is not finite and nonnegative"));
        }
        checkpoint(cx)?;
        Ok(value)
    }

    /// Apply a world impulse at the CURRENT barycentrically interpolated point.
    /// Its transpose distribution preserves resultant, moment and virtual work.
    pub fn point_impulse_candidate(
        &self,
        point: FiniteSolidPoint,
        impulse_world_kg_m_s: [f64; 3],
        origin_m: [f64; 3],
        controls: FiniteSolidImpulseControls,
        cx: &Cx<'_>,
    ) -> Result<(Self, FiniteSolidImpulse), FiniteSolidError> {
        checkpoint(cx)?;
        let impulse = self.point_distribution(point, impulse_world_kg_m_s)?;
        self.apply_impulse(&impulse, origin_m, controls, cx)
    }

    /// Apply nodal world impulses; repeated node entries accumulate in input order.
    /// Impulses on supported DOFs are retained in the support account.
    pub fn nodal_impulse_candidate(
        &self,
        impulses_kg_m_s: &[(usize, [f64; 3])],
        origin_m: [f64; 3],
        controls: FiniteSolidImpulseControls,
        cx: &Cx<'_>,
    ) -> Result<(Self, FiniteSolidImpulse), FiniteSolidError> {
        checkpoint(cx)?;
        let mut impulse = vec![0.0; self.mass.nrows()];
        for &(node, value) in impulses_kg_m_s {
            checkpoint(cx)?;
            finite(&value)?;
            if node >= self.problem.nodes_m.len() {
                return Err(invalid("impulse node is outside the original mesh"));
            }
            for c in 0..3 {
                impulse[3 * node + c] += value[c];
            }
        }
        self.apply_impulse(&impulse, origin_m, controls, cx)
    }

    fn point_nodes(&self, point: FiniteSolidPoint) -> Result<[usize; 4], FiniteSolidError> {
        self.problem
            .tetrahedra
            .get(point.tetrahedron)
            .copied()
            .ok_or_else(|| invalid("impulse point tetrahedron is outside the original mesh"))
    }

    fn point_distribution(
        &self,
        point: FiniteSolidPoint,
        value: [f64; 3],
    ) -> Result<Vec<f64>, FiniteSolidError> {
        finite(&value)?;
        let nodes = self.point_nodes(point)?;
        let mut impulse = vec![0.0; self.mass.nrows()];
        for (&node, weight) in nodes.iter().zip(point.weights) {
            for c in 0..3 {
                impulse[3 * node + c] += weight * value[c];
            }
        }
        finite(&impulse)?;
        Ok(impulse)
    }

    fn current_positions(&self, cx: &Cx<'_>) -> Result<Vec<[f64; 3]>, FiniteSolidError> {
        let u = inflate(&self.state.q, &self.free, &self.fixed);
        let mut positions = self.problem.nodes_m.to_vec();
        for (node, position) in positions.iter_mut().enumerate() {
            if node % 256 == 0 {
                checkpoint(cx)?;
            }
            for c in 0..3 {
                position[c] += u[3 * node + c];
            }
            finite(position)?;
        }
        Ok(positions)
    }

    fn impulse_response(
        &self,
        impulse: &[f64],
        controls: FiniteSolidInitialization,
        cx: &Cx<'_>,
    ) -> Result<(Vec<f64>, SolveReport), FiniteSolidError> {
        checkpoint(cx)?;
        finite(impulse)?;
        let forcing: Vec<f64> = self.free.iter().map(|&dof| impulse[dof]).collect();
        let zero = vec![0.0; self.free.len()];
        // The native initializer is a bounded M x = rhs solver. With zero
        // resistance, its unknown is a velocity jump rather than acceleration;
        // no time interval, constitutive force or new solver enters this map.
        let response = second_order_acceleration(
            &ImpulseMass {
                mass: &self.free_mass,
            },
            &zero,
            &zero,
            &forcing,
            &MassDiagonal {
                mass: &self.free_mass,
            },
            controls.solve,
            controls.maximum_workspace_components,
            &mut || cx.checkpoint().is_err(),
        )?;
        Ok((response.value, response.report))
    }

    fn apply_impulse(
        &self,
        impulse: &[f64],
        origin: [f64; 3],
        controls: FiniteSolidImpulseControls,
        cx: &Cx<'_>,
    ) -> Result<(Self, FiniteSolidImpulse), FiniteSolidError> {
        finite(&origin)?;
        if [
            controls.linear_tolerance_kg_m_s,
            controls.angular_tolerance_kg_m2_s,
            controls.work_tolerance_j,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
        {
            return Err(invalid(
                "impulse balance tolerances must be finite and positive",
            ));
        }
        let (change, solve) = self.impulse_response(impulse, controls.mass_solve, cx)?;
        let mut candidate = self.clone();
        let mut dv = vec![0.0; impulse.len()];
        let mut average_v = vec![0.0; impulse.len()];
        for (row, (&dof, change)) in self.free.iter().zip(change).enumerate() {
            if row % 256 == 0 {
                checkpoint(cx)?;
            }
            let old = self.state.v[row];
            let new = old + change;
            candidate.state.v[row] = new;
            // Account for the actual representable change, not the ideal solve.
            dv[dof] = new - old;
            average_v[dof] = 0.5 * old + 0.5 * new;
        }
        finite(&candidate.state.v)?;
        let mut delta_p = vec![0.0; impulse.len()];
        self.mass.spmv(&dv, &mut delta_p);
        finite(&delta_p)?;
        let mut support = vec![0.0; impulse.len()];
        let mut support_impulses_kg_m_s = Vec::with_capacity(self.problem.prescribed_m.len());
        for &(dof, _) in self.problem.prescribed_m {
            support[dof] = delta_p[dof] - impulse[dof];
            support_impulses_kg_m_s.push((dof, support[dof]));
        }
        finite(&support)?;
        let positions = self.current_positions(cx)?;
        let applied = moments(&positions, impulse, origin, cx)?;
        let support_moments = moments(&positions, &support, origin, cx)?;
        let residual: Vec<f64> = delta_p
            .iter()
            .zip(impulse)
            .zip(&support)
            .map(|((&actual, &applied), &support)| actual - applied - support)
            .collect();
        let defect = moments(&positions, &residual, origin, cx)?;
        let kinetic_change_j = super::dot(&average_v, &delta_p);
        let impulse_work_j = super::dot(&average_v, impulse);
        let work_residual_j = kinetic_change_j - impulse_work_j;
        finite(&[kinetic_change_j, impulse_work_j, work_residual_j])?;
        if magnitude(defect.linear_kg_m_s) > controls.linear_tolerance_kg_m_s
            || magnitude(defect.angular_kg_m2_s) > controls.angular_tolerance_kg_m2_s
            || work_residual_j.abs() > controls.work_tolerance_j
        {
            return Err(FiniteSolidError::ImpulseBalance {
                linear_residual_kg_m_s: defect.linear_kg_m_s,
                angular_residual_kg_m2_s: defect.angular_kg_m2_s,
                work_residual_j,
            });
        }
        let report = FiniteSolidImpulse {
            solve,
            before: self.momentum_about(origin, cx)?,
            after: candidate.momentum_about(origin, cx)?,
            applied,
            support: support_moments,
            support_impulses_kg_m_s,
            linear_residual_kg_m_s: defect.linear_kg_m_s,
            angular_residual_kg_m2_s: defect.angular_kg_m2_s,
            kinetic_before_j: self.kinetic_energy_j()?,
            kinetic_after_j: candidate.kinetic_energy_j()?,
            kinetic_change_j,
            impulse_work_j,
            work_residual_j,
        };
        checkpoint(cx)?;
        Ok((candidate, report))
    }
}

struct ImpulseMass<'a> {
    mass: &'a Csr,
}
impl SecondOrderProblem for ImpulseMass<'_> {
    fn dimension(&self) -> usize {
        self.mass.nrows()
    }
    fn mass_apply(&self, input: &[f64], output: &mut [f64]) {
        self.mass.spmv(input, output);
    }
    fn damping_apply(&self, _input: &[f64], output: &mut [f64]) {
        output.fill(0.0);
    }
    fn internal_force(&self, _q: &[f64], output: &mut [f64]) {
        output.fill(0.0);
    }
    fn tangent_apply(&self, _q: &[f64], _d: &[f64], output: &mut [f64]) {
        output.fill(0.0);
    }
    fn resistance(&self, _q: &[f64], _v: &[f64], output: &mut [f64]) {
        output.fill(0.0);
    }
}

fn invalid(message: &'static str) -> FiniteSolidError {
    FiniteSolidError::Invalid(message)
}
fn finite(values: &[f64]) -> Result<(), FiniteSolidError> {
    if values.iter().all(|v| v.is_finite()) {
        Ok(())
    } else {
        Err(invalid(
            "solid impulse, momentum or work is unrepresentable",
        ))
    }
}
fn magnitude(v: [f64; 3]) -> f64 {
    v[0].hypot(v[1]).hypot(v[2])
}

fn moments(
    positions: &[[f64; 3]],
    momentum: &[f64],
    origin_m: [f64; 3],
    cx: &Cx<'_>,
) -> Result<FiniteSolidMomentum, FiniteSolidError> {
    finite(&origin_m)?;
    let mut linear = [0.0; 3];
    let mut angular = [0.0; 3];
    for (node, position) in positions.iter().enumerate() {
        if node % 256 == 0 {
            checkpoint(cx)?;
        }
        let r: [f64; 3] = std::array::from_fn(|c| position[c] - origin_m[c]);
        let p = &momentum[3 * node..3 * node + 3];
        finite(&r)?;
        finite(p)?;
        for c in 0..3 {
            linear[c] += p[c];
            let j = (c + 1) % 3;
            let k = (c + 2) % 3;
            angular[c] += r[j].mul_add(p[k], -r[k] * p[j]);
        }
    }
    finite(&linear)?;
    finite(&angular)?;
    Ok(FiniteSolidMomentum {
        origin_m,
        linear_kg_m_s: linear,
        angular_kg_m2_s: angular,
    })
}
