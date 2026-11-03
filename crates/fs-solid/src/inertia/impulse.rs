//! Instantaneous momentum transfer through the original consistent mass.
//!
//! The constitutive owner keeps position, elapsed time and every material
//! history variable. This owner changes velocity only, at the same geometry,
//! using the native checked CG solver on M_ff. Held components keep their
//! current velocity; their reaction impulse and signed work remain explicit.

use fs_math::dd::Dd;
use fs_solver::{CheckedCgConfig, CheckedCgError, LinearOp, SolveReport, checked_cg};
use fs_sparse::precond::Precond;

use super::{Cx, InertiaError, ReferenceMass, checkpoint, dot, finite, invalid};

/// Explicit native solve, storage and physical balance allowances.
#[derive(Clone, Copy, Debug)]
pub struct SolidImpulseControls {
    pub solve: CheckedCgConfig,
    /// Conservative bound on simultaneous scalar/index scratch and output
    /// components, including native CG state and residual histories.
    pub maximum_workspace_components: usize,
    /// Bounds both the total linear residual norm and every free-DOF impulse
    /// residual after rounding the proposed velocity [kg m/s].
    pub linear_tolerance_kg_m_s: f64,
    pub angular_tolerance_kg_m2_s: f64,
    pub work_tolerance_j: f64,
}

/// Momentum integrated with the original mass at current nodal positions.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SolidMomentum {
    pub origin_m: [f64; 3],
    pub linear_kg_m_s: [f64; 3],
    pub angular_kg_m2_s: [f64; 3],
}

/// Physical impulse accounting for an admitted, fixed-position candidate.
#[derive(Clone, Debug)]
pub struct SolidImpulseReport {
    /// Original native solver report, with recomputed Euclidean residual.
    pub solve: SolveReport,
    pub before: SolidMomentum,
    pub after: SolidMomentum,
    pub applied: SolidMomentum,
    /// Reaction impulse ON the solid from its held velocity components.
    pub support: SolidMomentum,
    /// Original held-DOF order; includes zero reactions [kg m/s].
    pub support_impulses_kg_m_s: Vec<(usize, f64)>,
    /// Derived from M times the actual representable velocity change, without
    /// subtracting potentially much larger before/after momenta.
    pub linear_residual_kg_m_s: [f64; 3],
    pub angular_residual_kg_m2_s: [f64; 3],
    pub maximum_free_impulse_residual_kg_m_s: f64,
    pub kinetic_before_j: f64,
    pub kinetic_after_j: f64,
    /// Average velocity dotted with the original M times actual delta-v.
    pub kinetic_change_j: f64,
    /// Applied impulse dotted with average velocity.
    pub applied_work_j: f64,
    /// Signed reaction work; a held nonzero velocity can supply or remove work.
    pub support_work_j: f64,
    /// Average velocity dotted with the independently recomputed impulse
    /// residual, avoiding subtraction of large work totals. Never heat.
    pub work_residual_j: f64,
}

impl ReferenceMass {
    /// Preserve every accepted state component except instantaneous velocity.
    /// Repeated applied DOFs accumulate with compensation in caller order and
    /// round once per resultant; repeated held DOFs refuse.
    /// Held DOFs retain the accepted velocity, rather than being set to zero.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn impulse_candidate(
        &self,
        positions_m: &[[f64; 3]],
        velocity_mps: &[[f64; 3]],
        held_velocity_dofs: &[usize],
        nodal_impulses_kg_m_s: &[(usize, f64)],
        origin_m: [f64; 3],
        controls: SolidImpulseControls,
        cx: &Cx<'_>,
    ) -> Result<(Vec<[f64; 3]>, SolidImpulseReport), InertiaError> {
        checkpoint(cx)?;
        validate_motion(self.ndofs(), positions_m, velocity_mps, origin_m, cx)?;
        let (mut candidate, solve) =
            self.free_velocity_response(held_velocity_dofs, nodal_impulses_kg_m_s, controls, cx)?;
        for (node, (next, before)) in candidate.iter_mut().zip(velocity_mps).enumerate() {
            if node % 256 == 0 {
                checkpoint(cx)?;
            }
            for axis in 0..3 {
                next[axis] = if next[axis] == 0.0 {
                    before[axis]
                } else {
                    finite(before[axis] + next[axis])?
                };
            }
        }
        self.finish_impulse_candidate(
            positions_m,
            velocity_mps,
            held_velocity_dofs,
            nodal_impulses_kg_m_s,
            origin_m,
            controls,
            candidate,
            solve,
            cx,
        )
    }

    /// Supported velocity response without touching geometry or material history.
    /// The original sparse M_ff is borrowed, never lumped or explicitly inverted.
    pub(crate) fn free_velocity_response(
        &self,
        held_velocity_dofs: &[usize],
        nodal_impulses_kg_m_s: &[(usize, f64)],
        controls: SolidImpulseControls,
        cx: &Cx<'_>,
    ) -> Result<(Vec<[f64; 3]>, SolveReport), InertiaError> {
        checkpoint(cx)?;
        let n = self.ndofs();
        let (held, applied) =
            admit_impulses(n, held_velocity_dofs, nodal_impulses_kg_m_s, controls, cx)?;
        let mut free = Vec::with_capacity(n - held_velocity_dofs.len());
        let mut reduced = vec![usize::MAX; n];
        let mut diagonal = Vec::with_capacity(free.capacity());
        let mut rhs = Vec::with_capacity(free.capacity());
        for dof in 0..n {
            if dof % 256 == 0 {
                checkpoint(cx)?;
            }
            if !held[dof] {
                reduced[dof] = free.len();
                free.push(dof);
                diagonal.push(self.matrix.get(dof, dof));
                rhs.push(applied[dof]);
            }
        }
        let operator = FreeMass {
            mass: self,
            free: &free,
            reduced: &reduced,
            cx,
        };
        let solved = checked_cg(
            &operator,
            &MassDiagonal {
                values: diagonal,
                cx,
            },
            &rhs,
            controls.solve,
            |_| checkpoint(cx),
        )
        .map_err(|error| match error {
            CheckedCgError::InvalidInput(what) => invalid(what),
            CheckedCgError::Interrupted(error) => error,
        })?;
        if !solved.report.converged_euclidean() {
            return Err(invalid(&format!(
                "native impulse mass solve did not converge: {:?}",
                solved.report
            )));
        }
        let mut response = vec![[0.0; 3]; n / 3];
        for (row, (&dof, &change)) in free.iter().zip(&solved.x).enumerate() {
            if row % 256 == 0 {
                checkpoint(cx)?;
            }
            response[dof / 3][dof % 3] = finite(change)?;
        }
        checkpoint(cx)?;
        Ok((response, solved.report))
    }

    /// Audit the actual represented result of an existing native constrained
    /// solve. Pressure or other internal impulses must be included in the
    /// effective nodal right-hand side, with their separate physical account
    /// retained by that constraint owner. This method performs no second solve.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn finish_impulse_candidate(
        &self,
        positions_m: &[[f64; 3]],
        velocity_mps: &[[f64; 3]],
        held_velocity_dofs: &[usize],
        nodal_effective_impulses: &[(usize, f64)],
        origin_m: [f64; 3],
        controls: SolidImpulseControls,
        candidate: Vec<[f64; 3]>,
        solve: SolveReport,
        cx: &Cx<'_>,
    ) -> Result<(Vec<[f64; 3]>, SolidImpulseReport), InertiaError> {
        checkpoint(cx)?;
        let n = self.ndofs();
        validate_motion(n, positions_m, velocity_mps, origin_m, cx)?;
        validate_motion(n, positions_m, &candidate, origin_m, cx)?;
        if !solve.converged_euclidean() {
            return Err(invalid("impulse requires a converged native solve report"));
        }
        let (held, applied) = admit_impulses(
            n,
            held_velocity_dofs,
            nodal_effective_impulses,
            controls,
            cx,
        )?;
        let before_v = velocity_mps.as_flattened();
        let after_v = candidate.as_flattened();
        for (row, &dof) in held_velocity_dofs.iter().enumerate() {
            if row % 256 == 0 {
                checkpoint(cx)?;
            }
            if before_v[dof].to_bits() != after_v[dof].to_bits() {
                return Err(invalid("an impulse must preserve held velocities exactly"));
            }
        }
        let mut delta_v = vec![0.0; n];
        let mut average_v = vec![0.0; n];
        for dof in 0..n {
            if dof % 256 == 0 {
                checkpoint(cx)?;
            }
            // Gate what the returned velocity can actually represent, not the
            // ideal Krylov increment that may disappear at large velocities.
            delta_v[dof] = finite(after_v[dof] - before_v[dof])?;
            average_v[dof] = 0.5 * before_v[dof] + 0.5 * after_v[dof];
        }
        let delta_p = self.action(&delta_v, cx)?;
        let mut support = vec![0.0; n];
        let mut support_impulses_kg_m_s = Vec::with_capacity(held_velocity_dofs.len());
        for (row, &dof) in held_velocity_dofs.iter().enumerate() {
            if row % 256 == 0 {
                checkpoint(cx)?;
            }
            support[dof] = finite(delta_p[dof] - applied[dof])?;
            support_impulses_kg_m_s.push((dof, support[dof]));
        }
        let mut residual = vec![0.0; n];
        let mut maximum_free_impulse_residual_kg_m_s = 0.0_f64;
        for dof in 0..n {
            if dof % 256 == 0 {
                checkpoint(cx)?;
            }
            residual[dof] = finite((delta_p[dof] - applied[dof]) - support[dof])?;
            if !held[dof] {
                maximum_free_impulse_residual_kg_m_s =
                    maximum_free_impulse_residual_kg_m_s.max(residual[dof].abs());
            }
        }
        let defect = moments(positions_m, &residual, origin_m, cx)?;
        let kinetic_change_j = finite(dot(&average_v, &delta_p))?;
        let applied_work_j = finite(dot(&average_v, &applied))?;
        let support_work_j = finite(dot(&average_v, &support))?;
        let work_residual_j = finite(dot(&average_v, &residual))?;
        if magnitude(defect.linear_kg_m_s) > controls.linear_tolerance_kg_m_s
            || magnitude(defect.angular_kg_m2_s) > controls.angular_tolerance_kg_m2_s
            || maximum_free_impulse_residual_kg_m_s > controls.linear_tolerance_kg_m_s
            || work_residual_j.abs() > controls.work_tolerance_j
        {
            return Err(invalid(&format!(
                "represented impulse balance exceeds allowance: linear {:?}, angular {:?}, free {}, work {}",
                defect.linear_kg_m_s,
                defect.angular_kg_m2_s,
                maximum_free_impulse_residual_kg_m_s,
                work_residual_j,
            )));
        }
        let before_p = self.action(before_v, cx)?;
        let after_p = self.action(after_v, cx)?;
        let report = SolidImpulseReport {
            solve,
            before: moments(positions_m, &before_p, origin_m, cx)?,
            after: moments(positions_m, &after_p, origin_m, cx)?,
            applied: moments(positions_m, &applied, origin_m, cx)?,
            support: moments(positions_m, &support, origin_m, cx)?,
            support_impulses_kg_m_s,
            linear_residual_kg_m_s: defect.linear_kg_m_s,
            angular_residual_kg_m2_s: defect.angular_kg_m2_s,
            maximum_free_impulse_residual_kg_m_s,
            kinetic_before_j: kinetic(before_v, &before_p)?,
            kinetic_after_j: kinetic(after_v, &after_p)?,
            kinetic_change_j,
            applied_work_j,
            support_work_j,
            work_residual_j,
        };
        checkpoint(cx)?;
        Ok((candidate, report))
    }
}

fn validate_motion(
    n: usize,
    positions_m: &[[f64; 3]],
    velocity_mps: &[[f64; 3]],
    origin_m: [f64; 3],
    cx: &Cx<'_>,
) -> Result<(), InertiaError> {
    if n % 3 != 0
        || positions_m.len() != n / 3
        || velocity_mps.len() != n / 3
        || origin_m.iter().any(|x| !x.is_finite())
    {
        return Err(invalid(
            "complete finite positions/velocities and finite origin required",
        ));
    }
    for (node, (position, velocity)) in positions_m.iter().zip(velocity_mps).enumerate() {
        if node % 256 == 0 {
            checkpoint(cx)?;
        }
        if position.iter().chain(velocity).any(|x| !x.is_finite()) {
            return Err(invalid("finite positions and velocities required"));
        }
    }
    Ok(())
}

fn admit_impulses(
    n: usize,
    held_velocity_dofs: &[usize],
    nodal_impulses_kg_m_s: &[(usize, f64)],
    controls: SolidImpulseControls,
    cx: &Cx<'_>,
) -> Result<(Vec<bool>, Vec<f64>), InertiaError> {
    if n % 3 != 0
        || [
            controls.linear_tolerance_kg_m_s,
            controls.angular_tolerance_kg_m2_s,
            controls.work_tolerance_j,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
        || !controls.solve.tolerance.is_finite()
        || controls.solve.tolerance <= 0.0
        || controls.solve.tolerance >= 1.0
    {
        return Err(invalid(
            "positive finite impulse tolerances and a relative solve tolerance below one required",
        ));
    }
    // Include native solver histories as well as fixed vectors. Both the
    // response solve and post-solve account fit this conservative common cap.
    let required = n
        .checked_mul(48)
        .and_then(|vectors| {
            controls
                .solve
                .max_iterations
                .checked_add(1)
                .and_then(|history| history.checked_mul(4))
                .and_then(|history| vectors.checked_add(history))
        })
        .ok_or_else(|| invalid("impulse workspace size overflow"))?;
    if required > controls.maximum_workspace_components {
        return Err(invalid("impulse exceeds its workspace component cap"));
    }
    let mut held = vec![false; n];
    for (index, &dof) in held_velocity_dofs.iter().enumerate() {
        if index % 256 == 0 {
            checkpoint(cx)?;
        }
        if dof >= n || held[dof] {
            return Err(invalid("distinct in-range held velocity DOFs required"));
        }
        held[dof] = true;
    }
    // These entries are simultaneous loads. Preserve a small net impulse
    // between large opposing contributions before the native mass solve.
    let mut accumulated = vec![Dd::ZERO; n];
    for (index, &(dof, impulse)) in nodal_impulses_kg_m_s.iter().enumerate() {
        if index % 256 == 0 {
            checkpoint(cx)?;
        }
        if dof >= n || !impulse.is_finite() {
            return Err(invalid("finite in-range nodal impulses required"));
        }
        accumulated[dof] = accumulated[dof] + Dd::from_f64(impulse);
        finite(accumulated[dof].hi)?;
        finite(accumulated[dof].lo)?;
    }
    let mut applied = Vec::with_capacity(n);
    for (dof, resultant) in accumulated.into_iter().enumerate() {
        if dof % 256 == 0 {
            checkpoint(cx)?;
        }
        applied.push(finite(resultant.to_f64())?);
    }
    Ok((held, applied))
}

struct FreeMass<'a, 'cx> {
    mass: &'a ReferenceMass,
    free: &'a [usize],
    reduced: &'a [usize],
    cx: &'a Cx<'cx>,
}

impl LinearOp for FreeMass<'_, '_> {
    fn n(&self) -> usize {
        self.free.len()
    }

    fn apply(&self, input: &[f64], output: &mut [f64]) {
        for (row, &dof) in self.free.iter().enumerate() {
            if row % 256 == 0 && self.cx.checkpoint().is_err() {
                // checked_cg's next checkpoint preserves cancellation rather
                // than accepting this unfinished operator application.
                output.fill(f64::NAN);
                return;
            }
            let (columns, values) = self.mass.matrix.row(dof);
            let mut value = 0.0;
            for (&column, &entry) in columns.iter().zip(values) {
                let reduced = self.reduced[column];
                if reduced != usize::MAX {
                    value = entry.mul_add(input[reduced], value);
                }
            }
            output[row] = value;
        }
    }
}

struct MassDiagonal<'a, 'cx> {
    values: Vec<f64>,
    cx: &'a Cx<'cx>,
}
impl Precond for MassDiagonal<'_, '_> {
    fn apply(&self, input: &[f64], output: &mut [f64]) {
        for (row, (&rhs, &mass)) in input.iter().zip(&self.values).enumerate() {
            if row % 256 == 0 && self.cx.checkpoint().is_err() {
                output.fill(f64::NAN);
                return;
            }
            output[row] = rhs / mass;
        }
    }
}

fn kinetic(velocity: &[f64], momentum: &[f64]) -> Result<f64, InertiaError> {
    let energy = finite(
        velocity
            .iter()
            .zip(momentum)
            .fold(0.0, |sum, (&v, &p)| (0.5 * v).mul_add(p, sum)),
    )?;
    if energy < 0.0 {
        return Err(invalid("negative consistent-mass kinetic energy"));
    }
    Ok(energy)
}

fn magnitude(vector: [f64; 3]) -> f64 {
    vector[0].hypot(vector[1]).hypot(vector[2])
}

pub(crate) fn moments(
    positions: &[[f64; 3]],
    momentum: &[f64],
    origin_m: [f64; 3],
    cx: &Cx<'_>,
) -> Result<SolidMomentum, InertiaError> {
    if positions.len().checked_mul(3) != Some(momentum.len())
        || origin_m.iter().any(|value| !value.is_finite())
    {
        return Err(invalid(
            "complete momentum field and finite moment origin required",
        ));
    }
    let mut linear = [0.0; 3];
    let mut angular = [0.0; 3];
    for (node, position) in positions.iter().enumerate() {
        if node % 256 == 0 {
            checkpoint(cx)?;
        }
        let r: [f64; 3] = std::array::from_fn(|c| position[c] - origin_m[c]);
        if r.iter().any(|value| !value.is_finite()) {
            return Err(invalid("unrepresentable impulse moment arm"));
        }
        let p = &momentum[3 * node..3 * node + 3];
        for c in 0..3 {
            linear[c] += p[c];
            let j = (c + 1) % 3;
            let k = (c + 2) % 3;
            angular[c] += r[j].mul_add(p[k], -r[k] * p[j]);
        }
    }
    for value in linear.iter().chain(&angular) {
        finite(*value)?;
    }
    Ok(SolidMomentum {
        origin_m,
        linear_kg_m_s: linear,
        angular_kg_m2_s: angular,
    })
}
