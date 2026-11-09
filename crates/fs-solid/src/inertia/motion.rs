//! Original-mass momentum over native material and kinetic drift intervals.
use fs_math::dd::Dd;

use super::{Cx, InertiaError, ReferenceMass, SolidMomentum, checkpoint, finite, invalid};

/// Actual physical increments and endpoint loads about one declared origin.
/// Backward Euler's angular drift is a numerical defect, never a torque source.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SolidMotionMomentum {
    pub change: SolidMomentum,
    pub applied: SolidMomentum,
    pub support: SolidMomentum,
    /// Sum(delta_u cross M v_before), including stored-displacement rounding. For the
    /// exact BE displacement relation and symmetric mass this equals
    /// -dt sum(v_after cross M delta_v). It is not physical angular impulse.
    pub angular_integration_defect_kg_m2_s: [f64; 3],
    pub linear_residual_kg_m_s: [f64; 3],
    /// Actual angular change minus applied/support impulse and the separately
    /// reported numerical defect. Consumers must also bound that defect itself.
    pub angular_residual_kg_m2_s: [f64; 3],
}

/// Native constant-velocity phase: momentum inventories and stable increments
/// share one compensated sparse mass action. The same mass and stored velocity
/// imply an exactly zero kinetic change; stored-displacement rounding can leave
/// an angular error.
pub(crate) struct ConstantVelocityDrift {
    pub before: SolidMomentum,
    pub after: SolidMomentum,
    pub change: SolidMomentum,
    pub kinetic_j: f64,
}

impl ReferenceMass {
    pub(crate) fn constant_velocity_drift(
        &self,
        reference_nodes_m: &[[f64; 3]],
        displacement_before_m: &[[f64; 3]],
        displacement_after_m: &[[f64; 3]],
        velocity_mps: &[[f64; 3]],
        origin_m: [f64; 3],
        cx: &Cx<'_>,
    ) -> Result<ConstantVelocityDrift, InertiaError> {
        checkpoint(cx)?;
        let nodes = self.ndofs() / 3;
        if self.ndofs() % 3 != 0
            || reference_nodes_m.len() != nodes
            || displacement_before_m.len() != nodes
            || displacement_after_m.len() != nodes
            || velocity_mps.len() != nodes
            || origin_m.iter().any(|x| !x.is_finite())
        {
            return Err(invalid("complete native kinetic drift geometry required"));
        }
        let mut linear = [Dd::ZERO; 3];
        let mut before_angular = [Dd::ZERO; 3];
        let mut after_angular = [Dd::ZERO; 3];
        let mut angular_change = [Dd::ZERO; 3];
        let mut kinetic = Dd::ZERO;
        for node in 0..nodes {
            if node % 256 == 0 {
                checkpoint(cx)?;
            }
            let mut before_arm = [Dd::ZERO; 3];
            let mut after_arm = [Dd::ZERO; 3];
            let mut travel = [Dd::ZERO; 3];
            let mut momentum = [Dd::ZERO; 3];
            for axis in 0..3 {
                before_arm[axis] = current_arm(
                    reference_nodes_m[node][axis],
                    displacement_before_m[node][axis],
                    origin_m[axis],
                )?;
                after_arm[axis] = current_arm(
                    reference_nodes_m[node][axis],
                    displacement_after_m[node][axis],
                    origin_m[axis],
                )?;
                travel[axis] = Dd::from_f64(displacement_after_m[node][axis])
                    - Dd::from_f64(displacement_before_m[node][axis]);
                momentum[axis] = self
                    .motion_momentum_row(3 * node + axis, velocity_mps, None, cx)?
                    .0;
                linear[axis] = linear[axis] + momentum[axis];
                kinetic = kinetic
                    + Dd::from_f64(0.5)
                        * Dd::from_f64(finite(velocity_mps[node][axis])?)
                        * momentum[axis];
            }
            let old = cross(before_arm, momentum);
            let next = cross(after_arm, momentum);
            let change = cross(travel, momentum);
            for axis in 0..3 {
                before_angular[axis] = before_angular[axis] + old[axis];
                after_angular[axis] = after_angular[axis] + next[axis];
                angular_change[axis] = angular_change[axis] + change[axis];
            }
        }
        let kinetic_j = finite(kinetic.to_f64())?;
        if kinetic_j < 0.0 {
            return Err(invalid("nonnegative native drift kinetic energy required"));
        }
        checkpoint(cx)?;
        Ok(ConstantVelocityDrift {
            before: momentum(origin_m, linear, before_angular)?,
            after: momentum(origin_m, linear, after_angular)?,
            change: momentum(origin_m, [Dd::ZERO; 3], angular_change)?,
            kinetic_j,
        })
    }

    /// One sparse traversal forms old and incremental momentum with compensated
    /// products. No mass-sized scratch, mass reconstruction, history copy,
    /// constitutive evaluation or new linear solve is required.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn backward_euler_momentum(
        &self,
        reference_nodes_m: &[[f64; 3]],
        displacement_before_m: &[[f64; 3]],
        displacement_after_m: &[[f64; 3]],
        velocity_before_mps: &[[f64; 3]],
        velocity_after_mps: &[[f64; 3]],
        external_force_n: &[f64],
        reactions_n: &[(usize, f64)],
        duration_s: f64,
        origin_m: [f64; 3],
        cx: &Cx<'_>,
    ) -> Result<SolidMotionMomentum, InertiaError> {
        checkpoint(cx)?;
        let nodes = self.ndofs() / 3;
        if !duration_s.is_finite()
            || duration_s <= 0.0
            || origin_m.iter().any(|x| !x.is_finite())
            || self.ndofs() % 3 != 0
            || reference_nodes_m.len() != nodes
            || displacement_before_m.len() != nodes
            || displacement_after_m.len() != nodes
            || velocity_before_mps.len() != nodes
            || velocity_after_mps.len() != nodes
            || external_force_n.len() != self.ndofs()
        {
            return Err(invalid(
                "complete native motion and positive momentum interval required",
            ));
        }
        let mut linear = [Dd::ZERO; 3];
        let mut angular = [Dd::ZERO; 3];
        let mut applied_linear = [Dd::ZERO; 3];
        let mut applied_angular = [Dd::ZERO; 3];
        let mut support_linear = [Dd::ZERO; 3];
        let mut support_angular = [Dd::ZERO; 3];
        let mut drift = [Dd::ZERO; 3];
        let mut linear_residual = [Dd::ZERO; 3];
        let mut angular_residual = [Dd::ZERO; 3];
        for node in 0..nodes {
            if node % 256 == 0 {
                checkpoint(cx)?;
            }
            let mut arm = [Dd::ZERO; 3];
            let mut travel = [Dd::ZERO; 3];
            let mut impulse = [Dd::ZERO; 3];
            let mut applied = [Dd::ZERO; 3];
            let mut previous = [Dd::ZERO; 3];
            for axis in 0..3 {
                // Reference and displacement are separately retained native
                // coordinates. Rounding their sum would erase real motion
                // before the declared origin or the old endpoint is removed.
                current_arm(
                    reference_nodes_m[node][axis],
                    displacement_before_m[node][axis],
                    origin_m[axis],
                )?;
                arm[axis] = current_arm(
                    reference_nodes_m[node][axis],
                    displacement_after_m[node][axis],
                    origin_m[axis],
                )?;
                travel[axis] = Dd::from_f64(displacement_after_m[node][axis])
                    - Dd::from_f64(displacement_before_m[node][axis]);
                (previous[axis], impulse[axis]) = self.motion_momentum_row(
                    3 * node + axis,
                    velocity_before_mps,
                    Some(velocity_after_mps),
                    cx,
                )?;
                applied[axis] = Dd::from_f64(finite(external_force_n[3 * node + axis])?)
                    * Dd::from_f64(duration_s);
                linear[axis] = linear[axis] + impulse[axis];
                applied_linear[axis] = applied_linear[axis] + applied[axis];
                linear_residual[axis] = linear_residual[axis] + (impulse[axis] - applied[axis]);
            }
            let endpoint = cross(arm, impulse);
            let correction = cross(travel, previous);
            let applied_moment = cross(arm, applied);
            let residual_moment = cross(arm, std::array::from_fn(|a| impulse[a] - applied[a]));
            for axis in 0..3 {
                angular[axis] = angular[axis] + endpoint[axis] + correction[axis];
                drift[axis] = drift[axis] + correction[axis];
                applied_angular[axis] = applied_angular[axis] + applied_moment[axis];
                angular_residual[axis] = angular_residual[axis] + residual_moment[axis];
            }
        }
        for (i, &(dof, force)) in reactions_n.iter().enumerate() {
            if i % 256 == 0 {
                checkpoint(cx)?;
            }
            if dof >= self.ndofs() {
                return Err(invalid(
                    "native support reaction lies outside original mass",
                ));
            }
            let impulse = Dd::from_f64(finite(force)?) * Dd::from_f64(duration_s);
            let mut arm = [Dd::ZERO; 3];
            let mut applied = [Dd::ZERO; 3];
            applied[dof % 3] = impulse;
            for axis in 0..3 {
                arm[axis] = current_arm(
                    reference_nodes_m[dof / 3][axis],
                    displacement_after_m[dof / 3][axis],
                    origin_m[axis],
                )?;
            }
            let moment = cross(arm, applied);
            for axis in 0..3 {
                support_linear[axis] = support_linear[axis] + applied[axis];
                support_angular[axis] = support_angular[axis] + moment[axis];
                linear_residual[axis] = linear_residual[axis] - applied[axis];
                angular_residual[axis] = angular_residual[axis] - moment[axis];
            }
        }
        checkpoint(cx)?;
        Ok(SolidMotionMomentum {
            change: momentum(origin_m, linear, angular)?,
            applied: momentum(origin_m, applied_linear, applied_angular)?,
            support: momentum(origin_m, support_linear, support_angular)?,
            angular_integration_defect_kg_m2_s: rounded(drift)?,
            linear_residual_kg_m_s: rounded(linear_residual)?,
            angular_residual_kg_m2_s: rounded(angular_residual)?,
        })
    }
    /// Contract the original CSR row before rounding any nodal momentum or
    /// velocity difference. Opposing large row contributions must not erase
    /// a small physical total merely because their individual inventories
    /// round to the same f64. The native solve and its matrix remain unchanged.
    fn motion_momentum_row(
        &self,
        row: usize,
        before: &[[f64; 3]],
        after: Option<&[[f64; 3]]>,
        cx: &Cx<'_>,
    ) -> Result<(Dd, Dd), InertiaError> {
        let (columns, values) = self.matrix().row(row);
        let mut previous = Dd::ZERO;
        let mut change = Dd::ZERO;
        for (entry, (&column, &mass)) in columns.iter().zip(values).enumerate() {
            if entry % 256 == 0 {
                checkpoint(cx)?;
            }
            let old = Dd::from_f64(finite(before[column / 3][column % 3])?);
            let mass = Dd::from_f64(mass);
            previous = previous + mass * old;
            if let Some(after) = after {
                let new = Dd::from_f64(finite(after[column / 3][column % 3])?);
                change = change + mass * (new - old);
            }
        }
        finite(previous.hi)?;
        finite(previous.lo)?;
        finite(change.hi)?;
        finite(change.lo)?;
        Ok((previous, change))
    }
}

fn current_arm(reference: f64, displacement: f64, origin: f64) -> Result<Dd, InertiaError> {
    let position = Dd::from_f64(finite(reference)?) + Dd::from_f64(finite(displacement)?);
    finite(position.hi)?;
    finite(position.lo)?;
    let arm = position - Dd::from_f64(origin);
    finite(arm.hi)?;
    finite(arm.lo)?;
    Ok(arm)
}

fn cross(a: [Dd; 3], b: [Dd; 3]) -> [Dd; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn rounded(a: [Dd; 3]) -> Result<[f64; 3], InertiaError> {
    Ok([
        finite(a[0].to_f64())?,
        finite(a[1].to_f64())?,
        finite(a[2].to_f64())?,
    ])
}

fn momentum(
    origin_m: [f64; 3],
    linear: [Dd; 3],
    angular: [Dd; 3],
) -> Result<SolidMomentum, InertiaError> {
    Ok(SolidMomentum {
        origin_m,
        linear_kg_m_s: rounded(linear)?,
        angular_kg_m2_s: rounded(angular)?,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_alloc::{ArenaConfig, ArenaPool};
    use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
    use fs_sparse::Coo;

    fn with_cx(run: impl FnOnce(&Cx<'_>)) {
        let gate = CancelGate::new_clock_free();
        ArenaPool::new(ArenaConfig::default()).scope(|arena| {
            run(&Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 31,
                    kernel_id: 7,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ));
        });
    }

    // Symmetric positive consistent mass (I + 11^T) tensor I_3. Every
    // scalar row has mass n+1; cross-node entries remain present.
    fn mass(nodes: usize, cx: &Cx<'_>) -> ReferenceMass {
        let mut coo = Coo::new(3 * nodes, 3 * nodes);
        for a in 0..nodes {
            for b in 0..nodes {
                for axis in 0..3 {
                    coo.push(3 * a + axis, 3 * b + axis, if a == b { 2.0 } else { 1.0 });
                }
            }
        }
        ReferenceMass::from_matrix(coo.assemble(), (nodes * (nodes + 1)) as f64, 3 * nodes, cx)
            .unwrap()
    }

    /// G0/G3: a common world translation cannot erase retained deformation.
    #[test]
    fn drift_preserves_deformed_moment_under_large_common_translation() {
        with_cx(|cx| {
            let mass = mass(2, cx);
            let before = [[0.0, 0.5, 0.0]; 2];
            let after = [[0.25, 1.0, 0.0]; 2];
            let velocity = [[1.0, 2.0, 0.0]; 2];
            for shift in [[0.0; 3], [2.0_f64.powi(54), -2.0_f64.powi(54), 0.0]] {
                let reference = [shift, [shift[0] + 8.0, shift[1], shift[2]]];
                let drift = mass
                    .constant_velocity_drift(&reference, &before, &after, &velocity, shift, cx)
                    .unwrap();
                // p_a = [3,6,0]; sum((X_a+u_a) cross p_a)_z = 45.
                assert_eq!(drift.before.linear_kg_m_s, [6.0, 12.0, 0.0]);
                assert_eq!(drift.before.angular_kg_m2_s, [0.0, 0.0, 45.0]);
                assert_eq!(drift.after.angular_kg_m2_s, drift.before.angular_kg_m2_s);
                assert_eq!(drift.change.linear_kg_m_s, [0.0; 3]);
                assert_eq!(drift.change.angular_kg_m2_s, [0.0; 3]);
                assert_eq!(drift.kinetic_j, 15.0);
            }
        });
    }

    /// G0/G3: applied and support moments use the same unrounded endpoint arm.
    #[test]
    fn backward_euler_retains_load_support_and_travel_moments_after_translation() {
        with_cx(|cx| {
            let mass = mass(2, cx);
            let before = [[0.0, 0.5, 0.0]; 2];
            let after = [[0.25, 1.0, 0.0]; 2];
            for shift in [[0.0; 3], [2.0_f64.powi(54), -2.0_f64.powi(54), 0.0]] {
                let reference = [shift, [shift[0] + 8.0, shift[1], shift[2]]];
                let result = mass
                    .backward_euler_momentum(
                        &reference,
                        &before,
                        &after,
                        &[[1.0, 0.0, 0.0]; 2],
                        &[[1.0, 2.0, 0.0]; 2],
                        &[0.0, 12.0, 0.0, 0.0, 24.0, 0.0],
                        &[(1, 12.0)],
                        0.25,
                        shift,
                        cx,
                    )
                    .unwrap();
                assert_eq!(result.change.linear_kg_m_s, [0.0, 12.0, 0.0]);
                assert_eq!(result.change.angular_kg_m2_s, [0.0, 0.0, 48.0]);
                assert_eq!(result.applied.linear_kg_m_s, [0.0, 9.0, 0.0]);
                assert_eq!(result.applied.angular_kg_m2_s, [0.0, 0.0, 50.25]);
                assert_eq!(result.support.linear_kg_m_s, [0.0, 3.0, 0.0]);
                assert_eq!(result.support.angular_kg_m2_s, [0.0, 0.0, 0.75]);
                assert_eq!(result.angular_integration_defect_kg_m2_s, [0.0, 0.0, -3.0]);
                assert_eq!(result.linear_residual_kg_m_s, [0.0; 3]);
                assert_eq!(result.angular_residual_kg_m2_s, [0.0; 3]);
            }
        });
    }

    /// G0: signed sparse-row cancellation must retain the physical inventory.
    #[test]
    fn counterflow_retains_small_total_momentum_before_nodal_rounding() {
        with_cx(|cx| {
            let mass = mass(3, cx);
            let common = 2.0_f64.powi(54);
            let velocity = [[common, 0.0, 0.0], [1.0, 0.0, 0.0], [-common, 0.0, 0.0]];
            let drift = mass
                .constant_velocity_drift(
                    &[[0.0; 3]; 3],
                    &[[0.0; 3]; 3],
                    &[[0.0; 3]; 3],
                    &velocity,
                    [0.0; 3],
                    cx,
                )
                .unwrap();
            // Each node's scalar row sum is 4, so P = 4*sum(v) = [4,0,0].
            assert_eq!(drift.before.linear_kg_m_s, [4.0, 0.0, 0.0]);
            assert_eq!(drift.after, drift.before);
            assert_eq!(drift.change.linear_kg_m_s, [0.0; 3]);
        });
    }

    /// G0: an endpoint difference is taken before any f64 increment rounding.
    #[test]
    fn endpoint_velocity_difference_cannot_hide_unbalanced_momentum() {
        with_cx(|cx| {
            let mass = mass(2, cx);
            let common = 2.0_f64.powi(54);
            let after = [[1.0, 0.0, 0.0], [0.0; 3]];
            let result = mass
                .backward_euler_momentum(
                    &[[0.0; 3]; 2],
                    &[[0.0; 3]; 2],
                    &after,
                    &[[common, 0.0, 0.0], [-common, 0.0, 0.0]],
                    &after,
                    &[0.0; 6],
                    &[],
                    1.0,
                    [0.0; 3],
                    cx,
                )
                .unwrap();
            // Old total momentum is zero, new total momentum is 3. The
            // absent forcing must remain an actual nonzero balance residual.
            assert_eq!(result.change.linear_kg_m_s, [3.0, 0.0, 0.0]);
            assert_eq!(result.linear_residual_kg_m_s, [3.0, 0.0, 0.0]);
        });
    }

    /// G0: preserve real stored-displacement defects rather than imposing an
    /// ideal dt*v trajectory on already rounded accepted coordinates.
    #[test]
    fn actual_displacement_rounding_remains_a_numerical_angular_defect() {
        with_cx(|cx| {
            let mass = mass(2, cx);
            let common = 2.0_f64.powi(54);
            let before = [[common, 0.0, 0.0]; 2];
            // fl(common + 0.25) == common, while the y increment survives.
            let after = [[common + 0.25, 0.25, 0.0]; 2];
            let drift = mass
                .constant_velocity_drift(
                    &[[0.0; 3]; 2],
                    &before,
                    &after,
                    &[[1.0, 1.0, 0.0]; 2],
                    [0.0; 3],
                    cx,
                )
                .unwrap();
            assert_eq!(drift.change.angular_kg_m2_s, [0.0, 0.0, -1.5]);
        });
    }
}
