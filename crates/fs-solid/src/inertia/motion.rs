//! Original-mass momentum over a backward-Euler material interval.
use fs_math::dd::Dd;

use super::{Cx, InertiaError, ReferenceMass, SolidMomentum, checkpoint, finite, invalid};

/// Actual physical increments and endpoint loads about one declared origin.
/// Backward Euler's angular drift is a numerical defect, never a torque source.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SolidMotionMomentum {
    pub change: SolidMomentum,
    pub applied: SolidMomentum,
    pub support: SolidMomentum,
    /// Sum(delta_x cross M v_before), including coordinate rounding. For the
    /// exact BE displacement relation and symmetric mass this equals
    /// -dt sum(v_after cross M delta_v). It is not physical angular impulse.
    pub angular_integration_defect_kg_m2_s: [f64; 3],
    pub linear_residual_kg_m_s: [f64; 3],
    /// Actual angular change minus applied/support impulse and the separately
    /// reported numerical defect. Consumers must also bound that defect itself.
    pub angular_residual_kg_m2_s: [f64; 3],
}

impl ReferenceMass {
    /// Three mass-sized vectors and two sparse actions; no mass reconstruction,
    /// history copy, constitutive evaluation or new linear solve is required.
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
        let mut delta_v = Vec::with_capacity(self.ndofs());
        for (node, (before, after)) in velocity_before_mps
            .iter()
            .zip(velocity_after_mps)
            .enumerate()
        {
            if node % 256 == 0 {
                checkpoint(cx)?;
            }
            for axis in 0..3 {
                delta_v.push(finite(after[axis] - before[axis])?);
            }
        }
        let delta_p = self.action(&delta_v, cx)?;
        let old_p = self.action(velocity_before_mps.as_flattened(), cx)?;
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
                let before =
                    finite(reference_nodes_m[node][axis] + displacement_before_m[node][axis])?;
                let after =
                    finite(reference_nodes_m[node][axis] + displacement_after_m[node][axis])?;
                arm[axis] = Dd::from_f64(after) - Dd::from_f64(origin_m[axis]);
                travel[axis] = Dd::from_f64(after) - Dd::from_f64(before);
                impulse[axis] = Dd::from_f64(delta_p[3 * node + axis]);
                previous[axis] = Dd::from_f64(old_p[3 * node + axis]);
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
                arm[axis] = Dd::from_f64(finite(
                    reference_nodes_m[dof / 3][axis] + displacement_after_m[dof / 3][axis],
                )?) - Dd::from_f64(origin_m[axis]);
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
