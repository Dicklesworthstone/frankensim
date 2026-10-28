//! Prescribed supports as a projection of the existing rod equilibrium.
//!
//! Root node zero retains its original clamped position/frame. Additional
//! translational constraints are in the world frame. A rotation constraint
//! fixes the complete SO(3) frame, not selected quaternion components. No
//! penalty springs or stiffness-dependent artificial constraint coefficients.

use fs_exec::CancelGate;
use fs_ga::So3;
use super::{SpatialMaterialRod, SpatialRodError, checkpoint};
use crate::rod::TipLoad;

/// Additional ideal, fixed nodal support; geometry and loads remain independent
/// of material-card authority. None leaves that coordinate unconstrained.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct RodSupport {
    /// Node in 1..=segment_count; root zero is already clamped.
    pub node: usize,
    /// Prescribed world position components [m].
    pub position_m: [Option<f64>; 3],
    /// Complete prescribed body-to-world frame, or free rotation.
    pub frame: Option<So3>,
}

/// Reactions in the constrained coordinates of one support (including root).
#[derive(Clone, Debug, PartialEq)]
pub struct RodSupportReaction {
    /// Original node index.
    pub node: usize,
    /// Force exerted by the support on the rod, world components [N].
    pub force_n: [Option<f64>; 3],
    /// Moment exerted by a frame clamp, in that node's BODY frame [N m].
    pub moment_body_nm: Option<[f64; 3]>,
}

/// Numerical reaction estimates and the unbalanced free-coordinate residual.
/// Root reaction is available separately through `root_support_reaction`.
/// This report itself is not a whole-body force/moment balance.
#[derive(Clone, Debug, PartialEq)]
pub struct RodSupportReactions {
    /// Additional supports in canonical node order.
    pub nodes: Vec<RodSupportReaction>,
    /// Original numerical residual norm (mixed force/moment units) on free DOFs.
    /// Reactions represent equilibrium only when this residual is acceptable.
    pub free_residual_norm: f64,
}

impl SpatialMaterialRod {
    /// Return an independently supported static specimen, preserving materials.
    /// The support list replaces the old list; newly freed coordinates retain
    /// their old pose as a numerical initial guess. It is not a time-dependent
    /// imposed motion, transported history, or work-conserving support change.
    /// Constraints stay exact during every accepted trial and tangent probe.
    ///
    /// # Errors
    /// Refuses invalid/duplicate/root node indices, empty supports, nonfinite
    /// positions, material strain-domain failures, and cancellation. The input
    /// specimen remains unchanged, including on a late thermal strain refusal.
    pub fn with_supports(
        &self, supports: &[RodSupport], gate: &CancelGate,
    ) -> Result<Self, SpatialRodError> {
        checkpoint(gate)?;
        if supports.len() > self.sections.len() {
            return Err(SpatialRodError::Invalid("too many additional rod supports"));
        }
        let mut ordered = supports.to_vec();
        ordered.sort_by_key(|s| s.node);
        let mut candidate = self.clone();
        let mut fixed = vec![false; self.rod.ndof()];
        for (i, support) in ordered.iter().enumerate() {
            checkpoint(gate)?;
            if support.node == 0 || support.node > self.sections.len()
                || (i > 0 && ordered[i - 1].node == support.node)
                || (support.frame.is_none() && support.position_m.iter().all(Option::is_none))
                || support.position_m.iter().flatten().any(|x| !x.is_finite())
            {
                return Err(SpatialRodError::Invalid("rod support needs a unique nonroot node and finite prescribed coordinates"));
            }
            let k = 6 * (support.node - 1);
            for (c, position) in support.position_m.iter().enumerate() {
                if let Some(value) = position {
                    candidate.rod.positions[support.node][c] = *value;
                    fixed[k + c] = true;
                }
            }
            if let Some(frame) = support.frame {
                candidate.rod.frames[support.node] = frame;
                fixed[k + 3..k + 6].fill(true);
            }
        }
        candidate.free_dofs = (0..fixed.len()).filter(|&k| !fixed[k]).collect();
        candidate.supports = ordered;
        candidate.equilibrium_energy_for(&candidate.rod)?;
        checkpoint(gate)?;
        Ok(candidate)
    }

    /// Prescribed additional supports, sorted by node; root zero is implicit.
    #[must_use]
    pub fn supports(&self) -> &[RodSupport] { &self.supports }

    /// Recover `dE/dq - applied_load` using the SAME element-local body/right
    /// derivatives and analytical link/gravity loads as the Newton solve.
    /// Pass the load used to solve this pose. The free residual prevents a
    /// reaction-only report from hiding an unequilibrated configuration.
    /// Segment probes retain their source strain bounds; reactions are numerical
    /// estimates, not forced whole-body balance closures.
    ///
    /// # Errors
    /// Refuses nonfinite loads/residuals, source/kinematic failures and cancellation.
    pub fn support_reactions(
        &self, load: &TipLoad, gate: &CancelGate,
    ) -> Result<RodSupportReactions, SpatialRodError> {
        checkpoint(gate)?;
        if !load.force.iter().chain(&load.moment).all(|x| x.is_finite()) {
            return Err(SpatialRodError::Invalid("reaction loads must be finite"));
        }
        self.equilibrium_energy_for(&self.rod)?;
        let energy = super::derivatives::SpatialEnergy::new(self);
        let residual = self.rod.residual(load, 1.0, &energy,
            &|| checkpoint(gate), self.probe_translation_scale_m())?;
        let free_residual_norm = self.free_dofs.iter().map(|&k| residual[k] * residual[k]).sum::<f64>().sqrt();
        if !residual.iter().all(|v| v.is_finite()) || !free_residual_norm.is_finite() {
            return Err(SpatialRodError::Invalid("support reaction or free residual is unrepresentable"));
        }
        let nodes = self.supports.iter().map(|support| {
            let k = 6 * (support.node - 1);
            RodSupportReaction {
                node: support.node,
                force_n: [0, 1, 2].map(|c| support.position_m[c].map(|_| residual[k + c])),
                moment_body_nm: support.frame.map(|_| [residual[k + 3], residual[k + 4], residual[k + 5]]),
            }
        }).collect();
        checkpoint(gate)?;
        Ok(RodSupportReactions { nodes, free_residual_norm })
    }
}
