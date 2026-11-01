//! Bind the existing pure response laws to each native affine-load grid.
use super::{AffineResponseOptions3, admit_bound, admit_policy};
use crate::sdf3::response::limits::refinement::{
    ReferenceBoundedResponse3, ReferenceResponseLimit3,
};
use crate::sdf3::response::limits::{BoundedResponse3, ResponseLimit3};
use crate::sdf3::{AdaptiveSdf3Elasticity, Sdf3Elasticity, failure};
use crate::sdf3_goal::{GoalRefinementError3, controlled_goal3};
use crate::{EvaluationStop, SolveControl};

/// One unchanged ordered family of physical inequalities over the entire
/// affine load ball. Reuse the response-limit owner's observation and bound
/// types. Displacement laws are independently integrated on each native grid;
/// reaction modes remain live material-dependent observations during fitting.
#[derive(Clone, Copy)]
pub struct ReferenceAffineResponseLimits3<'a> {
    pub limits: &'a [ReferenceResponseLimit3<'a>],
    pub options: AffineResponseOptions3,
}

/// Admit the complete family before invoking any physical reference callback.
/// The same policy/bound helpers serve nodal evaluation and reference binding.
pub(in crate::sdf3::robust::envelope::affine_load_space) fn admit(
    operator: &impl Sdf3Elasticity,
    family: ReferenceAffineResponseLimits3<'_>,
) -> Result<(), EvaluationStop> {
    admit_policy(family.options, family.limits.len())?;
    for limit in family.limits {
        admit_bound(limit.bound, limit.scale)?;
        if matches!(limit.response, ReferenceBoundedResponse3::Reaction(_))
            && operator.embedded_dirichlet_penalty().is_none()
        {
            return Err(failure("sdf3-affine-response-embedded-support"));
        }
    }
    Ok(())
}

/// Own only the grid-dependent displacement observation vectors. The law
/// family and reaction callbacks stay borrowed; no density-dependent reaction
/// is converted into a frozen nodal row.
pub(in crate::sdf3::robust::envelope::affine_load_space) struct AssembledAffineResponseLimits3<'a> {
    family: ReferenceAffineResponseLimits3<'a>,
    observations: Vec<Option<Vec<f64>>>,
}
impl AssembledAffineResponseLimits3<'_> {
    pub(in crate::sdf3::robust::envelope::affine_load_space) fn limits(
        &self,
    ) -> Vec<ResponseLimit3<'_>> {
        self.family
            .limits
            .iter()
            .zip(&self.observations)
            .map(|(limit, observation)| ResponseLimit3 {
                response: match limit.response {
                    ReferenceBoundedResponse3::Displacement(_) => BoundedResponse3::Displacement(
                        observation
                            .as_ref()
                            .expect("assembled displacement observation"),
                    ),
                    ReferenceBoundedResponse3::Reaction(mode) => BoundedResponse3::Reaction(mode),
                },
                sense: limit.sense,
                bound: limit.bound,
                scale: limit.scale,
            })
            .collect()
    }
}

/// Integrate pure displacement laws on this actual grid, with the original
/// SolveControl budget/cancellation stop preserved by the shared native bridge.
pub(in crate::sdf3::robust::envelope::affine_load_space) fn assemble<'a>(
    operator: &impl AdaptiveSdf3Elasticity,
    family: ReferenceAffineResponseLimits3<'a>,
    control: &mut SolveControl<'_>,
) -> Result<AssembledAffineResponseLimits3<'a>, GoalRefinementError3> {
    admit(operator, family)?;
    let mut observations = Vec::with_capacity(family.limits.len());
    for limit in family.limits {
        control.checkpoint("sdf3-affine-reference-response")?;
        observations.push(match limit.response {
            ReferenceBoundedResponse3::Displacement(law) => Some(controlled_goal3(
                control,
                "sdf3-affine-reference-response",
                |poll| {
                    operator
                        .adaptive()
                        .reference_load(law, poll)
                        .map_err(Into::into)
                },
            )?),
            ReferenceBoundedResponse3::Reaction(_) => None,
        });
    }
    Ok(AssembledAffineResponseLimits3 {
        family,
        observations,
    })
}
