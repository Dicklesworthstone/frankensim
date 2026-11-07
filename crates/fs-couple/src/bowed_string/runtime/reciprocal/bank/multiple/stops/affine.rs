//! Exact coordinate substitution for a fixed condensed point displacement.
//! U(x+offset-gap) = U(x-(gap-offset)); the original source law remains in the
//! material bank. Every native secant, source check and preload sees this same
//! clearance, while force maps, storage ownership and constant-offset work stay
//! unchanged. A negative effective clearance represents initial interference.

use super::*;
use fs_dcontact::OpeningContactStep;

pub(super) fn lower(law: &Obstacle, offset_m: f64) -> Result<Obstacle, RenderError> {
    if !offset_m.is_finite() {
        return Err(error("residual contact displacement offset must be finite"));
    }
    if offset_m == 0.0 {
        return Ok(law.clone());
    }
    // Reuse the actual unit-point law admission before indexing raw parts.
    let contact_error = |source| RenderError::Coupled(ModalCouplingError::ContactLaw(source));
    OpeningContactStep::new(law, 0.0).map_err(contact_error)?;
    Obstacle::new(
        law.collocation().to_vec(),
        1,
        1,
        vec![law.gaps()[0] - offset_m],
        law.weights().to_vec(),
        law.stiffness(),
        law.alpha(),
        law.provenance().to_owned(),
    )
    .and_then(|lowered| lowered.with_internal_loss(law.internal_loss()))
    .map_err(contact_error)
}

#[cfg(test)]
mod tests;
