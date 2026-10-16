//! Reciprocal felt support on a supplied mass-normalized mechanical body.
//!
//! A backing may translate, rock, or bend. Its point-motion row comes from the
//! SAME rigid/elastic basis as its inertia and stiffness; it is never normalized
//! to unit length. Pad history and the support are solved by ImpactSystem, not
//! independently. Preparation is cold and does not add per-step allocations.
use super::super::FeltPad;
use super::super::super::{ImpactBody, ImpactError, MAX_IMPACT_MODES, invalid};

/// One physical patch, referred to the original structure and the new support.
#[derive(Debug, Clone)]
pub struct ModalMountSite {
    pub pad: FeltPad,
    /// Signed d(compression)/d(support q) in the supplied mass-normalized basis.
    /// Rotational and bending rows include their actual lever/shape coefficients.
    pub support_weights: Vec<f64>,
}

/// Append this body after the original structure and REPLACE its fixed pads.
/// Every unrelated force/contact/volume/observer row must be zero-extended.
#[derive(Debug)]
pub struct ModalMountedFelt {
    pub body: ImpactBody,
    pub pads: Vec<FeltPad>,
}
impl ModalMountedFelt {
    /// Compose an existing mechanical body without changing its potential,
    /// initial conditions, modal signs, damping, or pad/material history.
    /// ImpactSystem still owns timestep, force, energy and history admission.
    pub fn new(structural_modes: usize, body: ImpactBody, sites: &[ModalMountSite])
        -> Result<Self, ImpactError>
    {
        let support_modes = body.potential.count();
        let total = structural_modes.checked_add(support_modes)
            .ok_or_else(|| invalid("felt support extent overflow"))?;
        if structural_modes == 0 || support_modes == 0 || total > MAX_IMPACT_MODES
            || sites.is_empty() || sites.len() > 16
            || body.initial.len() != support_modes || body.damping_per_s.len() != support_modes
            || body.potential.omegas().iter().any(|w| !w.is_finite() || *w < 0.0 || !(w*w).is_finite())
            || body.damping_per_s.iter().any(|r| !r.is_finite() || *r < 0.0)
            || body.initial.iter().any(|s| !s.displacement_m_sqrt_kg.is_finite()
                || !s.velocity_m_sqrt_kg_per_s.is_finite()) {
            return Err(invalid("felt support needs a bounded passive body and complete initial/modal rows"));
        }
        let mut pads = Vec::with_capacity(sites.len());
        let mut coupled = false;
        for site in sites {
            site.pad.validate(structural_modes)?;
            if site.support_weights.len() != support_modes
                || site.support_weights.iter().any(|b| !b.is_finite())
                || site.pad.weights.iter().all(|b| *b == 0.0) {
                return Err(invalid("felt support needs actual structural participation and complete finite support rows"));
            }
            coupled |= site.support_weights.iter().any(|b| *b != 0.0);
            let mut pad = site.pad.clone();
            pad.weights.extend_from_slice(&site.support_weights);
            pad.validate(total)?;
            pads.push(pad);
        }
        if !coupled { return Err(invalid("all felt patches are decoupled from their support")); }
        Ok(Self { body, pads })
    }
}

#[cfg(test)]
mod tests;
