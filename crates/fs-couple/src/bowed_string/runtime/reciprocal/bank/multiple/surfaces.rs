//! A material string rubbing a dynamically resolved receiving body.
//!
//! This is the missing opposite side of the existing friction port, not a new
//! exciter law. A supplied compliant bow, finger or other body may receive the
//! reaction through its original geometry/material basis. No hair construction,
//! rosin behavior, material identity or measurement authority is inferred.

use super::*;
use crate::render::schedule::force::coupled::contact::friction::multiple::surfaces::FrictionSurfaceAttachment;

/// Tangential surface motion of an already retained physical body.
#[derive(Clone, Debug, PartialEq)]
pub struct BowBodySurface {
    /// Existing bow/contact index in original drive order.
    pub bow: usize,
    /// Original body index (before strings are prepended) and its physical
    /// displacement map along the string contact's positive tangent.
    /// The body producer must derive this row from the same basis as its mass,
    /// stiffness and material loss; a material swap must reproject this row.
    pub body: ModalAttachment,
}

impl MultiBowStringBankSource {
    /// Original body surfaces, retained alongside the string/material sources.
    #[must_use]
    pub fn body_surfaces(&self) -> &[BowBodySurface] {
        &self.body_surfaces
    }
}

impl MultiBowStringBank {
    /// Return each contact's tangential reaction to its actual moving receiver.
    ///
    /// The receiving body is already part of the original bank. Its inertia,
    /// material loss and compliance enter the SAME implicit solve as every
    /// string, bridge, cavity, gas-memory coordinate and other contact. Multiple
    /// contacts may share a receiver; their cross-reactions are recomputed.
    /// Existing displacement, velocity, pressure maps and source bindings stay
    /// intact. Cavity-neck, fluid-memory and longitudinal helper coordinates cannot be
    /// relabeled as physical bodies.
    ///
    /// Bind before normal stops and playback. Every contact may bind once;
    /// its original normal-force and interface laws remain in effect. Gesture
    /// speed is an additional prescribed relative actuator speed. Set it to
    /// zero when all surface motion should arise from the solved mechanics and
    /// physical external forces. No acoustic normalization or pitch change is
    /// introduced by this adapter.
    pub fn with_bow_body_surfaces(
        mut self,
        surfaces: Vec<BowBodySurface>,
        gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        poll(gate)?;
        if surfaces.is_empty() {
            return Ok(self);
        }
        if surfaces.len() > self.source.bows.len() {
            return Err(invalid("bow surfaces exceed the retained contact count"));
        }
        let internal = self
            .source
            .bank
            .fluid_loading
            .as_ref()
            .map_or(0, |f| f.strings().len())
            + self
                .source
                .bank
                .longitudinal_waves
                .as_ref()
                .map_or(0, |l| l.strings().len())
            + usize::from(self.source.bank.material_cavity.is_some());
        let internal = internal + self.source.fiber_ribbon.as_ref()
            .map_or(0, |ribbon| ribbon.model().resolution().cells);
        let body_count = self
            .system
            .components()
            .len()
            .checked_sub(self.cards.len())
            .and_then(|count| count.checked_sub(internal))
            .ok_or_else(|| invalid("bow surface lost the original physical-body layout"))?;
        let mut attachments = Vec::with_capacity(surfaces.len());
        for surface in &surfaces {
            poll(gate)?;
            if surface.bow >= self.source.bows.len() || surface.body.component >= body_count {
                return Err(invalid(
                    "bow surface must name an existing contact and original physical body",
                ));
            }
            let mut body = surface.body.clone();
            body.component += self.cards.len();
            attachments.push(FrictionSurfaceAttachment {
                port: surface.bow,
                surface: body,
            });
        }
        self.system = self
            .system
            .with_surface_attachments(attachments, gate)
            .map_err(ReciprocalStringError::Friction)?;
        self.source.body_surfaces.extend(surfaces);
        Ok(self)
    }
}
