//! Interface system cards (bead 5hmy, PR-3 of 5).
//!
//! Friction, wetting, contact conductance, wear, and adhesion are
//! SYSTEM + HISTORY properties: they belong to an ORDERED pair of
//! named surfaces, the medium between them, any third body, the
//! environment, and the named history state — never to an unordered
//! pair of bulk materials. "Steel on PTFE, dry, run-in, in air" and
//! "PTFE on steel, oil-lubricated, virgin, in argon" are DIFFERENT
//! systems with different cards. Wetting is a solid–liquid–gas system:
//! the liquid is the medium, the gas is the environment.

use fs_blake3::{ContentHash, hash_domain};

use crate::cards::{ConstitutiveModelCard, MaterialStateId};
use crate::{ClaimId, ClaimSet, MatDbError, PropertyClaim};

/// Hash domain for interface-system-card canonical identity.
const INTERFACE_HASH_DOMAIN: &str = "org.frankensim.fs-matdb.interface-system-card.v1";

/// One side of an interface: a named material state plus its surface
/// texture frame. The texture id is OPAQUE at L1 — roughness spectra,
/// lay, and coating stacks are property claims against the system, and
/// the frame id only names which measured texture they refer to.
#[derive(Debug, Clone, PartialEq)]
pub struct SurfaceSpec {
    /// The named bulk material state this surface is cut from.
    pub material: MaterialStateId,
    /// Opaque texture-frame identifier (empty refuses: an interface
    /// claim against an unnamed surface state is unreproducible).
    pub texture_frame: String,
}

impl SurfaceSpec {
    fn validate(&self) -> Result<(), MatDbError> {
        if self.texture_frame.trim().is_empty() {
            return Err(MatDbError::MissingTextureFrame {
                material: self.material.clone(),
            });
        }
        Ok(())
    }
}

/// The system context an interface lives in: medium, optional third
/// body, environment, and the NAMED history state. Grouped because the
/// four travel together — a context with a blank member is an
/// incomplete system identity and refuses at assembly.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SystemContext {
    /// The intervening medium ("dry", "air gap", a liquid for wetting).
    pub medium: String,
    /// Lubricant or third body, when present.
    pub third_body: Option<String>,
    /// The surrounding environment ("air", "argon", "vacuum").
    pub environment: String,
    /// The NAMED history state ("virgin", "run-in-1000-cycles",
    /// "aged-500h-85C"). History is identity-bearing, not a footnote.
    pub history: String,
}

impl SystemContext {
    fn validate(&self) -> Result<(), MatDbError> {
        for (field, value) in [
            ("medium", self.medium.as_str()),
            ("environment", self.environment.as_str()),
            ("history", self.history.as_str()),
        ] {
            if value.trim().is_empty() {
                return Err(MatDbError::MissingSystemField { field });
            }
        }
        if let Some(third_body) = &self.third_body
            && third_body.trim().is_empty()
        {
            return Err(MatDbError::MissingSystemField {
                field: "third_body",
            });
        }
        Ok(())
    }
}

/// The immutable card for one ORDERED interface system. Constructed
/// only by [`InterfaceSystemCard::assemble`]; no mutable access exists
/// afterward.
#[derive(Debug, Clone, PartialEq)]
pub struct InterfaceSystemCard {
    surface_a: SurfaceSpec,
    surface_b: SurfaceSpec,
    context: SystemContext,
    claims: ClaimSet,
    models: Vec<ConstitutiveModelCard>,
}

impl InterfaceSystemCard {
    /// Assemble an interface system card. Order is semantic: `(a, b)`
    /// and `(b, a)` are different systems (directional friction, lay
    /// orientation), and callers claiming symmetry must say so with two
    /// cards or a symmetric claim on each.
    ///
    /// # Errors
    /// [`MatDbError::MissingTextureFrame`] for an unnamed surface
    /// state; [`MatDbError::MissingSystemField`] for a blank
    /// medium/environment/history/third-body; model-card refusals as in
    /// [`crate::MaterialCard::assemble`].
    pub fn assemble(
        surface_a: SurfaceSpec,
        surface_b: SurfaceSpec,
        context: SystemContext,
        claims: ClaimSet,
        models: Vec<ConstitutiveModelCard>,
    ) -> Result<InterfaceSystemCard, MatDbError> {
        surface_a.validate()?;
        surface_b.validate()?;
        context.validate()?;
        for model in &models {
            model.validate()?;
        }
        Ok(InterfaceSystemCard {
            surface_a,
            surface_b,
            context,
            claims,
            models,
        })
    }

    /// The first (ordered) surface.
    #[must_use]
    pub fn surface_a(&self) -> &SurfaceSpec {
        &self.surface_a
    }

    /// The second (ordered) surface.
    #[must_use]
    pub fn surface_b(&self) -> &SurfaceSpec {
        &self.surface_b
    }

    /// The system context (medium, third body, environment, history).
    #[must_use]
    pub fn context(&self) -> &SystemContext {
        &self.context
    }

    /// The intervening medium.
    #[must_use]
    pub fn medium(&self) -> &str {
        &self.context.medium
    }

    /// The lubricant/third body, when present.
    #[must_use]
    pub fn third_body(&self) -> Option<&str> {
        self.context.third_body.as_deref()
    }

    /// The surrounding environment.
    #[must_use]
    pub fn environment(&self) -> &str {
        &self.context.environment
    }

    /// The named history state.
    #[must_use]
    pub fn history(&self) -> &str {
        &self.context.history
    }

    /// The system's property claims (friction, conductance, wetting
    /// hysteresis, wear — claimed against the SYSTEM, not the bulks).
    #[must_use]
    pub fn claims(&self) -> &ClaimSet {
        &self.claims
    }

    /// EVERY claim for a property name, in insertion order.
    #[must_use]
    pub fn claims_for(&self, name: &str) -> Vec<(ClaimId, &PropertyClaim)> {
        self.claims.claims_for(name)
    }

    /// The system's constitutive model cards (e.g. a friction law).
    #[must_use]
    pub fn models(&self) -> &[ConstitutiveModelCard] {
        &self.models
    }

    /// Canonical content identity: both ordered surfaces (material
    /// state id + texture frame), medium, third body, environment,
    /// history, every claim/observation content id, every model hash.
    #[must_use]
    pub fn content_hash(&self) -> ContentHash {
        let mut payload = Vec::new();
        let mut push = |part: &[u8]| {
            payload.extend_from_slice(&u64::try_from(part.len()).unwrap_or(u64::MAX).to_le_bytes());
            payload.extend_from_slice(part);
        };
        for surface in [&self.surface_a, &self.surface_b] {
            push(surface.material.chemistry.as_bytes());
            push(surface.material.phase.as_bytes());
            push(surface.material.process.as_bytes());
            push(&surface.material.revision.to_le_bytes());
            push(surface.texture_frame.as_bytes());
        }
        push(self.context.medium.as_bytes());
        match &self.context.third_body {
            Some(third_body) => push(third_body.as_bytes()),
            None => push(b"none"),
        }
        push(self.context.environment.as_bytes());
        push(self.context.history.as_bytes());
        for (claim_id, _) in self.claims.claims_ordered() {
            push(&claim_id.0.0);
        }
        for observation_id in self.claims.observation_ids() {
            push(&observation_id.0.0);
        }
        for model in &self.models {
            push(&model.content_hash().0);
        }
        hash_domain(INTERFACE_HASH_DOMAIN, &payload)
    }
}
