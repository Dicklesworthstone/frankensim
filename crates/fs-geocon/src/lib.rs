//! fs-geocon — geometric constraint primitives (plan §7.6). Layer: L4.
//!
//! The manufacturability and design-intent constraints ASCENT enforces,
//! each FIRST-CLASS: a value, a derivative story, a certificate story,
//! and its fs-constraint kind mapping ([`GeoPrimitive::descriptor`]) —
//! geometry supplies values and derivatives; ASCENT owns the constraint
//! SEMANTICS.
//!
//! - [`min_thickness_soft`] — the anti-paperclip constraint: a smooth
//!   p-norm aggregation of fs-query thickness samples (C¹ where walls
//!   are smooth), with violations LOCALIZED to their samples;
//! - [`draft_violations`] — casting/molding: normals versus the pull
//!   direction, smooth hinge² penalty plus EXACT violating regions
//!   (undercuts flagged separately and harder);
//! - [`QuotientChart`] — symmetry ENFORCED BY CONSTRUCTION: the
//!   evaluation point folds into the group's fundamental domain, so
//!   the shape is invariant for ARBITRARY design levers — violation is
//!   structurally impossible, not penalized (the plan's preferred
//!   mechanism);
//! - [`envelope_violation`] — containment/keep-out via SDF composition
//!   over boundary samples, with a smooth softmax aggregate for
//!   derivatives and the sampled-max reported alongside;
//! - [`volume_certified`] / [`volume_smooth`] — a RIGOROUS volume
//!   enclosure (sure-inside cells vs Lipschitz uncertainty band) next
//!   to a smoothed-Heaviside estimate whose lever derivative matches
//!   the Hadamard formula on fixtures.

mod primitives;
mod quotient;

pub use primitives::{
    DraftReport, EnvelopeReport, ThicknessReport, VolumeEnclosure, draft_violations,
    envelope_violation, min_thickness_soft, volume_certified, volume_smooth,
};
pub use quotient::{QuotientChart, SymmetryGroup};

use fs_constraint::{ConstraintKind, PenaltyLaw, ProofKind};
use fs_opt::Class;

/// Crate version, re-exported for provenance stamping.
pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Certificate story a primitive ships with.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CertKind {
    /// A rigorous two-sided enclosure (volume).
    Enclosure,
    /// Smooth estimate with measured/declared error (aggregations).
    SmoothEstimate,
    /// The property holds BY CONSTRUCTION (symmetry quotients).
    ExactByConstruction,
}

/// The five primitives, as an enumerable table.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeoPrimitive {
    /// Minimum wall thickness (medial oracle aggregation).
    MinThickness,
    /// Draft angle for casting/molding.
    DraftAngle,
    /// Symmetry group quotient.
    Symmetry,
    /// Bounding envelope / keep-out containment.
    Envelope,
    /// Volume / mass.
    Volume,
}

/// A primitive's declared story: differentiability, certificate, and
/// the ASCENT constraint kind it maps to.
#[derive(Debug, Clone, PartialEq)]
pub struct PrimitiveDescriptor {
    /// Which primitive.
    pub primitive: GeoPrimitive,
    /// Differentiability class of its smooth form.
    pub class: Class,
    /// Certificate story.
    pub certificate: CertKind,
    /// The fs-constraint kind ASCENT should host it under.
    pub kind: ConstraintKind,
}

impl GeoPrimitive {
    /// The declared table row (the bead's "every primitive declares"
    /// requirement, in code rather than prose).
    #[must_use]
    pub fn descriptor(self) -> PrimitiveDescriptor {
        match self {
            GeoPrimitive::MinThickness => PrimitiveDescriptor {
                primitive: self,
                class: Class::C1,
                certificate: CertKind::SmoothEstimate,
                kind: ConstraintKind::Fabrication {
                    process: "min-wall-thickness".to_string(),
                },
            },
            GeoPrimitive::DraftAngle => PrimitiveDescriptor {
                primitive: self,
                class: Class::C1,
                certificate: CertKind::SmoothEstimate,
                kind: ConstraintKind::Fabrication {
                    process: "casting-draft".to_string(),
                },
            },
            GeoPrimitive::Symmetry => PrimitiveDescriptor {
                primitive: self,
                class: Class::Smooth,
                certificate: CertKind::ExactByConstruction,
                // Structurally impossible to violate: hosted as Hard,
                // but never active.
                kind: ConstraintKind::Hard,
            },
            GeoPrimitive::Envelope => PrimitiveDescriptor {
                primitive: self,
                class: Class::C1,
                certificate: CertKind::SmoothEstimate,
                kind: ConstraintKind::Hard,
            },
            GeoPrimitive::Volume => PrimitiveDescriptor {
                primitive: self,
                class: Class::C1,
                certificate: CertKind::Enclosure,
                kind: ConstraintKind::Soft(PenaltyLaw::Quadratic { weight: 1.0 }),
            },
        }
    }

    /// Primitives whose rigorous form could escalate to a proof
    /// obligation (interval/SOS) in a certification deliverable.
    #[must_use]
    pub fn proof_escalation(self) -> Option<ProofKind> {
        match self {
            GeoPrimitive::Envelope | GeoPrimitive::Volume => Some(ProofKind::Interval),
            _ => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn version_is_stamped() {
        assert!(!VERSION.is_empty());
    }

    #[test]
    fn descriptor_table_is_total() {
        for p in [
            GeoPrimitive::MinThickness,
            GeoPrimitive::DraftAngle,
            GeoPrimitive::Symmetry,
            GeoPrimitive::Envelope,
            GeoPrimitive::Volume,
        ] {
            let d = p.descriptor();
            assert_eq!(d.primitive, p);
        }
        assert_eq!(
            GeoPrimitive::Symmetry.descriptor().certificate,
            CertKind::ExactByConstruction
        );
    }
}
