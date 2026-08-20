//! AeroEffectOwners record contract (bead wf-root-guzez.5.6.2, E4.3-ii).
//! Plan §5.2 (Round-2 revision): every physical aero effect has EXACTLY
//! ONE owner — the record makes double-counting (the classic added-mass /
//! circulatory overlap) and orphaned effects ADMISSION-TIME refusals
//! instead of silent physics bugs. The record's content digest is a
//! ModelId ingredient, canonical in effect order (assignment list order
//! must not move the digest).
//!
//! The load-bearing class rule: the NONCIRCULATORY owner must be of class
//! `AddedMassOnly` — the indicial section channels (fs-airfoil::unsteady)
//! carry no apparent-mass term, and this admission is what makes that
//! division checkable rather than folklore.

use crate::Refusal;
use fs_blake3::hash_domain;

/// The six owned physical effects (plan §5.2 effect-ownership graph).
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum AeroEffect {
    /// Circulatory response to section motion (Wagner class).
    MotionCirculatory,
    /// Circulatory response to incident gusts (Küssner class).
    IncidentGust,
    /// Apparent-mass (noncirculatory) reaction.
    Noncirculatory,
    /// Trailing-edge separation / dynamic stall state.
    Separation,
    /// Finite-span 3-D induction.
    Induction3d,
    /// Far-wake memory beyond the bound system.
    FarWake,
}

/// All effects in canonical (digest) order.
pub const ALL_EFFECTS: [AeroEffect; 6] = [
    AeroEffect::MotionCirculatory,
    AeroEffect::IncidentGust,
    AeroEffect::Noncirculatory,
    AeroEffect::Separation,
    AeroEffect::Induction3d,
    AeroEffect::FarWake,
];

/// Owner component class (the admission's type system).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OwnerClass {
    /// Indicial kernel state machinery (fs-airfoil::unsteady).
    IndicialKernel,
    /// Generalized added mass ONLY — no circulatory content.
    AddedMassOnly,
    /// Bound lifting-surface solve (fs-wing).
    LiftingSurface,
    /// Lagged separation coordinate.
    SeparationLag,
    /// Wake model (trailing system now; fs-vpm hybrid at E4.7).
    WakeModel,
}

impl OwnerClass {
    fn discriminant(self) -> u8 {
        match self {
            OwnerClass::IndicialKernel => 0,
            OwnerClass::AddedMassOnly => 1,
            OwnerClass::LiftingSurface => 2,
            OwnerClass::SeparationLag => 3,
            OwnerClass::WakeModel => 4,
        }
    }
}

/// One effect→owner assignment.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct OwnerAssignment {
    /// The owned effect.
    pub effect: AeroEffect,
    /// Registered owner component id.
    pub owner_id: &'static str,
    /// Owner class (checked for the noncirculatory rule).
    pub class: OwnerClass,
}

/// The full record. Build, then [`AeroEffectOwners::admit`] — every
/// downstream consumer takes the [`AdmittedOwners`] token, so an
/// unadmitted record cannot flow into a model identity.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AeroEffectOwners {
    /// Assignment list (any order; admission canonicalizes).
    pub assignments: Vec<OwnerAssignment>,
}

/// Proof-of-admission wrapper: the only path to the digest.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct AdmittedOwners {
    canonical: Vec<OwnerAssignment>,
}

/// The registered Wright v1 record.
#[must_use]
pub fn wright_owners_v1() -> AeroEffectOwners {
    AeroEffectOwners {
        assignments: vec![
            OwnerAssignment {
                effect: AeroEffect::MotionCirculatory,
                owner_id: "fs-airfoil.unsteady.wagner-jones-2pole-v1",
                class: OwnerClass::IndicialKernel,
            },
            OwnerAssignment {
                effect: AeroEffect::IncidentGust,
                owner_id: "fs-airfoil.unsteady.kussner-2pole-v1",
                class: OwnerClass::IndicialKernel,
            },
            OwnerAssignment {
                effect: AeroEffect::Noncirculatory,
                owner_id: "fs-flyer.addedmass.joint-6nc-v1",
                class: OwnerClass::AddedMassOnly,
            },
            OwnerAssignment {
                effect: AeroEffect::Separation,
                owner_id: "fs-airfoil.unsteady.separation-lag-kirchhoff-v1",
                class: OwnerClass::SeparationLag,
            },
            OwnerAssignment {
                effect: AeroEffect::Induction3d,
                owner_id: "fs-wing.weissinger-l-multisurface-v1",
                class: OwnerClass::LiftingSurface,
            },
            OwnerAssignment {
                effect: AeroEffect::FarWake,
                owner_id: "fs-wing.trailing-horseshoe-far-v1",
                class: OwnerClass::WakeModel,
            },
        ],
    }
}

impl AeroEffectOwners {
    /// Admit the record: every effect owned EXACTLY once, and the
    /// noncirculatory owner is `AddedMassOnly`.
    ///
    /// # Errors
    /// `effect-owner-missing`, `effect-owner-duplicate`,
    /// `noncirculatory-owner-not-added-mass`, `owner-id-empty`.
    pub fn admit(&self) -> Result<AdmittedOwners, Refusal> {
        for a in &self.assignments {
            if a.owner_id.is_empty() {
                return Err(Refusal {
                    code: "owner-id-empty",
                    message: format!("effect {:?} has an empty owner id", a.effect),
                    ranked_repairs: vec!["use a registered component id".into()],
                });
            }
        }
        let mut canonical = Vec::with_capacity(ALL_EFFECTS.len());
        for effect in ALL_EFFECTS {
            let owners: Vec<&OwnerAssignment> = self
                .assignments
                .iter()
                .filter(|a| a.effect == effect)
                .collect();
            match owners.len() {
                0 => {
                    return Err(Refusal {
                        code: "effect-owner-missing",
                        message: format!("effect {effect:?} has NO owner — orphaned physics"),
                        ranked_repairs: vec![
                            "assign exactly one registered owner per effect".into(),
                            "if the effect is out of scope, that is a MODEL change, not an \
                             omission — record it in the model identity"
                                .into(),
                        ],
                    });
                }
                1 => canonical.push(*owners[0]),
                n => {
                    return Err(Refusal {
                        code: "effect-owner-duplicate",
                        message: format!(
                            "effect {effect:?} claimed by {n} owners ({:?}) — double-counted \
                             physics",
                            owners.iter().map(|o| o.owner_id).collect::<Vec<_>>()
                        ),
                        ranked_repairs: vec![
                            "one owner per effect; the others must consume, not re-add".into(),
                        ],
                    });
                }
            }
        }
        let nc = canonical
            .iter()
            .find(|a| a.effect == AeroEffect::Noncirculatory)
            .expect("canonical covers all effects");
        if nc.class != OwnerClass::AddedMassOnly {
            return Err(Refusal {
                code: "noncirculatory-owner-not-added-mass",
                message: format!(
                    "noncirculatory owner {} has class {:?} — only AddedMassOnly may own \
                     apparent mass (the double-count guard)",
                    nc.owner_id, nc.class
                ),
                ranked_repairs: vec![
                    "route apparent mass through the generalized added-mass solve".into(),
                ],
            });
        }
        Ok(AdmittedOwners { canonical })
    }
}

impl AdmittedOwners {
    /// Canonical assignments (effect order).
    #[must_use]
    pub fn assignments(&self) -> &[OwnerAssignment] {
        &self.canonical
    }

    /// Content digest (ModelId ingredient); canonical in effect order, so
    /// assignment-list permutations do not move it.
    #[must_use]
    pub fn digest(&self) -> String {
        let mut p = Vec::new();
        for (i, a) in self.canonical.iter().enumerate() {
            p.push(u8::try_from(i).expect("6 effects"));
            p.extend_from_slice(
                &u32::try_from(a.owner_id.len())
                    .expect("short id")
                    .to_le_bytes(),
            );
            p.extend_from_slice(a.owner_id.as_bytes());
            p.push(a.class.discriminant());
        }
        hash_domain("org.frankensim.fs-flyer.aero-effect-owners.v1", &p).to_hex()
    }
}
