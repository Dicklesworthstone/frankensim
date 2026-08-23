//! Binding-constraint witnesses for evidence composition
//! (bead frankensim-h0vur, core algebra slice).
//!
//! [`super::color::compose`] returns the meet's VALUE and discards its
//! WITNESS: a six-stage spine's composed `Estimated` tells a user the
//! answer's colour but not WHICH stage made it so. The sentence worth
//! paying for is "Estimated BECAUSE the convection correlation is outside
//! its validated Reynolds domain; every other stage is Verified" — the
//! meet's argmin, machine-readable, not prose. This module carries it.
//!
//! Rules, matching the bead's canonical scope:
//! - The witness set is TYPED, NONEMPTY, and CANONICALLY ORDERED. Where
//!   several operands tie at the weakest colour, ALL of them are retained
//!   in canonical order — the binding constraint is a set, and silently
//!   picking one would be fresh laundering.
//! - The discriminant is machine-readable: a stable `(stage, detail)`
//!   identity pair plus a [`BindingCause`], never a message string. (The
//!   ProvenanceClass defect taught this crate's consumers what prose in a
//!   string costs; bead ead4c2f6.)
//! - Composition can land BELOW both operands (disjoint validated regimes
//!   degrade to `Estimated` with infinite dispersion). That is not either
//!   operand's fault alone: both witnesses are retained with
//!   [`BindingCause::CompositionDegraded`], because the honest explanation
//!   is the incompatibility, not one side of it.
//! - Anti-laundering by construction: witnesses only ever flow from the
//!   operands whose rank equals the composed rank (or from both, on
//!   degradation). A stage stronger than the composed colour cannot appear
//!   as its explanation.
//!
//! BOUNDARY: this is the fs-evidence layer's algebra only. Threading the
//! witness through the solve spine's stage receipts and rendering it as
//! the report headline are the owning bead's remaining scope, blocked on
//! the conduction/QoI stages executing; the e2e discrimination test named
//! there cannot be non-vacuous until then, and this module does not claim
//! it. Serialization to ledger rows is deliberately absent until the
//! spine consumer fixes the row schema.

use std::collections::BTreeSet;

use crate::color::{Color, IntervalOp, compose};

/// Bounded witness-set cardinality: a spine has single-digit stages; a set
/// approaching this cap indicates a misuse (e.g. folding a whole campaign
/// into one claim), which refuses rather than truncating silently.
pub const MAX_WITNESSES: usize = 64;

/// Machine-readable identity of one binding constraint: the stage that
/// bound the composed colour and the model/card discriminant within it.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct WitnessId {
    /// Stage identity (e.g. `"flow-network"`), stable across runs.
    pub stage: String,
    /// Finer discriminant within the stage (e.g. the correlation card id);
    /// empty is admitted for stages with a single claim source.
    pub detail: String,
}

impl core::fmt::Display for WitnessId {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        if self.detail.is_empty() {
            write!(f, "{}", self.stage)
        } else {
            write!(f, "{}:{}", self.stage, self.detail)
        }
    }
}

/// WHY this witness binds the composed colour.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum BindingCause {
    /// The witness's own claim is the (or a tied) weakest operand.
    WeakestOperand,
    /// The composition itself landed below both operands (e.g. disjoint
    /// validated regimes); both sides' witnesses carry this cause.
    CompositionDegraded,
}

/// One retained binding constraint.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct BindingWitness {
    /// Who bound it.
    pub id: WitnessId,
    /// Why it binds.
    pub cause: BindingCause,
}

/// Typed refusals of the witness algebra.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum WitnessError {
    /// The canonical witness set exceeded [`MAX_WITNESSES`].
    Cardinality {
        /// Attempted size.
        attempted: usize,
    },
    /// A leaf was declared with an empty stage identity.
    EmptyStage,
}

impl core::fmt::Display for WitnessError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            WitnessError::Cardinality { attempted } => write!(
                f,
                "FS-EVIDENCE-WITNESS-CARDINALITY: {attempted} witnesses exceed the \
                 {MAX_WITNESSES} cap; a claim folding this much is a misuse, not a set to truncate"
            ),
            WitnessError::EmptyStage => write!(
                f,
                "FS-EVIDENCE-WITNESS-EMPTY-STAGE: a witness needs a nonempty stage identity"
            ),
        }
    }
}

impl std::error::Error for WitnessError {}

/// A colour together with its nonempty, canonically ordered binding-witness
/// set. Constructed only through [`WitnessedColor::leaf`] and
/// [`compose_witnessed`], so the nonemptiness and ordering invariants hold
/// by construction.
#[derive(Debug, Clone, PartialEq)]
pub struct WitnessedColor {
    color: Color,
    witnesses: Vec<BindingWitness>,
}

impl WitnessedColor {
    /// A leaf claim: its own stage is its binding constraint.
    ///
    /// # Errors
    /// [`WitnessError::EmptyStage`].
    pub fn leaf(
        color: Color,
        stage: impl Into<String>,
        detail: impl Into<String>,
    ) -> Result<WitnessedColor, WitnessError> {
        let stage = stage.into();
        if stage.is_empty() {
            return Err(WitnessError::EmptyStage);
        }
        Ok(WitnessedColor {
            color,
            witnesses: vec![BindingWitness {
                id: WitnessId {
                    stage,
                    detail: detail.into(),
                },
                cause: BindingCause::WeakestOperand,
            }],
        })
    }

    /// The composed colour.
    #[must_use]
    pub fn color(&self) -> &Color {
        &self.color
    }

    /// The binding-constraint set: nonempty, canonical `(stage, detail)`
    /// order, duplicate-free.
    #[must_use]
    pub fn witnesses(&self) -> &[BindingWitness] {
        &self.witnesses
    }
}

fn canonicalize(mut witnesses: Vec<BindingWitness>) -> Result<Vec<BindingWitness>, WitnessError> {
    // Dedup by identity, keeping the STRONGER cause claim out of the way:
    // if the same id appears as both WeakestOperand and CompositionDegraded
    // (possible across fold steps), retain CompositionDegraded — it is the
    // more specific explanation and never weaker.
    let mut by_id: BTreeSet<(WitnessId, BindingCause)> = BTreeSet::new();
    witnesses.sort();
    let mut out: Vec<BindingWitness> = Vec::new();
    for witness in witnesses {
        if by_id.iter().any(|(id, _)| *id == witness.id) {
            if witness.cause == BindingCause::CompositionDegraded
                && let Some(existing) = out.iter_mut().find(|w| w.id == witness.id)
            {
                existing.cause = BindingCause::CompositionDegraded;
            }
            continue;
        }
        by_id.insert((witness.id.clone(), witness.cause));
        out.push(witness);
    }
    if out.len() > MAX_WITNESSES {
        return Err(WitnessError::Cardinality {
            attempted: out.len(),
        });
    }
    Ok(out)
}

/// Compose two witnessed colours: the value is exactly
/// [`compose`]`(a, b, op)`, and the witness set is the meet's argmin —
/// the operands whose rank equals the composed rank, or BOTH operands
/// (as [`BindingCause::CompositionDegraded`]) when the composition lands
/// strictly below both.
///
/// # Errors
/// [`WitnessError::Cardinality`].
pub fn compose_witnessed(
    a: &WitnessedColor,
    b: &WitnessedColor,
    op: IntervalOp,
) -> Result<WitnessedColor, WitnessError> {
    let color = compose(&a.color, &b.color, op);
    let composed_rank = color.rank();
    let degraded = composed_rank < a.color.rank().min(b.color.rank());
    let mut witnesses = Vec::new();
    if degraded {
        for side in [a, b] {
            witnesses.extend(side.witnesses.iter().map(|w| BindingWitness {
                id: w.id.clone(),
                cause: BindingCause::CompositionDegraded,
            }));
        }
    } else {
        for side in [a, b] {
            if side.color.rank() == composed_rank {
                witnesses.extend(side.witnesses.iter().cloned());
            }
        }
    }
    Ok(WitnessedColor {
        color,
        witnesses: canonicalize(witnesses)?,
    })
}

/// Left fold of [`compose_witnessed`] over a nonempty sequence.
///
/// # Errors
/// [`WitnessError::Cardinality`]; `None`-like misuse (empty input) is a
/// caller bug surfaced as `EmptyStage` to stay in this error family.
pub fn compose_all_witnessed(
    items: &[WitnessedColor],
    op: IntervalOp,
) -> Result<WitnessedColor, WitnessError> {
    let Some((first, rest)) = items.split_first() else {
        return Err(WitnessError::EmptyStage);
    };
    let mut acc = first.clone();
    for item in rest {
        acc = compose_witnessed(&acc, item, op)?;
    }
    Ok(acc)
}
