//! Material and constitutive-model cards (bead 5hmy, PR-2 of 5).
//!
//! A [`MaterialCard`] is the immutable, content-addressed identity of a
//! NAMED MATERIAL STATE — chemistry + phase + temper/process + revision
//! — carrying its property claims, its constitutive model cards, and
//! explicit supersedes lineage. Revision is never an edit: a successor
//! card links its predecessor's content hash, and both remain
//! retrievable forever.
//!
//! A [`ConstitutiveModelCard`] names a LAW (id + version), its
//! canonical parameter block, the state schema its internal variables
//! follow, how initial state is obtained, and where its parameters are
//! valid. The card stores DATA about the law — the executable law-node
//! protocol lives in L3 fs-material (bead kagp), never here.

use std::collections::BTreeMap;
use std::fmt;

use fs_blake3::{ContentHash, hash_domain};
use fs_evidence::ValidityDomain;
use fs_qty::Dims;

use crate::{ClaimId, ClaimSet, MatDbError, PropertyClaim, Provenance};

/// Hash domain for constitutive-model-card canonical identity.
const MODEL_HASH_DOMAIN: &str = "org.frankensim.fs-matdb.constitutive-model-card.v1";
/// Semantic version of the canonical constitutive parameter-block identity.
pub const CANONICAL_PARAMETER_BLOCK_IDENTITY_VERSION: u32 = 1;
/// BLAKE3 domain for the canonical constitutive parameter-block identity.
pub const CANONICAL_PARAMETER_BLOCK_IDENTITY_DOMAIN: &str =
    "org.frankensim.fs-matdb.canonical-parameter-block.v1";
/// Hash domain for material-card canonical identity.
const MATERIAL_HASH_DOMAIN: &str = "org.frankensim.fs-matdb.material-card.v1";

/// The fs-matdb card schema version (bumped only with a migration note
/// in CONTRACT.md).
pub const MATDB_SCHEMA_VERSION: u32 = 1;

/// A constitutive law's stable identity (e.g. `"j2-plasticity"`,
/// `"neo-hookean"`, `"norton-bailey-creep"`). Free-form name; the
/// (id, version) pair is what parameter blocks bind to.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct LawId(pub String);

/// How a law's internal state vector starts. Data only — the L3
/// executable protocol interprets it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InitialStatePolicy {
    /// All internal variables start at zero (virgin material).
    ZeroInternalState,
    /// The consumer MUST supply an explicit initial state; the card
    /// refuses to imply one.
    RequiresDeclaredState,
}

/// One named, dimensioned law parameter.
#[derive(Debug, Clone, PartialEq)]
pub struct LawParameter {
    /// SI base-unit value (finite).
    pub value: f64,
    /// The parameter's dimensions.
    pub dims: Dims,
}

/// A constitutive model card: law identity, canonical parameter block,
/// state schema, initial-state policy, validity, and source hashes.
#[derive(Debug, Clone, PartialEq)]
pub struct ConstitutiveModelCard {
    /// The law this parameter block instantiates.
    pub law: LawId,
    /// The law's semantic version (parameter meaning may change across
    /// versions; the pair is the binding identity).
    pub law_version: u32,
    /// Canonical parameter block, keyed by parameter name (BTreeMap =
    /// one canonical order for hashing).
    pub parameters: BTreeMap<String, LawParameter>,
    /// The internal-state schema version the law's variables follow.
    pub state_schema_version: u32,
    /// How initial state is obtained.
    pub initial_state: InitialStatePolicy,
    /// Where the parameter block is valid (THE shared validity type).
    pub validity: ValidityDomain,
    /// Content hashes of the calibration/source artifacts.
    pub sources: Vec<ContentHash>,
    /// Where the card came from and under what license.
    pub provenance: Provenance,
}

impl ConstitutiveModelCard {
    /// Validate the card's structural gates (finite parameters, usable
    /// validity, load-bearing provenance, nonempty parameter block).
    ///
    /// # Errors
    /// Typed [`MatDbError`] refusals; nothing partial.
    pub fn validate(&self) -> Result<(), MatDbError> {
        self.provenance.validate()?;
        if self.parameters.is_empty() {
            return Err(MatDbError::EmptyParameterBlock {
                law: self.law.clone(),
            });
        }
        for (name, parameter) in &self.parameters {
            if !parameter.value.is_finite() {
                return Err(MatDbError::NonFiniteParameter {
                    law: self.law.clone(),
                    parameter: name.clone(),
                    bits: parameter.value.to_bits(),
                });
            }
        }
        for (axis, &(lo, hi)) in self.validity.bounds() {
            if lo.is_nan() || hi.is_nan() {
                return Err(MatDbError::UnusableValidity { axis: axis.clone() });
            }
        }
        Ok(())
    }

    /// Canonical identity of only the ordered, dimensioned parameter block.
    ///
    /// The law id, law version, state schema, validity, sources, and
    /// provenance deliberately do not enter this digest. Consumers must bind
    /// this hash alongside those separate semantic fields rather than treating
    /// a parameter block as a complete model identity.
    ///
    /// # Errors
    /// Refuses to mint authority for a structurally invalid model card.
    pub fn canonical_parameters_hash(&self) -> Result<ContentHash, MatDbError> {
        self.validate()?;
        Ok(self.canonical_parameters_hash_with_schema(
            CANONICAL_PARAMETER_BLOCK_IDENTITY_VERSION,
            CANONICAL_PARAMETER_BLOCK_IDENTITY_DOMAIN,
        ))
    }

    fn canonical_parameters_hash_with_schema(&self, version: u32, domain: &str) -> ContentHash {
        let mut payload = Vec::new();
        payload.extend_from_slice(&version.to_le_bytes());
        payload.extend_from_slice(
            &u64::try_from(self.parameters.len())
                .unwrap_or(u64::MAX)
                .to_le_bytes(),
        );
        let mut push = |part: &[u8]| {
            payload.extend_from_slice(&u64::try_from(part.len()).unwrap_or(u64::MAX).to_le_bytes());
            payload.extend_from_slice(part);
        };
        for (name, parameter) in &self.parameters {
            push(name.as_bytes());
            push(&parameter.value.to_bits().to_le_bytes());
            push(&crate::dims_bytes(parameter.dims));
        }
        hash_domain(domain, &payload)
    }

    /// Canonical content identity over every semantic field.
    #[must_use]
    pub fn content_hash(&self) -> ContentHash {
        let mut payload = Vec::new();
        let mut push = |part: &[u8]| {
            payload.extend_from_slice(&u64::try_from(part.len()).unwrap_or(u64::MAX).to_le_bytes());
            payload.extend_from_slice(part);
        };
        push(self.law.0.as_bytes());
        push(&self.law_version.to_le_bytes());
        for (name, parameter) in &self.parameters {
            push(name.as_bytes());
            push(&parameter.value.to_bits().to_le_bytes());
            push(&crate::dims_bytes(parameter.dims));
        }
        push(&self.state_schema_version.to_le_bytes());
        push(match self.initial_state {
            InitialStatePolicy::ZeroInternalState => b"zero-internal-state".as_slice(),
            InitialStatePolicy::RequiresDeclaredState => b"requires-declared-state",
        });
        for (axis, &(lo, hi)) in self.validity.bounds() {
            push(axis.as_bytes());
            push(&lo.to_bits().to_le_bytes());
            push(&hi.to_bits().to_le_bytes());
        }
        for source in &self.sources {
            push(&source.0);
        }
        push(self.provenance.source.as_bytes());
        push(self.provenance.license.as_bytes());
        if let Some(artifact) = &self.provenance.artifact {
            push(&artifact.0);
        }
        hash_domain(MODEL_HASH_DOMAIN, &payload)
    }
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use fs_evidence::ValidityDomain;

    use super::*;
    use crate::Provenance;

    fn model() -> ConstitutiveModelCard {
        let mut parameters = BTreeMap::new();
        parameters.insert(
            "p".to_string(),
            LawParameter {
                value: 1.0,
                dims: Dims([0; 6]),
            },
        );
        ConstitutiveModelCard {
            law: LawId("identity-test".to_string()),
            law_version: 1,
            parameters,
            state_schema_version: 0,
            initial_state: InitialStatePolicy::ZeroInternalState,
            validity: ValidityDomain::unconstrained(),
            sources: Vec::new(),
            provenance: Provenance {
                source: "identity unit test".to_string(),
                license: "test-only".to_string(),
                artifact: None,
            },
        }
    }

    #[test]
    fn canonical_parameter_identity_schema_moves_version_and_domain() {
        let model = model();
        let canonical = model
            .canonical_parameters_hash()
            .expect("valid model mints parameter identity");
        assert_ne!(
            canonical,
            model.canonical_parameters_hash_with_schema(
                CANONICAL_PARAMETER_BLOCK_IDENTITY_VERSION + 1,
                CANONICAL_PARAMETER_BLOCK_IDENTITY_DOMAIN,
            )
        );
        assert_ne!(
            canonical,
            model.canonical_parameters_hash_with_schema(
                CANONICAL_PARAMETER_BLOCK_IDENTITY_VERSION,
                "org.frankensim.fs-matdb.canonical-parameter-block.foreign",
            )
        );
    }
}

/// The identity of a NAMED MATERIAL STATE. "AA6061" is not a material
/// state; "AA6061, wrought FCC matrix, T6, revision 3" is.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct MaterialStateId {
    /// Chemistry / grade designation (e.g. "AA6061", "PA12+CF15").
    pub chemistry: String,
    /// Phase / microstructural family (e.g. "wrought", "as-printed").
    pub phase: String,
    /// Temper / process state (e.g. "T6", "annealed-2h-350C").
    pub process: String,
    /// Revision of THIS named state's data card. Advances only through
    /// [`MaterialCard::supersede`].
    pub revision: u32,
}

impl fmt::Display for MaterialStateId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}/{}/{} rev {}",
            self.chemistry, self.phase, self.process, self.revision
        )
    }
}

/// The immutable card for one material state: claims, models, lineage.
/// Constructed ONLY by [`MaterialCard::assemble`] (revision 0) or
/// [`MaterialCard::supersede`] (revision n+1 linking the predecessor's
/// content hash) — both validate everything; there is no mutable
/// access afterward.
#[derive(Debug, Clone, PartialEq)]
pub struct MaterialCard {
    id: MaterialStateId,
    schema_version: u32,
    supersedes: Option<ContentHash>,
    claims: ClaimSet,
    models: Vec<ConstitutiveModelCard>,
}

impl MaterialCard {
    /// Assemble a REVISION-0 card from a claim set and model cards.
    ///
    /// # Errors
    /// [`MatDbError::RevisionNotZero`] when the id claims a nonzero
    /// revision (lineage must start at 0 and only `supersede` advances
    /// it), plus every model-card validation refusal.
    pub fn assemble(
        id: MaterialStateId,
        claims: ClaimSet,
        models: Vec<ConstitutiveModelCard>,
    ) -> Result<MaterialCard, MatDbError> {
        if id.revision != 0 {
            return Err(MatDbError::RevisionNotZero {
                offered: id.revision,
            });
        }
        for model in &models {
            model.validate()?;
        }
        Ok(MaterialCard {
            id,
            schema_version: MATDB_SCHEMA_VERSION,
            supersedes: None,
            claims,
            models,
        })
    }

    /// Build the successor card: same chemistry/phase/process, revision
    /// exactly one higher, supersedes link bound to the predecessor's
    /// content hash. The predecessor is untouched and stays valid.
    ///
    /// # Errors
    /// [`MatDbError::SupersedesMismatch`] when the named state differs;
    /// model validation refusals as in [`MaterialCard::assemble`].
    pub fn supersede(
        predecessor: &MaterialCard,
        claims: ClaimSet,
        models: Vec<ConstitutiveModelCard>,
    ) -> Result<MaterialCard, MatDbError> {
        for model in &models {
            model.validate()?;
        }
        let id = MaterialStateId {
            revision: predecessor.id.revision.checked_add(1).ok_or(
                MatDbError::SupersedesMismatch {
                    reason: "revision counter exhausted",
                },
            )?,
            ..predecessor.id.clone()
        };
        Ok(MaterialCard {
            id,
            schema_version: MATDB_SCHEMA_VERSION,
            supersedes: Some(predecessor.content_hash()),
            claims,
            models,
        })
    }

    /// The named material state this card describes.
    #[must_use]
    pub fn id(&self) -> &MaterialStateId {
        &self.id
    }

    /// The card schema version.
    #[must_use]
    pub fn schema_version(&self) -> u32 {
        self.schema_version
    }

    /// The predecessor card's content hash, when this is a successor.
    #[must_use]
    pub fn supersedes(&self) -> Option<ContentHash> {
        self.supersedes
    }

    /// The card's claims (read-only; the by-key index lives on
    /// [`ClaimSet::claims_for`]).
    #[must_use]
    pub fn claims(&self) -> &ClaimSet {
        &self.claims
    }

    /// EVERY claim for a property name, in insertion order (the by-key
    /// index the bead requires, delegated to the claim set).
    #[must_use]
    pub fn claims_for(&self, name: &str) -> Vec<(ClaimId, &PropertyClaim)> {
        self.claims.claims_for(name)
    }

    /// The constitutive model cards.
    #[must_use]
    pub fn models(&self) -> &[ConstitutiveModelCard] {
        &self.models
    }

    /// Model cards instantiating one law id, in card order.
    #[must_use]
    pub fn models_for(&self, law: &LawId) -> Vec<&ConstitutiveModelCard> {
        self.models.iter().filter(|m| &m.law == law).collect()
    }

    /// Canonical content identity: id, schema version, supersedes link,
    /// every claim id (content-ordered), every observation id, every
    /// model-card hash. Claim/observation CONTENT is already
    /// content-addressed, so hashing the id sets binds the full
    /// transitive content.
    #[must_use]
    pub fn content_hash(&self) -> ContentHash {
        let mut payload = Vec::new();
        let mut push = |part: &[u8]| {
            payload.extend_from_slice(&u64::try_from(part.len()).unwrap_or(u64::MAX).to_le_bytes());
            payload.extend_from_slice(part);
        };
        push(self.id.chemistry.as_bytes());
        push(self.id.phase.as_bytes());
        push(self.id.process.as_bytes());
        push(&self.id.revision.to_le_bytes());
        push(&self.schema_version.to_le_bytes());
        match &self.supersedes {
            Some(hash) => push(&hash.0),
            None => push(b"genesis"),
        }
        for (claim_id, _) in self.claims.claims_ordered() {
            push(&claim_id.0.0);
        }
        for observation_id in self.claims.observation_ids() {
            push(&observation_id.0.0);
        }
        for model in &self.models {
            push(&model.content_hash().0);
        }
        hash_domain(MATERIAL_HASH_DOMAIN, &payload)
    }
}
