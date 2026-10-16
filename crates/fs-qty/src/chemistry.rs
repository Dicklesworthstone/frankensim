//! Exact chemistry bookkeeping for the six-base quantity system.
//!
//! Elemental composition, stoichiometry, and integer charge are immutable
//! artifacts with canonical identifier ordering. Conservation is proved with
//! checked `i128` arithmetic: a successful [`verify_conservation`] call is an
//! exact G0 proof of `A * N = 0` and `z^T * N = 0`, not a floating tolerance.
//! Identifiers are opaque canonical labels: this module does not claim periodic-
//! table membership, chemical-formula parsing, or bracket-balance validation.
//! A conservation certificate proves only those two bookkeeping identities;
//! zero reaction columns are admitted, and no kinetic, thermodynamic, phase, or
//! chemical-meaningfulness claim follows from conservation alone.

use crate::{Amount, Mass, MolarMass};
use core::fmt;
use fs_blake3::{ContentHash, hash_domain};

const ELEMENTAL_MATRIX_ID_DOMAIN: &str = "frankensim.fs-qty.elemental-matrix.v1";
const STOICHIOMETRIC_MATRIX_ID_DOMAIN: &str = "frankensim.fs-qty.stoichiometric-matrix.v1";
const CHARGE_VECTOR_ID_DOMAIN: &str = "frankensim.fs-qty.charge-vector.v1";
const CONSERVATION_CERTIFICATE_ID_DOMAIN: &str = "frankensim.fs-qty.chemistry-conservation.v1";

/// Identifier category used by structured diagnostics.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum IdKind {
    /// Chemical species identifier.
    Species,
    /// Element symbol.
    Element,
    /// Reaction identifier.
    Reaction,
}

/// Artifact category used by shape and basis diagnostics.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ArtifactKind {
    /// Element-by-species matrix `A`.
    ElementalMatrix,
    /// Species-by-reaction matrix `N`.
    StoichiometricMatrix,
    /// Species charge vector `z`.
    ChargeVector,
}

/// Named matrix axis.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum AxisKind {
    /// Element rows.
    Elements,
    /// Species rows or columns.
    Species,
    /// Reaction columns.
    Reactions,
}

/// Conservation law whose exact arithmetic overflowed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ConservationLaw {
    /// Elemental balance `A * N = 0`.
    Elements,
    /// Charge balance `z^T * N = 0`.
    Charge,
}

/// Source quantity used to establish a mass/amount basis record.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BasisSource {
    /// Amount was derived from a supplied mass and molar mass.
    Mass,
    /// Mass was derived from a supplied amount and molar mass.
    Amount,
}

/// Quantity field named by a basis-validation error.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum BasisField {
    /// Mass value.
    Mass,
    /// Molar-mass value.
    MolarMass,
    /// Amount-of-substance value.
    Amount,
    /// Amount derived from mass and molar mass.
    DerivedAmount,
    /// Mass derived from amount and molar mass.
    DerivedMass,
}

/// Validated, byte-stable chemical species identifier.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SpeciesId(String);

impl SpeciesId {
    /// Construct an ASCII species id.
    ///
    /// The first byte must be alphanumeric. Remaining bytes may additionally
    /// use `_`, `+`, `-`, `.`, parentheses, or square brackets, which covers
    /// compact formula/charge labels without admitting whitespace or Unicode
    /// normalization ambiguity.
    ///
    /// # Errors
    /// Returns [`ChemistryError::MalformedId`] for an empty, overlong, or
    /// non-canonical identifier.
    pub fn new(value: impl Into<String>) -> Result<Self, ChemistryError> {
        let value = value.into();
        validate_general_id(IdKind::Species, &value)?;
        Ok(Self(value))
    }

    /// Borrow the canonical identifier text.
    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl fmt::Display for SpeciesId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

/// Validated element symbol.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ElementId(String);

impl ElementId {
    /// Construct a canonical one-to-three-letter element symbol.
    ///
    /// # Errors
    /// Returns [`ChemistryError::MalformedId`] unless the symbol begins with
    /// one ASCII uppercase letter followed by at most two lowercase letters.
    pub fn new(value: impl Into<String>) -> Result<Self, ChemistryError> {
        let value = value.into();
        let bytes = value.as_bytes();
        if !(1..=3).contains(&bytes.len())
            || !bytes[0].is_ascii_uppercase()
            || !bytes[1..].iter().all(u8::is_ascii_lowercase)
        {
            return Err(ChemistryError::MalformedId {
                kind: IdKind::Element,
                value,
                reason: "expected one uppercase ASCII letter followed by at most two lowercase letters",
            });
        }
        Ok(Self(value))
    }

    /// Borrow the canonical element symbol.
    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl fmt::Display for ElementId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

/// Validated reaction identifier.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ReactionId(String);

impl ReactionId {
    /// Construct a deterministic reaction id under the same compact ASCII
    /// grammar as [`SpeciesId`].
    ///
    /// # Errors
    /// Returns [`ChemistryError::MalformedId`] for malformed text.
    pub fn new(value: impl Into<String>) -> Result<Self, ChemistryError> {
        let value = value.into();
        validate_general_id(IdKind::Reaction, &value)?;
        Ok(Self(value))
    }

    /// Borrow the canonical identifier text.
    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

impl fmt::Display for ReactionId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.0)
    }
}

fn validate_general_id(kind: IdKind, value: &str) -> Result<(), ChemistryError> {
    if value.is_empty() {
        return Err(ChemistryError::MalformedId {
            kind,
            value: value.to_string(),
            reason: "identifier must not be empty",
        });
    }
    if value.len() > 128 {
        return Err(ChemistryError::MalformedId {
            kind,
            value: value.to_string(),
            reason: "identifier exceeds the 128-byte limit",
        });
    }
    let bytes = value.as_bytes();
    if !bytes[0].is_ascii_alphanumeric() {
        return Err(ChemistryError::MalformedId {
            kind,
            value: value.to_string(),
            reason: "identifier must begin with an ASCII letter or digit",
        });
    }
    if !bytes.iter().all(|byte| {
        byte.is_ascii_alphanumeric()
            || matches!(byte, b'_' | b'+' | b'-' | b'.' | b'(' | b')' | b'[' | b']')
    }) {
        return Err(ChemistryError::MalformedId {
            kind,
            value: value.to_string(),
            reason: "identifier contains whitespace, non-ASCII, or unsupported punctuation",
        });
    }
    Ok(())
}

/// Typed construction or exact-conservation failure.
#[derive(Debug, Clone, PartialEq)]
pub enum ChemistryError {
    /// Identifier text is not canonical.
    MalformedId {
        /// Identifier category.
        kind: IdKind,
        /// Rejected text.
        value: String,
        /// Stable guidance.
        reason: &'static str,
    },
    /// An axis contains the same identifier more than once.
    DuplicateId {
        /// Identifier category.
        kind: IdKind,
        /// Duplicate canonical text.
        value: String,
    },
    /// A required axis is empty.
    EmptyAxis {
        /// Artifact under construction.
        artifact: ArtifactKind,
        /// Empty axis.
        axis: AxisKind,
    },
    /// Matrix outer row count does not match its row identifiers.
    RowCountMismatch {
        /// Artifact under construction.
        artifact: ArtifactKind,
        /// Required row count.
        expected: usize,
        /// Supplied row count.
        actual: usize,
    },
    /// A matrix row has the wrong number of columns.
    ColumnCountMismatch {
        /// Artifact under construction.
        artifact: ArtifactKind,
        /// Zero-based input row.
        row: usize,
        /// Required column count.
        expected: usize,
        /// Supplied column count.
        actual: usize,
    },
    /// A vector length does not match its species axis.
    VectorLengthMismatch {
        /// Artifact under construction.
        artifact: ArtifactKind,
        /// Required length.
        expected: usize,
        /// Supplied length.
        actual: usize,
    },
    /// An elemental count is negative.
    NegativeElementCount {
        /// Element row.
        element: ElementId,
        /// Species column.
        species: SpeciesId,
        /// Rejected count.
        count: i128,
    },
    /// Canonical species axes disagree between artifacts.
    SpeciesBasisMismatch {
        /// Artifact whose axis failed to match the elemental matrix.
        artifact: ArtifactKind,
        /// Expected canonical species ids.
        expected: Vec<SpeciesId>,
        /// Actual canonical species ids.
        actual: Vec<SpeciesId>,
    },
    /// Exact `i128` multiplication or accumulation overflowed.
    ArithmeticOverflow {
        /// Law being checked.
        law: ConservationLaw,
        /// Element row for an elemental check.
        element: Option<ElementId>,
        /// Reaction column.
        reaction: ReactionId,
        /// Species term whose multiplication/accumulation overflowed.
        species: SpeciesId,
    },
    /// `A * N` has a nonzero entry.
    ElementImbalance {
        /// Unbalanced element.
        element: ElementId,
        /// Unbalanced reaction.
        reaction: ReactionId,
        /// Exact signed residual.
        residual: i128,
    },
    /// `z^T * N` has a nonzero entry.
    ChargeImbalance {
        /// Unbalanced reaction.
        reaction: ReactionId,
        /// Exact signed residual in elementary-charge units.
        residual: i128,
    },
    /// Typed mass/amount basis input is invalid.
    InvalidBasisValue {
        /// Field that failed validation.
        field: BasisField,
        /// Rejected or derived coherent-SI value.
        value: f64,
        /// Stable guidance.
        reason: &'static str,
    },
}

impl fmt::Display for ChemistryError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MalformedId {
                kind,
                value,
                reason,
            } => write!(f, "malformed {kind:?} id {value:?}: {reason}"),
            Self::DuplicateId { kind, value } => {
                write!(f, "duplicate {kind:?} id {value:?}")
            }
            Self::EmptyAxis { artifact, axis } => {
                write!(f, "{artifact:?} requires a nonempty {axis:?} axis")
            }
            Self::RowCountMismatch {
                artifact,
                expected,
                actual,
            } => write!(
                f,
                "{artifact:?} row count mismatch: expected {expected}, got {actual}"
            ),
            Self::ColumnCountMismatch {
                artifact,
                row,
                expected,
                actual,
            } => write!(
                f,
                "{artifact:?} row {row} column count mismatch: expected {expected}, got {actual}"
            ),
            Self::VectorLengthMismatch {
                artifact,
                expected,
                actual,
            } => write!(
                f,
                "{artifact:?} length mismatch: expected {expected}, got {actual}"
            ),
            Self::NegativeElementCount {
                element,
                species,
                count,
            } => write!(
                f,
                "elemental matrix count A[{element},{species}] must be nonnegative, got {count}"
            ),
            Self::SpeciesBasisMismatch {
                artifact,
                expected,
                actual,
            } => write!(
                f,
                "{artifact:?} species basis {actual:?} does not match elemental basis {expected:?}"
            ),
            Self::ArithmeticOverflow {
                law,
                element,
                reaction,
                species,
            } => write!(
                f,
                "checked {law:?} conservation arithmetic overflowed at element {element:?}, reaction {reaction}, species {species}"
            ),
            Self::ElementImbalance {
                element,
                reaction,
                residual,
            } => write!(
                f,
                "element {element} is imbalanced in reaction {reaction}: exact residual {residual}"
            ),
            Self::ChargeImbalance { reaction, residual } => write!(
                f,
                "charge is imbalanced in reaction {reaction}: exact residual {residual}"
            ),
            Self::InvalidBasisValue {
                field,
                value,
                reason,
            } => write!(f, "invalid {field:?} basis value {value}: {reason}"),
        }
    }
}

impl core::error::Error for ChemistryError {}

fn canonical_order<T: Ord>(ids: &[T]) -> Result<Vec<usize>, usize> {
    let mut order: Vec<usize> = (0..ids.len()).collect();
    order.sort_by(|left, right| ids[*left].cmp(&ids[*right]));
    for pair in order.windows(2) {
        if ids[pair[0]] == ids[pair[1]] {
            return Err(pair[1]);
        }
    }
    Ok(order)
}

fn reorder_owned<T>(values: Vec<T>, order: &[usize]) -> Vec<T> {
    debug_assert_eq!(values.len(), order.len());
    let mut slots: Vec<Option<T>> = values.into_iter().map(Some).collect();
    order
        .iter()
        .map(|&source| {
            slots[source]
                .take()
                .expect("canonical order must contain each source index exactly once")
        })
        .collect()
}

fn push_usize(bytes: &mut Vec<u8>, value: usize) {
    bytes.extend_from_slice(&(value as u128).to_le_bytes());
}

fn push_text(bytes: &mut Vec<u8>, value: &str) {
    push_usize(bytes, value.len());
    bytes.extend_from_slice(value.as_bytes());
}

fn push_species(bytes: &mut Vec<u8>, species: &[SpeciesId]) {
    push_usize(bytes, species.len());
    for id in species {
        push_text(bytes, id.as_str());
    }
}

/// Immutable canonical element-by-species matrix `A`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ElementalMatrix {
    elements: Vec<ElementId>,
    species: Vec<SpeciesId>,
    values: Vec<i128>,
}

impl ElementalMatrix {
    /// Construct and canonicalize an element-by-species matrix.
    ///
    /// Input rows correspond to `elements`; columns correspond to `species`.
    /// Both axes are sorted by identifier and the values are permuted with
    /// them, so logical equality and identity do not depend on caller order.
    ///
    /// # Errors
    /// Returns a typed error for empty/duplicate axes, shape defects, or a
    /// negative elemental count.
    pub fn new(
        elements: Vec<ElementId>,
        species: Vec<SpeciesId>,
        rows: Vec<Vec<i128>>,
    ) -> Result<Self, ChemistryError> {
        require_nonempty(
            ArtifactKind::ElementalMatrix,
            AxisKind::Elements,
            elements.len(),
        )?;
        require_nonempty(
            ArtifactKind::ElementalMatrix,
            AxisKind::Species,
            species.len(),
        )?;
        validate_matrix_shape(
            ArtifactKind::ElementalMatrix,
            elements.len(),
            species.len(),
            &rows,
        )?;
        let element_order =
            canonical_order(&elements).map_err(|index| ChemistryError::DuplicateId {
                kind: IdKind::Element,
                value: elements[index].as_str().to_string(),
            })?;
        let species_order =
            canonical_order(&species).map_err(|index| ChemistryError::DuplicateId {
                kind: IdKind::Species,
                value: species[index].as_str().to_string(),
            })?;

        let canonical_elements = reorder_owned(elements, &element_order);
        let canonical_species = reorder_owned(species, &species_order);
        let canonical_rows = reorder_owned(rows, &element_order);
        let mut values = Vec::new();
        for (row, element) in canonical_rows.into_iter().zip(&canonical_elements) {
            for (count, species_id) in reorder_owned(row, &species_order)
                .into_iter()
                .zip(&canonical_species)
            {
                if count < 0 {
                    return Err(ChemistryError::NegativeElementCount {
                        element: element.clone(),
                        species: species_id.clone(),
                        count,
                    });
                }
                values.push(count);
            }
        }
        Ok(Self {
            elements: canonical_elements,
            species: canonical_species,
            values,
        })
    }

    /// Canonically ordered element rows.
    #[must_use]
    pub fn elements(&self) -> &[ElementId] {
        &self.elements
    }

    /// Canonically ordered species columns.
    #[must_use]
    pub fn species(&self) -> &[SpeciesId] {
        &self.species
    }

    /// Number of element rows.
    #[must_use]
    pub fn row_count(&self) -> usize {
        self.elements.len()
    }

    /// Number of species columns.
    #[must_use]
    pub fn column_count(&self) -> usize {
        self.species.len()
    }

    /// Matrix entry by canonical row/column index.
    #[must_use]
    pub fn get(&self, row: usize, column: usize) -> Option<i128> {
        flat_index(row, column, self.column_count())
            .and_then(|index| self.values.get(index).copied())
    }

    /// Matrix entry by typed ids.
    #[must_use]
    pub fn count(&self, element: &ElementId, species: &SpeciesId) -> Option<i128> {
        let row = self.elements.binary_search(element).ok()?;
        let column = self.species.binary_search(species).ok()?;
        self.get(row, column)
    }

    /// Deterministic BLAKE3 identity over canonical axes and exact entries.
    #[must_use]
    pub fn identity(&self) -> ContentHash {
        let mut bytes = Vec::new();
        push_usize(&mut bytes, self.elements.len());
        for element in &self.elements {
            push_text(&mut bytes, element.as_str());
        }
        push_species(&mut bytes, &self.species);
        push_usize(&mut bytes, self.values.len());
        for value in &self.values {
            bytes.extend_from_slice(&value.to_le_bytes());
        }
        hash_domain(ELEMENTAL_MATRIX_ID_DOMAIN, &bytes)
    }
}

/// Immutable canonical species-by-reaction stoichiometric matrix `N`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StoichiometricMatrix {
    species: Vec<SpeciesId>,
    reactions: Vec<ReactionId>,
    values: Vec<i128>,
}

impl StoichiometricMatrix {
    /// Construct and canonicalize a species-by-reaction matrix.
    ///
    /// Reactants use negative coefficients and products use positive
    /// coefficients. Both axes are sorted with the entries permuted.
    ///
    /// # Errors
    /// Returns a typed error for empty/duplicate axes or shape defects.
    pub fn new(
        species: Vec<SpeciesId>,
        reactions: Vec<ReactionId>,
        rows: Vec<Vec<i128>>,
    ) -> Result<Self, ChemistryError> {
        require_nonempty(
            ArtifactKind::StoichiometricMatrix,
            AxisKind::Species,
            species.len(),
        )?;
        require_nonempty(
            ArtifactKind::StoichiometricMatrix,
            AxisKind::Reactions,
            reactions.len(),
        )?;
        validate_matrix_shape(
            ArtifactKind::StoichiometricMatrix,
            species.len(),
            reactions.len(),
            &rows,
        )?;
        let species_order =
            canonical_order(&species).map_err(|index| ChemistryError::DuplicateId {
                kind: IdKind::Species,
                value: species[index].as_str().to_string(),
            })?;
        let reaction_order =
            canonical_order(&reactions).map_err(|index| ChemistryError::DuplicateId {
                kind: IdKind::Reaction,
                value: reactions[index].as_str().to_string(),
            })?;

        let canonical_species = reorder_owned(species, &species_order);
        let canonical_reactions = reorder_owned(reactions, &reaction_order);
        let canonical_rows = reorder_owned(rows, &species_order);
        let mut values = Vec::new();
        for row in canonical_rows {
            values.extend(reorder_owned(row, &reaction_order));
        }
        Ok(Self {
            species: canonical_species,
            reactions: canonical_reactions,
            values,
        })
    }

    /// Canonically ordered species rows.
    #[must_use]
    pub fn species(&self) -> &[SpeciesId] {
        &self.species
    }

    /// Canonically ordered reaction columns.
    #[must_use]
    pub fn reactions(&self) -> &[ReactionId] {
        &self.reactions
    }

    /// Number of species rows.
    #[must_use]
    pub fn row_count(&self) -> usize {
        self.species.len()
    }

    /// Number of reaction columns.
    #[must_use]
    pub fn column_count(&self) -> usize {
        self.reactions.len()
    }

    /// Matrix entry by canonical row/column index.
    #[must_use]
    pub fn get(&self, row: usize, column: usize) -> Option<i128> {
        flat_index(row, column, self.column_count())
            .and_then(|index| self.values.get(index).copied())
    }

    /// Deterministic BLAKE3 identity over canonical axes and exact entries.
    #[must_use]
    pub fn identity(&self) -> ContentHash {
        let mut bytes = Vec::new();
        push_species(&mut bytes, &self.species);
        push_usize(&mut bytes, self.reactions.len());
        for reaction in &self.reactions {
            push_text(&mut bytes, reaction.as_str());
        }
        push_usize(&mut bytes, self.values.len());
        for value in &self.values {
            bytes.extend_from_slice(&value.to_le_bytes());
        }
        hash_domain(STOICHIOMETRIC_MATRIX_ID_DOMAIN, &bytes)
    }
}

/// Immutable canonical species charge vector `z` in elementary-charge units.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ChargeVector {
    species: Vec<SpeciesId>,
    values: Vec<i128>,
}

impl ChargeVector {
    /// Construct and canonicalize a species charge vector.
    ///
    /// # Errors
    /// Returns a typed error for an empty/duplicate species axis or a length
    /// mismatch.
    pub fn new(species: Vec<SpeciesId>, values: Vec<i128>) -> Result<Self, ChemistryError> {
        require_nonempty(ArtifactKind::ChargeVector, AxisKind::Species, species.len())?;
        if values.len() != species.len() {
            return Err(ChemistryError::VectorLengthMismatch {
                artifact: ArtifactKind::ChargeVector,
                expected: species.len(),
                actual: values.len(),
            });
        }
        let species_order =
            canonical_order(&species).map_err(|index| ChemistryError::DuplicateId {
                kind: IdKind::Species,
                value: species[index].as_str().to_string(),
            })?;
        let canonical_species = reorder_owned(species, &species_order);
        let canonical_values = reorder_owned(values, &species_order);
        Ok(Self {
            species: canonical_species,
            values: canonical_values,
        })
    }

    /// Canonically ordered species axis.
    #[must_use]
    pub fn species(&self) -> &[SpeciesId] {
        &self.species
    }

    /// Charge by canonical species index.
    #[must_use]
    pub fn get(&self, index: usize) -> Option<i128> {
        self.values.get(index).copied()
    }

    /// Charge by typed species id.
    #[must_use]
    pub fn charge(&self, species: &SpeciesId) -> Option<i128> {
        let index = self.species.binary_search(species).ok()?;
        self.get(index)
    }

    /// Deterministic BLAKE3 identity over canonical species and charges.
    #[must_use]
    pub fn identity(&self) -> ContentHash {
        let mut bytes = Vec::new();
        push_species(&mut bytes, &self.species);
        push_usize(&mut bytes, self.values.len());
        for value in &self.values {
            bytes.extend_from_slice(&value.to_le_bytes());
        }
        hash_domain(CHARGE_VECTOR_ID_DOMAIN, &bytes)
    }
}

fn require_nonempty(
    artifact: ArtifactKind,
    axis: AxisKind,
    length: usize,
) -> Result<(), ChemistryError> {
    if length == 0 {
        Err(ChemistryError::EmptyAxis { artifact, axis })
    } else {
        Ok(())
    }
}

fn validate_matrix_shape(
    artifact: ArtifactKind,
    expected_rows: usize,
    expected_columns: usize,
    rows: &[Vec<i128>],
) -> Result<(), ChemistryError> {
    if rows.len() != expected_rows {
        return Err(ChemistryError::RowCountMismatch {
            artifact,
            expected: expected_rows,
            actual: rows.len(),
        });
    }
    for (row, values) in rows.iter().enumerate() {
        if values.len() != expected_columns {
            return Err(ChemistryError::ColumnCountMismatch {
                artifact,
                row,
                expected: expected_columns,
                actual: values.len(),
            });
        }
    }
    Ok(())
}

fn flat_index(row: usize, column: usize, columns: usize) -> Option<usize> {
    if column >= columns {
        return None;
    }
    row.checked_mul(columns)?.checked_add(column)
}

/// Immutable exact G0 conservation certificate.
///
/// The private constructor is reachable only after all entries of `A * N`
/// and `z^T * N` have been checked equal to zero.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ConservationCertificate {
    elemental_matrix: ContentHash,
    stoichiometric_matrix: ContentHash,
    charge_vector: ContentHash,
    identity: ContentHash,
}

impl ConservationCertificate {
    /// Elemental-matrix artifact identity proved by this certificate.
    #[must_use]
    pub const fn elemental_matrix(&self) -> ContentHash {
        self.elemental_matrix
    }

    /// Stoichiometric-matrix artifact identity proved by this certificate.
    #[must_use]
    pub const fn stoichiometric_matrix(&self) -> ContentHash {
        self.stoichiometric_matrix
    }

    /// Charge-vector artifact identity proved by this certificate.
    #[must_use]
    pub const fn charge_vector(&self) -> ContentHash {
        self.charge_vector
    }

    /// Domain-separated identity binding all three proved artifacts.
    #[must_use]
    pub const fn identity(&self) -> ContentHash {
        self.identity
    }

    /// Whether these immutable artifacts are exactly those named by the proof.
    #[must_use]
    pub fn matches(
        &self,
        elemental: &ElementalMatrix,
        stoichiometric: &StoichiometricMatrix,
        charge: &ChargeVector,
    ) -> bool {
        self.elemental_matrix == elemental.identity()
            && self.stoichiometric_matrix == stoichiometric.identity()
            && self.charge_vector == charge.identity()
    }
}

/// Prove exact elemental and charge conservation.
///
/// Checks proceed in canonical element/reaction/species order, making the
/// first returned failure deterministic. Every multiply and add is checked;
/// overflow refuses rather than wrapping or weakening the proof.
///
/// # Errors
/// Returns a typed species-basis mismatch, checked-arithmetic overflow, or
/// exact imbalance naming the element/reaction (or charge/reaction).
pub fn verify_conservation(
    elemental: &ElementalMatrix,
    stoichiometric: &StoichiometricMatrix,
    charge: &ChargeVector,
) -> Result<ConservationCertificate, ChemistryError> {
    if elemental.species != stoichiometric.species {
        return Err(ChemistryError::SpeciesBasisMismatch {
            artifact: ArtifactKind::StoichiometricMatrix,
            expected: elemental.species.clone(),
            actual: stoichiometric.species.clone(),
        });
    }
    if elemental.species != charge.species {
        return Err(ChemistryError::SpeciesBasisMismatch {
            artifact: ArtifactKind::ChargeVector,
            expected: elemental.species.clone(),
            actual: charge.species.clone(),
        });
    }

    verify_element_balance(elemental, stoichiometric)?;
    verify_charge_balance(elemental, stoichiometric, charge)?;

    let elemental_matrix = elemental.identity();
    let stoichiometric_matrix = stoichiometric.identity();
    let charge_vector = charge.identity();
    let mut material = Vec::with_capacity(96);
    material.extend_from_slice(elemental_matrix.as_bytes());
    material.extend_from_slice(stoichiometric_matrix.as_bytes());
    material.extend_from_slice(charge_vector.as_bytes());
    Ok(ConservationCertificate {
        elemental_matrix,
        stoichiometric_matrix,
        charge_vector,
        identity: hash_domain(CONSERVATION_CERTIFICATE_ID_DOMAIN, &material),
    })
}

fn verify_element_balance(
    elemental: &ElementalMatrix,
    stoichiometric: &StoichiometricMatrix,
) -> Result<(), ChemistryError> {
    for (element_index, element) in elemental.elements.iter().enumerate() {
        for (reaction_index, reaction) in stoichiometric.reactions.iter().enumerate() {
            let mut residual = 0i128;
            for (species_index, species) in elemental.species.iter().enumerate() {
                let count =
                    elemental.values[element_index * elemental.column_count() + species_index];
                let coefficient = stoichiometric.values
                    [species_index * stoichiometric.column_count() + reaction_index];
                let term = count.checked_mul(coefficient).ok_or_else(|| {
                    arithmetic_overflow(ConservationLaw::Elements, Some(element), reaction, species)
                })?;
                residual = residual.checked_add(term).ok_or_else(|| {
                    arithmetic_overflow(ConservationLaw::Elements, Some(element), reaction, species)
                })?;
            }
            if residual != 0 {
                return Err(ChemistryError::ElementImbalance {
                    element: element.clone(),
                    reaction: reaction.clone(),
                    residual,
                });
            }
        }
    }
    Ok(())
}

fn verify_charge_balance(
    elemental: &ElementalMatrix,
    stoichiometric: &StoichiometricMatrix,
    charge: &ChargeVector,
) -> Result<(), ChemistryError> {
    for (reaction_index, reaction) in stoichiometric.reactions.iter().enumerate() {
        let mut residual = 0i128;
        for (species_index, species) in elemental.species.iter().enumerate() {
            let coefficient = stoichiometric.values
                [species_index * stoichiometric.column_count() + reaction_index];
            let term = charge.values[species_index]
                .checked_mul(coefficient)
                .ok_or_else(|| {
                    arithmetic_overflow(ConservationLaw::Charge, None, reaction, species)
                })?;
            residual = residual.checked_add(term).ok_or_else(|| {
                arithmetic_overflow(ConservationLaw::Charge, None, reaction, species)
            })?;
        }
        if residual != 0 {
            return Err(ChemistryError::ChargeImbalance {
                reaction: reaction.clone(),
                residual,
            });
        }
    }
    Ok(())
}

fn arithmetic_overflow(
    law: ConservationLaw,
    element: Option<&ElementId>,
    reaction: &ReactionId,
    species: &SpeciesId,
) -> ChemistryError {
    ChemistryError::ArithmeticOverflow {
        law,
        element: element.cloned(),
        reaction: reaction.clone(),
        species: species.clone(),
    }
}

/// Recorded typed relationship among mass, molar mass, and amount.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MassAmountBasis {
    source: BasisSource,
    mass: Mass,
    molar_mass: MolarMass,
    amount: Amount,
}

impl MassAmountBasis {
    /// Record a mass basis and derive amount as `mass / molar_mass`.
    ///
    /// # Errors
    /// Refuses non-finite/negative mass, non-positive molar mass, or a
    /// non-finite derived amount.
    pub fn from_mass(mass: Mass, molar_mass: MolarMass) -> Result<Self, ChemistryError> {
        validate_nonnegative(BasisField::Mass, mass.value())?;
        validate_positive(BasisField::MolarMass, molar_mass.value())?;
        let amount = Amount::new(mass.value() / molar_mass.value());
        validate_nonnegative(BasisField::DerivedAmount, amount.value())?;
        Ok(Self {
            source: BasisSource::Mass,
            mass,
            molar_mass,
            amount,
        })
    }

    /// Record an amount basis and derive mass as `amount * molar_mass`.
    ///
    /// # Errors
    /// Refuses non-finite/negative amount, non-positive molar mass, or a
    /// non-finite derived mass.
    pub fn from_amount(amount: Amount, molar_mass: MolarMass) -> Result<Self, ChemistryError> {
        validate_nonnegative(BasisField::Amount, amount.value())?;
        validate_positive(BasisField::MolarMass, molar_mass.value())?;
        let mass = Mass::new(amount.value() * molar_mass.value());
        validate_nonnegative(BasisField::DerivedMass, mass.value())?;
        Ok(Self {
            source: BasisSource::Amount,
            mass,
            molar_mass,
            amount,
        })
    }

    /// Quantity from which the other basis value was derived.
    #[must_use]
    pub const fn source(&self) -> BasisSource {
        self.source
    }

    /// Recorded mass in coherent SI units.
    #[must_use]
    pub const fn mass(&self) -> Mass {
        self.mass
    }

    /// Recorded molar mass in coherent SI units.
    #[must_use]
    pub const fn molar_mass(&self) -> MolarMass {
        self.molar_mass
    }

    /// Recorded amount in coherent SI units.
    #[must_use]
    pub const fn amount(&self) -> Amount {
        self.amount
    }
}

fn validate_nonnegative(field: BasisField, value: f64) -> Result<(), ChemistryError> {
    if !value.is_finite() || value < 0.0 {
        Err(ChemistryError::InvalidBasisValue {
            field,
            value,
            reason: "expected a finite nonnegative coherent-SI value",
        })
    } else {
        Ok(())
    }
}

fn validate_positive(field: BasisField, value: f64) -> Result<(), ChemistryError> {
    if !value.is_finite() || value <= 0.0 {
        Err(ChemistryError::InvalidBasisValue {
            field,
            value,
            reason: "expected a finite positive coherent-SI value",
        })
    } else {
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn species(value: &str) -> SpeciesId {
        SpeciesId::new(value).expect("valid species")
    }

    fn element(value: &str) -> ElementId {
        ElementId::new(value).expect("valid element")
    }

    fn reaction(value: &str) -> ReactionId {
        ReactionId::new(value).expect("valid reaction")
    }

    fn assert_close(actual: f64, expected: f64) {
        let tolerance = 8.0 * f64::EPSILON * expected.abs().max(1.0);
        assert!(
            (actual - expected).abs() <= tolerance,
            "expected {expected}, got {actual} (tolerance {tolerance})"
        );
    }

    #[test]
    fn g0_balanced_neutral_reaction_is_canonical_and_certified() {
        // 2 H2 + O2 -> 2 H2O, deliberately supplied in noncanonical order.
        let elemental = ElementalMatrix::new(
            vec![element("O"), element("H")],
            vec![species("O2"), species("H2O"), species("H2")],
            vec![vec![2, 1, 0], vec![0, 2, 2]],
        )
        .expect("elemental matrix");
        let stoichiometric = StoichiometricMatrix::new(
            vec![species("H2O"), species("O2"), species("H2")],
            vec![reaction("water")],
            vec![vec![2], vec![-1], vec![-2]],
        )
        .expect("stoichiometry");
        let charge = ChargeVector::new(
            vec![species("O2"), species("H2"), species("H2O")],
            vec![0, 0, 0],
        )
        .expect("charge");

        let certificate = verify_conservation(&elemental, &stoichiometric, &charge)
            .expect("neutral reaction balances exactly");
        assert!(certificate.matches(&elemental, &stoichiometric, &charge));
        assert_eq!(
            elemental
                .elements()
                .iter()
                .map(ElementId::as_str)
                .collect::<Vec<_>>(),
            vec!["H", "O"]
        );
        assert_eq!(
            elemental
                .species()
                .iter()
                .map(SpeciesId::as_str)
                .collect::<Vec<_>>(),
            vec!["H2", "H2O", "O2"]
        );

        let same_elemental = ElementalMatrix::new(
            vec![element("H"), element("O")],
            vec![species("H2"), species("H2O"), species("O2")],
            vec![vec![2, 2, 0], vec![0, 1, 2]],
        )
        .expect("canonical matrix");
        assert_eq!(elemental, same_elemental);
        assert_eq!(elemental.identity(), same_elemental.identity());
    }

    #[test]
    fn g0_balanced_charged_reaction_is_certified() {
        // Fe2+ -> Fe3+ + e-
        let elemental = ElementalMatrix::new(
            vec![element("Fe")],
            vec![species("e-"), species("Fe3+"), species("Fe2+")],
            vec![vec![0, 1, 1]],
        )
        .expect("elemental matrix");
        let stoichiometric = StoichiometricMatrix::new(
            vec![species("Fe2+"), species("e-"), species("Fe3+")],
            vec![reaction("oxidation")],
            vec![vec![-1], vec![1], vec![1]],
        )
        .expect("stoichiometry");
        let charge = ChargeVector::new(
            vec![species("Fe3+"), species("Fe2+"), species("e-")],
            vec![3, 2, -1],
        )
        .expect("charge");
        verify_conservation(&elemental, &stoichiometric, &charge)
            .expect("charged reaction balances exactly");
    }

    #[test]
    fn g0_element_imbalance_names_element_and_reaction() {
        let elemental = ElementalMatrix::new(
            vec![element("H"), element("O")],
            vec![species("H2"), species("H2O")],
            vec![vec![2, 2], vec![0, 1]],
        )
        .unwrap();
        let stoichiometric = StoichiometricMatrix::new(
            vec![species("H2"), species("H2O")],
            vec![reaction("hydrate")],
            vec![vec![-1], vec![1]],
        )
        .unwrap();
        let charge = ChargeVector::new(vec![species("H2"), species("H2O")], vec![0, 0]).unwrap();
        assert!(matches!(
            verify_conservation(&elemental, &stoichiometric, &charge),
            Err(ChemistryError::ElementImbalance {
                element,
                reaction,
                residual: 1,
            }) if element.as_str() == "O" && reaction.as_str() == "hydrate"
        ));
    }

    #[test]
    fn g0_charge_imbalance_names_reaction() {
        let elemental = ElementalMatrix::new(
            vec![element("Fe")],
            vec![species("Fe2+"), species("Fe3+")],
            vec![vec![1, 1]],
        )
        .unwrap();
        let stoichiometric = StoichiometricMatrix::new(
            vec![species("Fe2+"), species("Fe3+")],
            vec![reaction("bad-oxidation")],
            vec![vec![-1], vec![1]],
        )
        .unwrap();
        let charge = ChargeVector::new(vec![species("Fe2+"), species("Fe3+")], vec![2, 3]).unwrap();
        assert!(matches!(
            verify_conservation(&elemental, &stoichiometric, &charge),
            Err(ChemistryError::ChargeImbalance {
                reaction,
                residual: 1,
            }) if reaction.as_str() == "bad-oxidation"
        ));
    }

    #[test]
    fn malformed_and_duplicate_ids_are_refused() {
        assert!(matches!(
            SpeciesId::new(""),
            Err(ChemistryError::MalformedId {
                kind: IdKind::Species,
                ..
            })
        ));
        assert!(SpeciesId::new("H 2").is_err());
        assert!(ElementId::new("fe").is_err());
        assert!(ReactionId::new("water/reverse").is_err());

        let duplicate = ElementalMatrix::new(
            vec![element("H")],
            vec![species("H2"), species("H2")],
            vec![vec![2, 2]],
        );
        assert!(matches!(
            duplicate,
            Err(ChemistryError::DuplicateId {
                kind: IdKind::Species,
                value,
            }) if value == "H2"
        ));
    }

    #[test]
    fn malformed_shapes_and_negative_counts_are_refused() {
        assert!(matches!(
            ElementalMatrix::new(
                vec![element("H"), element("O")],
                vec![species("H2")],
                vec![vec![2]],
            ),
            Err(ChemistryError::RowCountMismatch {
                artifact: ArtifactKind::ElementalMatrix,
                expected: 2,
                actual: 1,
            })
        ));
        assert!(matches!(
            StoichiometricMatrix::new(
                vec![species("H2"), species("O2")],
                vec![reaction("r")],
                vec![vec![-1], vec![]],
            ),
            Err(ChemistryError::ColumnCountMismatch {
                artifact: ArtifactKind::StoichiometricMatrix,
                row: 1,
                expected: 1,
                actual: 0,
            })
        ));
        assert!(matches!(
            ChargeVector::new(vec![species("H2"), species("O2")], vec![0]),
            Err(ChemistryError::VectorLengthMismatch {
                artifact: ArtifactKind::ChargeVector,
                expected: 2,
                actual: 1,
            })
        ));
        assert!(matches!(
            ElementalMatrix::new(vec![element("H")], vec![species("H2")], vec![vec![-2]],),
            Err(ChemistryError::NegativeElementCount { count: -2, .. })
        ));
    }

    #[test]
    fn matrix_accessors_reject_both_out_of_bounds_axes() {
        let elemental = ElementalMatrix::new(
            vec![element("H"), element("O")],
            vec![species("H2"), species("O2")],
            vec![vec![2, 0], vec![0, 2]],
        )
        .unwrap();
        assert_eq!(elemental.get(0, elemental.column_count()), None);
        assert_eq!(elemental.get(elemental.row_count(), 0), None);

        let stoichiometric = StoichiometricMatrix::new(
            vec![species("H2"), species("O2")],
            vec![reaction("forward"), reaction("reverse")],
            vec![vec![-2, 2], vec![-1, 1]],
        )
        .unwrap();
        assert_eq!(stoichiometric.get(0, stoichiometric.column_count()), None);
        assert_eq!(stoichiometric.get(stoichiometric.row_count(), 0), None);
    }

    #[test]
    fn empty_and_duplicate_named_axes_are_refused() {
        assert!(matches!(
            ElementalMatrix::new(Vec::new(), vec![species("H2")], Vec::new()),
            Err(ChemistryError::EmptyAxis {
                artifact: ArtifactKind::ElementalMatrix,
                axis: AxisKind::Elements,
            })
        ));
        assert!(matches!(
            StoichiometricMatrix::new(vec![species("H2")], Vec::new(), vec![Vec::new()]),
            Err(ChemistryError::EmptyAxis {
                artifact: ArtifactKind::StoichiometricMatrix,
                axis: AxisKind::Reactions,
            })
        ));
        assert!(matches!(
            ChargeVector::new(Vec::new(), Vec::new()),
            Err(ChemistryError::EmptyAxis {
                artifact: ArtifactKind::ChargeVector,
                axis: AxisKind::Species,
            })
        ));
        assert!(matches!(
            ElementalMatrix::new(
                vec![element("H"), element("H")],
                vec![species("H2")],
                vec![vec![2], vec![2]],
            ),
            Err(ChemistryError::DuplicateId {
                kind: IdKind::Element,
                value,
            }) if value == "H"
        ));
        assert!(matches!(
            StoichiometricMatrix::new(
                vec![species("H2")],
                vec![reaction("r"), reaction("r")],
                vec![vec![-1, 1]],
            ),
            Err(ChemistryError::DuplicateId {
                kind: IdKind::Reaction,
                value,
            }) if value == "r"
        ));
    }

    #[test]
    fn species_basis_mismatch_names_the_inconsistent_artifact() {
        let elemental =
            ElementalMatrix::new(vec![element("H")], vec![species("H2")], vec![vec![2]]).unwrap();
        let stoichiometric =
            StoichiometricMatrix::new(vec![species("O2")], vec![reaction("r")], vec![vec![0]])
                .unwrap();
        let charge = ChargeVector::new(vec![species("H2")], vec![0]).unwrap();
        assert!(matches!(
            verify_conservation(&elemental, &stoichiometric, &charge),
            Err(ChemistryError::SpeciesBasisMismatch {
                artifact: ArtifactKind::StoichiometricMatrix,
                ..
            })
        ));
    }

    #[test]
    fn stoichiometric_and_charge_identities_ignore_input_axis_order() {
        let first_n = StoichiometricMatrix::new(
            vec![species("B"), species("A")],
            vec![reaction("reverse"), reaction("forward")],
            vec![vec![-2, 2], vec![1, -1]],
        )
        .unwrap();
        let second_n = StoichiometricMatrix::new(
            vec![species("A"), species("B")],
            vec![reaction("forward"), reaction("reverse")],
            vec![vec![-1, 1], vec![2, -2]],
        )
        .unwrap();
        assert_eq!(first_n, second_n);
        assert_eq!(first_n.identity(), second_n.identity());

        let first_z = ChargeVector::new(vec![species("B"), species("A")], vec![1, -2]).unwrap();
        let second_z = ChargeVector::new(vec![species("A"), species("B")], vec![-2, 1]).unwrap();
        assert_eq!(first_z, second_z);
        assert_eq!(first_z.identity(), second_z.identity());
    }

    #[test]
    fn checked_i128_overflow_refuses_the_proof() {
        let elemental = ElementalMatrix::new(
            vec![element("X")],
            vec![species("X")],
            vec![vec![i128::MAX]],
        )
        .unwrap();
        let stoichiometric = StoichiometricMatrix::new(
            vec![species("X")],
            vec![reaction("overflow")],
            vec![vec![2]],
        )
        .unwrap();
        let charge = ChargeVector::new(vec![species("X")], vec![0]).unwrap();
        assert!(matches!(
            verify_conservation(&elemental, &stoichiometric, &charge),
            Err(ChemistryError::ArithmeticOverflow {
                law: ConservationLaw::Elements,
                reaction,
                species,
                ..
            }) if reaction.as_str() == "overflow" && species.as_str() == "X"
        ));
    }

    #[test]
    fn checked_i128_accumulation_overflow_refuses_both_laws() {
        let elemental = ElementalMatrix::new(
            vec![element("X")],
            vec![species("A"), species("B")],
            vec![vec![i128::MAX, 1]],
        )
        .unwrap();
        let stoichiometric = StoichiometricMatrix::new(
            vec![species("A"), species("B")],
            vec![reaction("element-overflow")],
            vec![vec![1], vec![1]],
        )
        .unwrap();
        let charge = ChargeVector::new(vec![species("A"), species("B")], vec![0, 0]).unwrap();
        assert!(matches!(
            verify_conservation(&elemental, &stoichiometric, &charge),
            Err(ChemistryError::ArithmeticOverflow {
                law: ConservationLaw::Elements,
                species,
                ..
            }) if species.as_str() == "B"
        ));

        let zero_elements = ElementalMatrix::new(
            vec![element("X")],
            vec![species("A"), species("B")],
            vec![vec![0, 0]],
        )
        .unwrap();
        let charge =
            ChargeVector::new(vec![species("A"), species("B")], vec![i128::MAX, 1]).unwrap();
        assert!(matches!(
            verify_conservation(&zero_elements, &stoichiometric, &charge),
            Err(ChemistryError::ArithmeticOverflow {
                law: ConservationLaw::Charge,
                species,
                ..
            }) if species.as_str() == "B"
        ));
    }

    #[test]
    fn typed_mass_amount_basis_records_its_source() {
        let molar_mass = MolarMass::new(0.018);
        let from_mass = MassAmountBasis::from_mass(Mass::new(18.0), molar_mass).unwrap();
        assert_eq!(from_mass.source(), BasisSource::Mass);
        assert_close(from_mass.amount().value(), 1000.0);

        let from_amount = MassAmountBasis::from_amount(Amount::new(1000.0), molar_mass).unwrap();
        assert_eq!(from_amount.source(), BasisSource::Amount);
        assert_close(from_amount.mass().value(), 18.0);
        assert!(MassAmountBasis::from_mass(Mass::new(-1.0), molar_mass).is_err());
        assert!(MassAmountBasis::from_amount(Amount::new(1.0), MolarMass::new(0.0)).is_err());
        assert!(MassAmountBasis::from_mass(Mass::new(f64::NAN), molar_mass).is_err());
        assert!(MassAmountBasis::from_amount(Amount::new(f64::INFINITY), molar_mass).is_err());
        assert!(MassAmountBasis::from_mass(Mass::new(1.0), MolarMass::new(f64::NAN)).is_err());
    }
}
