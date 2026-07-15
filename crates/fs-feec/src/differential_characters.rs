//! Typed schema for finite relative differential characters.
//!
//! This module fixes the object, coefficient, degree, relative-boundary, and
//! map conventions consumed by the constructive algebra in later RA.2
//! slices.  It deliberately does **not** claim that a finite cochain model
//! approximates smooth differential cohomology.  Exact-sequence entries are
//! therefore explicit schemas whose kernel/image obligations remain pending
//! until an admitted constructive checker supplies witnesses.
//!
//! The degree convention is Cheeger--Simons/Hopkins--Singer: an object in
//! `DiffChar^k` has a degree-`k` curvature and characteristic class, and its
//! holonomy pairs with `(k - 1)`-cycles.  Cellular coboundary is written
//! `delta`; de Rham differentiation is written `d`.  They are never exposed
//! through one ambiguous operator.

use core::fmt;
use core::num::NonZeroU16;
use fs_qty::Dims;
use std::collections::BTreeSet;

/// Canonical schema version for the RA.2a object algebra.
pub const DIFFERENTIAL_CHARACTER_SCHEMA_VERSION: u32 = 1;

/// A deterministic 256-bit content address for an algebra schema.
///
/// The digest is a stable, length-framed four-lane FNV construction.  It is
/// suitable for deterministic replay and accidental-corruption detection, but
/// it is **not** a cryptographic authenticity proof.  Ledger authority must
/// wrap it in the workspace's admitted identity/receipt machinery.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct AlgebraId([u8; 32]);

impl AlgebraId {
    /// Construct an identity from an already content-addressed immutable
    /// artifact.  Authentication remains the caller's ledger concern.
    #[must_use]
    pub const fn from_bytes(bytes: [u8; 32]) -> Self {
        Self(bytes)
    }

    /// Raw digest bytes.
    #[must_use]
    pub const fn as_bytes(&self) -> &[u8; 32] {
        &self.0
    }

    /// Lower-case hexadecimal representation.
    #[must_use]
    pub fn to_hex(self) -> String {
        const HEX: &[u8; 16] = b"0123456789abcdef";
        let mut out = String::with_capacity(64);
        for byte in self.0 {
            out.push(char::from(HEX[usize::from(byte >> 4)]));
            out.push(char::from(HEX[usize::from(byte & 0x0f)]));
        }
        out
    }
}

impl fmt::Display for AlgebraId {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(&self.to_hex())
    }
}

/// Primal and dual complexes are separate nominal lanes.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ComplexLane {
    /// The oriented primal cell complex.
    Primal,
    /// The oriented dual cell complex.
    Dual,
}

/// Chosen orientation of a finite complex.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum Orientation {
    /// Reference orientation.
    Positive,
    /// Reversed reference orientation.
    Negative,
}

/// Cohomological degree `k`.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct CohomologicalDegree(u8);

impl CohomologicalDegree {
    /// Construct a degree.  Compatibility with a complex dimension is checked
    /// when the degree is attached to an object.
    #[must_use]
    pub const fn new(value: u8) -> Self {
        Self(value)
    }

    /// Numeric degree.
    #[must_use]
    pub const fn get(self) -> u8 {
        self.0
    }

    /// Degree immediately below this one.
    #[must_use]
    pub const fn predecessor(self) -> Option<Self> {
        match self.0.checked_sub(1) {
            Some(value) => Some(Self(value)),
            None => None,
        }
    }

    /// Degree immediately above this one.
    #[must_use]
    pub const fn successor(self) -> Option<Self> {
        match self.0.checked_add(1) {
            Some(value) => Some(Self(value)),
            None => None,
        }
    }
}

/// Explicit finite-schema construction limits.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct AlgebraBudget {
    /// Maximum cells across the ambient and relative complexes.
    max_cells: u64,
    /// Maximum named boundary components.
    max_boundary_components: u32,
    /// Maximum bytes emitted by canonical schema serialization.
    max_canonical_bytes: u64,
}

impl AlgebraBudget {
    /// Checked budget constructor.
    pub fn new(
        max_cells: u64,
        max_boundary_components: u32,
        max_canonical_bytes: u64,
    ) -> Result<Self, CharacterError> {
        if max_cells == 0 || max_boundary_components == 0 || max_canonical_bytes == 0 {
            return Err(CharacterError::ZeroBudget);
        }
        Ok(Self {
            max_cells,
            max_boundary_components,
            max_canonical_bytes,
        })
    }

    /// Cell-count limit.
    #[must_use]
    pub const fn max_cells(self) -> u64 {
        self.max_cells
    }

    /// Boundary-component limit.
    #[must_use]
    pub const fn max_boundary_components(self) -> u32 {
        self.max_boundary_components
    }

    /// Canonical-byte limit.
    #[must_use]
    pub const fn max_canonical_bytes(self) -> u64 {
        self.max_canonical_bytes
    }
}

/// Finite oriented complex metadata used by the schema layer.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FiniteComplexSchema {
    id: String,
    revision: u32,
    dimension: u8,
    lane: ComplexLane,
    orientation: Orientation,
    cell_counts: Vec<u64>,
}

impl FiniteComplexSchema {
    /// Construct a finite complex schema with one cell count for every degree
    /// `0..=dimension`.
    pub fn new(
        id: impl Into<String>,
        revision: u32,
        dimension: u8,
        lane: ComplexLane,
        orientation: Orientation,
        cell_counts: Vec<u64>,
    ) -> Result<Self, CharacterError> {
        let id = nonempty(id.into(), "complex id")?;
        if revision == 0 {
            return Err(CharacterError::ZeroRevision { object: "complex" });
        }
        if cell_counts.len() != usize::from(dimension) + 1 {
            return Err(CharacterError::CellCountArity {
                dimension,
                actual: cell_counts.len(),
            });
        }
        Ok(Self {
            id,
            revision,
            dimension,
            lane,
            orientation,
            cell_counts,
        })
    }

    /// Stable caller-supplied identity label.
    #[must_use]
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Semantic revision of this finite complex.
    #[must_use]
    pub const fn revision(&self) -> u32 {
        self.revision
    }

    /// Top cell degree.
    #[must_use]
    pub const fn dimension(&self) -> u8 {
        self.dimension
    }

    /// Primal or dual lane.
    #[must_use]
    pub const fn lane(&self) -> ComplexLane {
        self.lane
    }

    /// Chosen orientation.
    #[must_use]
    pub const fn orientation(&self) -> Orientation {
        self.orientation
    }

    /// Cell counts by degree.
    #[must_use]
    pub fn cell_counts(&self) -> &[u64] {
        &self.cell_counts
    }

    /// Whether every represented chain group has rank zero.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.cell_counts.iter().all(|count| *count == 0)
    }

    fn total_cells(&self) -> Result<u64, CharacterError> {
        self.cell_counts.iter().try_fold(0_u64, |sum, count| {
            sum.checked_add(*count)
                .ok_or(CharacterError::CellCountOverflow)
        })
    }

    fn encode(&self, encoder: &mut CanonicalEncoder) {
        encoder.string(&self.id);
        encoder.u32(self.revision);
        encoder.u8(self.dimension);
        encoder.lane(self.lane);
        encoder.orientation(self.orientation);
        encoder.u64(usize_to_u64(self.cell_counts.len()));
        for count in &self.cell_counts {
            encoder.u64(*count);
        }
    }
}

/// Role of a named component of the relative subcomplex.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum BoundaryRole {
    /// Ordinary relative boundary on which the character is trivialized.
    Relative,
    /// Physical terminal that requires an explicit `(k - 1)` trivialization.
    Terminal,
    /// Material/interface component, retained separately from terminals.
    Interface,
}

/// Orientation convention for an embedded boundary component.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum BoundaryOrientation {
    /// Boundary orientation induced by the ambient orientation.
    Induced,
    /// Opposite of the induced boundary orientation.
    Opposite,
}

/// Finite model selected for relative differential cohomology.
///
/// RA.2a admits the mapping-cone model only: an object on `(X, A)` carries a
/// degree-`k - 1` trivialization of its restriction on the **entire**
/// subcomplex `A`.  The nominal tag is part of the pair identity so a future
/// inequivalent relative model cannot silently reuse artifacts.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum RelativeModel {
    /// Hopkins--Singer mapping-cone relative object.
    MappingCone,
}

/// Named boundary component metadata.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BoundaryComponent {
    id: String,
    dimension: u8,
    lane: ComplexLane,
    orientation: BoundaryOrientation,
    role: BoundaryRole,
}

impl BoundaryComponent {
    /// Construct a component.  Codimension and lane compatibility are checked
    /// by [`RelativePairSchema::new`].
    pub fn new(
        id: impl Into<String>,
        dimension: u8,
        lane: ComplexLane,
        orientation: BoundaryOrientation,
        role: BoundaryRole,
    ) -> Result<Self, CharacterError> {
        Ok(Self {
            id: nonempty(id.into(), "boundary id")?,
            dimension,
            lane,
            orientation,
            role,
        })
    }

    /// Component identity.
    #[must_use]
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Cell dimension.
    #[must_use]
    pub const fn dimension(&self) -> u8 {
        self.dimension
    }

    /// Primal or dual lane.
    #[must_use]
    pub const fn lane(&self) -> ComplexLane {
        self.lane
    }

    /// Boundary orientation convention.
    #[must_use]
    pub const fn orientation(&self) -> BoundaryOrientation {
        self.orientation
    }

    /// Semantic boundary role.
    #[must_use]
    pub const fn role(&self) -> BoundaryRole {
        self.role
    }

    fn encode(&self, encoder: &mut CanonicalEncoder) {
        encoder.string(&self.id);
        encoder.u8(self.dimension);
        encoder.lane(self.lane);
        encoder.boundary_orientation(self.orientation);
        encoder.boundary_role(self.role);
    }
}

/// A finite relative pair `(X, A)` with explicit boundary decomposition.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RelativePairSchema {
    id: String,
    revision: u32,
    ambient: FiniteComplexSchema,
    relative: FiniteComplexSchema,
    relative_model: RelativeModel,
    boundaries: Vec<BoundaryComponent>,
    budget: AlgebraBudget,
    algebra_id: AlgebraId,
}

impl RelativePairSchema {
    /// Construct and validate `(X, A)`.  `A` must use the same primal/dual lane
    /// and have no more cells than `X` in any represented degree.  Ambient
    /// dimension 255 is refused because this schema must represent the
    /// Hopkins--Singer ceiling `dim(X) + 1` in [`CohomologicalDegree`].
    pub fn new(
        id: impl Into<String>,
        revision: u32,
        ambient: FiniteComplexSchema,
        relative: FiniteComplexSchema,
        relative_model: RelativeModel,
        mut boundaries: Vec<BoundaryComponent>,
        budget: AlgebraBudget,
    ) -> Result<Self, CharacterError> {
        let id = nonempty(id.into(), "relative-pair id")?;
        if revision == 0 {
            return Err(CharacterError::ZeroRevision {
                object: "relative pair",
            });
        }
        if ambient.dimension == u8::MAX {
            return Err(CharacterError::CharacterDegreeRepresentationOverflow {
                dimension: ambient.dimension,
            });
        }
        if ambient.lane != relative.lane {
            return Err(CharacterError::LaneMismatch {
                expected: ambient.lane,
                actual: relative.lane,
                object: "relative subcomplex",
            });
        }
        if relative.dimension > ambient.dimension {
            return Err(CharacterError::RelativeDimension {
                ambient: ambient.dimension,
                relative: relative.dimension,
            });
        }
        for (degree, relative_count) in relative.cell_counts.iter().enumerate() {
            let ambient_count = ambient.cell_counts.get(degree).copied().unwrap_or(0);
            if *relative_count > ambient_count {
                return Err(CharacterError::RelativeCellCount {
                    degree: usize_to_u8(degree),
                    ambient: ambient_count,
                    relative: *relative_count,
                });
            }
        }
        let boundary_count =
            u32::try_from(boundaries.len()).map_err(|_| CharacterError::BoundaryBudgetExceeded)?;
        if boundary_count > budget.max_boundary_components {
            return Err(CharacterError::BoundaryBudgetExceeded);
        }
        let total_cells = ambient
            .total_cells()?
            .checked_add(relative.total_cells()?)
            .ok_or(CharacterError::CellCountOverflow)?;
        if total_cells > budget.max_cells {
            return Err(CharacterError::CellBudgetExceeded {
                requested: total_cells,
                limit: budget.max_cells,
            });
        }

        let mut ids = BTreeSet::new();
        for boundary in &boundaries {
            if boundary.lane != ambient.lane {
                return Err(CharacterError::LaneMismatch {
                    expected: ambient.lane,
                    actual: boundary.lane,
                    object: "boundary component",
                });
            }
            if boundary.dimension > relative.dimension {
                return Err(CharacterError::BoundaryDimension {
                    boundary: boundary.dimension,
                    relative: relative.dimension,
                });
            }
            if boundary.role == BoundaryRole::Terminal
                && boundary.dimension.checked_add(1) != Some(ambient.dimension)
            {
                return Err(CharacterError::TerminalCodimension {
                    ambient: ambient.dimension,
                    terminal: boundary.dimension,
                });
            }
            if !ids.insert(boundary.id.clone()) {
                return Err(CharacterError::DuplicateBoundary {
                    id: boundary.id.clone(),
                });
            }
        }
        boundaries.sort_by(|left, right| left.id.cmp(&right.id));

        let mut pair = Self {
            id,
            revision,
            ambient,
            relative,
            relative_model,
            boundaries,
            budget,
            algebra_id: AlgebraId([0; 32]),
        };
        let canonical = pair.canonical_bytes_without_id()?;
        pair.algebra_id = stable_digest(&canonical);
        Ok(pair)
    }

    /// Human-readable pair identity.
    #[must_use]
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Pair semantic revision.
    #[must_use]
    pub const fn revision(&self) -> u32 {
        self.revision
    }

    /// Ambient complex `X`.
    #[must_use]
    pub const fn ambient(&self) -> &FiniteComplexSchema {
        &self.ambient
    }

    /// Relative subcomplex `A`.
    #[must_use]
    pub const fn relative(&self) -> &FiniteComplexSchema {
        &self.relative
    }

    /// Selected finite relative-cohomology model.
    #[must_use]
    pub const fn relative_model(&self) -> RelativeModel {
        self.relative_model
    }

    /// Named decomposition of the relative boundary.
    #[must_use]
    pub fn boundaries(&self) -> &[BoundaryComponent] {
        &self.boundaries
    }

    /// Construction budget committed by this pair.
    #[must_use]
    pub const fn budget(&self) -> AlgebraBudget {
        self.budget
    }

    /// Deterministic content identity.
    #[must_use]
    pub const fn algebra_id(&self) -> AlgebraId {
        self.algebra_id
    }

    /// Locate a named boundary component.
    #[must_use]
    pub fn boundary(&self, id: &str) -> Option<&BoundaryComponent> {
        self.boundaries.iter().find(|boundary| boundary.id == id)
    }

    /// Canonical, versioned bytes committed by [`Self::algebra_id`].
    #[must_use]
    pub fn canonical_bytes(&self) -> Result<Vec<u8>, CharacterError> {
        self.canonical_bytes_without_id()
    }

    fn canonical_bytes_without_id(&self) -> Result<Vec<u8>, CharacterError> {
        let mut encoder = CanonicalEncoder::new(
            b"frankensim.fs-feec.relative-pair",
            self.budget.max_canonical_bytes,
        );
        encoder.u32(DIFFERENTIAL_CHARACTER_SCHEMA_VERSION);
        encoder.string(&self.id);
        encoder.u32(self.revision);
        self.ambient.encode(&mut encoder);
        self.relative.encode(&mut encoder);
        encoder.relative_model(self.relative_model);
        encoder.u64(usize_to_u64(self.boundaries.len()));
        for boundary in &self.boundaries {
            boundary.encode(&mut encoder);
        }
        encoder.budget(self.budget);
        encoder.finish()
    }
}

/// Versioned lattice embedded in a real coefficient space.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CoefficientLattice {
    id: String,
    revision: u32,
    units: Dims,
    rank: NonZeroU16,
    generator_scale_bits: Vec<u64>,
}

impl CoefficientLattice {
    /// Construct a diagonalized, named lattice.  Generator scales are explicit
    /// positive finite SI-normalization factors; a general basis change is a
    /// later constructive-algebra concern and must not be smuggled into these
    /// scalars.
    pub fn new(
        id: impl Into<String>,
        revision: u32,
        units: Dims,
        generator_scales: &[f64],
    ) -> Result<Self, CharacterError> {
        let id = nonempty(id.into(), "coefficient-lattice id")?;
        if revision == 0 {
            return Err(CharacterError::ZeroRevision {
                object: "coefficient lattice",
            });
        }
        if generator_scales.is_empty() || generator_scales.len() > usize::from(u16::MAX) {
            return Err(CharacterError::InvalidLatticeRank {
                rank: generator_scales.len(),
            });
        }
        let rank = u16::try_from(generator_scales.len())
            .ok()
            .and_then(NonZeroU16::new)
            .ok_or(CharacterError::InvalidLatticeRank {
                rank: generator_scales.len(),
            })?;
        let mut bits = Vec::with_capacity(generator_scales.len());
        for (index, scale) in generator_scales.iter().copied().enumerate() {
            if !scale.is_finite() || scale <= 0.0 {
                return Err(CharacterError::InvalidLatticeScale { index });
            }
            bits.push(scale.to_bits());
        }
        Ok(Self {
            id,
            revision,
            units,
            rank,
            generator_scale_bits: bits,
        })
    }

    /// Canonical rank-one dimensionless integer coefficients for cellular
    /// chains/cycles.  Physical flux/charge lattices belong to character
    /// values, not to the cycles on which holonomy is evaluated.
    #[must_use]
    pub fn dimensionless_integer_cycles() -> Self {
        Self {
            id: "frankensim.integral-cycles.z.v1".to_owned(),
            revision: 1,
            units: Dims::NONE,
            rank: NonZeroU16::MIN,
            generator_scale_bits: vec![1.0_f64.to_bits()],
        }
    }

    /// Lattice identity.
    #[must_use]
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Semantic revision.
    #[must_use]
    pub const fn revision(&self) -> u32 {
        self.revision
    }

    /// Coefficient rank.
    #[must_use]
    pub const fn rank(&self) -> NonZeroU16 {
        self.rank
    }

    /// Physical dimensions of every normalized generator.
    #[must_use]
    pub const fn units(&self) -> Dims {
        self.units
    }

    /// Normalization scale for one generator.
    #[must_use]
    pub fn generator_scale(&self, index: usize) -> Option<f64> {
        self.generator_scale_bits
            .get(index)
            .copied()
            .map(f64::from_bits)
    }

    fn encode(&self, encoder: &mut CanonicalEncoder) {
        encoder.string(&self.id);
        encoder.u32(self.revision);
        encoder.dims(self.units);
        encoder.u16(self.rank.get());
        encoder.u64(usize_to_u64(self.generator_scale_bits.len()));
        for bits in &self.generator_scale_bits {
            encoder.u64(*bits);
        }
    }
}

/// Valid modulus for a finite cyclic coefficient group.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct CyclicModulus(u64);

impl CyclicModulus {
    /// Construct `Z/n`; `n < 2` is not a torsion coefficient group.
    pub fn new(modulus: u64) -> Result<Self, CharacterError> {
        if modulus < 2 {
            return Err(CharacterError::InvalidCyclicModulus { modulus });
        }
        Ok(Self(modulus))
    }

    /// Numeric modulus.
    #[must_use]
    pub const fn get(self) -> u64 {
        self.0
    }
}

/// Mutually exclusive coefficient sectors.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum CoefficientSector {
    /// Integral lattice sector.
    Integral,
    /// Finite torsion sector.
    Torsion,
    /// Real vector-space sector.
    Real,
    /// Real vector space modulo a named integral lattice.
    RealModuloLattice,
    /// Real coefficients constrained to have periods in a named lattice.
    RealWithLatticePeriods,
}

/// Coefficient group attached to a cochain, class, or character.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CoefficientSystem {
    /// Integral coefficients in a named lattice `Lambda`.
    Integral(CoefficientLattice),
    /// Finite cyclic coefficients `Z/n`.
    Torsion(CyclicModulus),
    /// Real coefficients of the stated rank and units.
    Real {
        /// Nonzero vector-space rank.
        rank: NonZeroU16,
        /// Physical dimensions of the coefficient normalization.
        units: Dims,
    },
    /// Circle/torus-valued coefficients `R^r / Lambda`.
    RealModuloLattice(CoefficientLattice),
    /// Real coefficients whose closed cohomology class lies in the image of a
    /// named lattice.  The lattice is retained so curvature codomains with
    /// equal rank/units but different normalization remain nominally distinct.
    RealWithLatticePeriods(CoefficientLattice),
}

impl CoefficientSystem {
    /// Coarse sector, deliberately preserving integral/torsion/real nominal
    /// distinctions.
    #[must_use]
    pub const fn sector(&self) -> CoefficientSector {
        match self {
            Self::Integral(_) => CoefficientSector::Integral,
            Self::Torsion(_) => CoefficientSector::Torsion,
            Self::Real { .. } => CoefficientSector::Real,
            Self::RealModuloLattice(_) => CoefficientSector::RealModuloLattice,
            Self::RealWithLatticePeriods(_) => CoefficientSector::RealWithLatticePeriods,
        }
    }

    /// Realification of a lattice-valued coefficient system.
    pub fn realification(&self) -> Result<Self, CharacterError> {
        match self {
            Self::Integral(lattice)
            | Self::RealModuloLattice(lattice)
            | Self::RealWithLatticePeriods(lattice) => Ok(Self::Real {
                rank: lattice.rank(),
                units: lattice.units(),
            }),
            Self::Real { .. } => Ok(self.clone()),
            Self::Torsion(_) => Err(CharacterError::TorsionHasNoCurvature),
        }
    }

    /// Integral/torsion coefficient system targeted by a characteristic map.
    pub fn characteristic_coefficients(&self) -> Result<Self, CharacterError> {
        match self {
            Self::Integral(lattice) | Self::RealModuloLattice(lattice) => {
                Ok(Self::Integral(lattice.clone()))
            }
            Self::Torsion(modulus) => Ok(Self::Torsion(*modulus)),
            Self::Real { .. } => Err(CharacterError::RealHasNoCharacteristicClass),
            Self::RealWithLatticePeriods(_) => {
                Err(CharacterError::CurvatureSpaceHasNoCharacteristicClass)
            }
        }
    }

    fn encode(&self, encoder: &mut CanonicalEncoder) {
        match self {
            Self::Integral(lattice) => {
                encoder.u8(0);
                lattice.encode(encoder);
            }
            Self::Torsion(modulus) => {
                encoder.u8(1);
                encoder.u64(modulus.get());
            }
            Self::Real { rank, units } => {
                encoder.u8(2);
                encoder.u16(rank.get());
                encoder.dims(*units);
            }
            Self::RealModuloLattice(lattice) => {
                encoder.u8(3);
                lattice.encode(encoder);
            }
            Self::RealWithLatticePeriods(lattice) => {
                encoder.u8(4);
                lattice.encode(encoder);
            }
        }
    }
}

/// Concrete representative family.  Variants are nominally coupled to their
/// coefficient sector; cross-sector coercions are rejected.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum RepresentativeKind {
    /// Discrete Hopkins--Singer triple `(c, h, omega)`.
    HopkinsSingerTriple,
    /// Integral cocycle without real quotient data.
    IntegralCocycle,
    /// Flat finite-torsion cocycle.
    FlatTorsionCocycle,
    /// Ordinary real cochain (not a differential character claim).
    RealCochain,
}

/// A relative or terminal boundary trivialization lives one degree below its
/// character.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BoundaryTrivialization {
    boundary_id: String,
    degree: CohomologicalDegree,
    lane: ComplexLane,
    coefficients: CoefficientSystem,
    representative_id: String,
}

impl BoundaryTrivialization {
    /// Construct boundary trivialization metadata.  The relative-pair
    /// constructor decides whether the named component admits one.
    pub fn new(
        boundary_id: impl Into<String>,
        degree: CohomologicalDegree,
        lane: ComplexLane,
        coefficients: CoefficientSystem,
        representative_id: impl Into<String>,
    ) -> Result<Self, CharacterError> {
        Ok(Self {
            boundary_id: nonempty(boundary_id.into(), "boundary id")?,
            degree,
            lane,
            coefficients,
            representative_id: nonempty(representative_id.into(), "boundary representative id")?,
        })
    }

    /// Named relative or terminal boundary.
    #[must_use]
    pub fn boundary_id(&self) -> &str {
        &self.boundary_id
    }

    /// Trivialization cochain degree.
    #[must_use]
    pub const fn degree(&self) -> CohomologicalDegree {
        self.degree
    }

    /// Primal or dual lane.
    #[must_use]
    pub const fn lane(&self) -> ComplexLane {
        self.lane
    }

    /// Coefficient group.
    #[must_use]
    pub const fn coefficients(&self) -> &CoefficientSystem {
        &self.coefficients
    }

    /// Content/provenance reference for the actual representative.
    #[must_use]
    pub fn representative_id(&self) -> &str {
        &self.representative_id
    }

    fn encode(&self, encoder: &mut CanonicalEncoder) {
        encoder.string(&self.boundary_id);
        encoder.u8(self.degree.get());
        encoder.lane(self.lane);
        self.coefficients.encode(encoder);
        encoder.string(&self.representative_id);
    }
}

/// Terminal-specialized spelling retained for call sites that are declaring a
/// physical terminal rather than an ordinary relative boundary.
pub type TerminalTrivialization = BoundaryTrivialization;

/// Mapping-cone trivialization of the character restriction on the entire
/// relative subcomplex `A`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RelativeTrivialization {
    degree: CohomologicalDegree,
    lane: ComplexLane,
    coefficients: CoefficientSystem,
    representative_id: String,
}

impl RelativeTrivialization {
    /// Construct whole-subcomplex trivialization metadata.
    pub fn new(
        degree: CohomologicalDegree,
        lane: ComplexLane,
        coefficients: CoefficientSystem,
        representative_id: impl Into<String>,
    ) -> Result<Self, CharacterError> {
        Ok(Self {
            degree,
            lane,
            coefficients,
            representative_id: nonempty(
                representative_id.into(),
                "relative-subcomplex representative id",
            )?,
        })
    }

    /// Degree `k - 1` of the mapping-cone trivialization.
    #[must_use]
    pub const fn degree(&self) -> CohomologicalDegree {
        self.degree
    }

    /// Primal/dual lane.
    #[must_use]
    pub const fn lane(&self) -> ComplexLane {
        self.lane
    }

    /// Coefficient group.
    #[must_use]
    pub const fn coefficients(&self) -> &CoefficientSystem {
        &self.coefficients
    }

    /// Immutable representative/provenance reference.
    #[must_use]
    pub fn representative_id(&self) -> &str {
        &self.representative_id
    }

    fn encode(&self, encoder: &mut CanonicalEncoder) {
        encoder.u8(self.degree.get());
        encoder.lane(self.lane);
        self.coefficients.encode(encoder);
        encoder.string(&self.representative_id);
    }
}

/// Validated differential-character object schema.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct RelativeDifferentialCharacter {
    pair: RelativePairSchema,
    degree: CohomologicalDegree,
    coefficients: CoefficientSystem,
    representative: RepresentativeKind,
    relative_trivialization: Option<RelativeTrivialization>,
    boundary_trivializations: Vec<BoundaryTrivialization>,
    algebra_id: AlgebraId,
}

impl RelativeDifferentialCharacter {
    /// Construct a relative object and enforce degree, coefficient,
    /// primal/dual, and terminal-boundary compatibility.
    pub fn new(
        pair: RelativePairSchema,
        degree: CohomologicalDegree,
        coefficients: CoefficientSystem,
        representative: RepresentativeKind,
        relative_trivialization: Option<RelativeTrivialization>,
        mut boundary_trivializations: Vec<BoundaryTrivialization>,
    ) -> Result<Self, CharacterError> {
        let maximum_degree = match representative {
            RepresentativeKind::HopkinsSingerTriple => pair.ambient.dimension + 1,
            RepresentativeKind::IntegralCocycle
            | RepresentativeKind::FlatTorsionCocycle
            | RepresentativeKind::RealCochain => pair.ambient.dimension,
        };
        if degree.get() > maximum_degree {
            return Err(CharacterError::DegreeOutOfRange {
                degree: degree.get(),
                maximum: maximum_degree,
            });
        }
        let expected_sector = match representative {
            RepresentativeKind::HopkinsSingerTriple => CoefficientSector::RealModuloLattice,
            RepresentativeKind::IntegralCocycle => CoefficientSector::Integral,
            RepresentativeKind::FlatTorsionCocycle => CoefficientSector::Torsion,
            RepresentativeKind::RealCochain => CoefficientSector::Real,
        };
        if coefficients.sector() != expected_sector {
            return Err(CharacterError::RepresentativeCoefficientMismatch {
                representative,
                expected: expected_sector,
                actual: coefficients.sector(),
            });
        }

        let requires_relative_trivialization = !pair.relative.is_empty();
        let requires_component_trivialization = pair
            .boundaries
            .iter()
            .any(|boundary| boundary.role == BoundaryRole::Terminal);
        let trivialization_degree = if !requires_relative_trivialization
            && !requires_component_trivialization
            && relative_trivialization.is_none()
            && boundary_trivializations.is_empty()
        {
            None
        } else {
            Some(
                degree
                    .predecessor()
                    .ok_or(CharacterError::DegreeZeroHasNoBoundaryTrivialization)?,
            )
        };

        match (
            requires_relative_trivialization,
            relative_trivialization.as_ref(),
        ) {
            (true, None) => return Err(CharacterError::MissingRelativeSubcomplexTrivialization),
            (false, Some(_)) => {
                return Err(CharacterError::UnexpectedRelativeSubcomplexTrivialization);
            }
            (false, None) => {}
            (true, Some(trivialization)) => {
                let expected_degree = trivialization_degree
                    .ok_or(CharacterError::DegreeZeroHasNoBoundaryTrivialization)?;
                if trivialization.degree != expected_degree {
                    return Err(CharacterError::BoundaryTrivializationDegreeMismatch {
                        expected: expected_degree.get(),
                        actual: trivialization.degree.get(),
                    });
                }
                if trivialization.lane != pair.ambient.lane {
                    return Err(CharacterError::LaneMismatch {
                        expected: pair.ambient.lane,
                        actual: trivialization.lane,
                        object: "relative-subcomplex trivialization",
                    });
                }
                if trivialization.coefficients != coefficients {
                    return Err(CharacterError::CoefficientMismatch {
                        object: "relative-subcomplex trivialization",
                    });
                }
            }
        }

        let mut presented = BTreeSet::new();
        for trivialization in &boundary_trivializations {
            let Some(boundary) = pair.boundary(&trivialization.boundary_id) else {
                return Err(CharacterError::UnknownBoundary {
                    id: trivialization.boundary_id.clone(),
                });
            };
            if boundary.role == BoundaryRole::Interface {
                return Err(CharacterError::BoundaryDoesNotAdmitTrivialization {
                    id: boundary.id.clone(),
                });
            }
            let Some(expected_trivialization_degree) = trivialization_degree else {
                return Err(CharacterError::DegreeZeroHasNoBoundaryTrivialization);
            };
            if trivialization.degree != expected_trivialization_degree {
                return Err(CharacterError::BoundaryTrivializationDegreeMismatch {
                    expected: expected_trivialization_degree.get(),
                    actual: trivialization.degree.get(),
                });
            }
            if trivialization.lane != pair.ambient.lane {
                return Err(CharacterError::LaneMismatch {
                    expected: pair.ambient.lane,
                    actual: trivialization.lane,
                    object: "boundary trivialization",
                });
            }
            if trivialization.coefficients != coefficients {
                return Err(CharacterError::CoefficientMismatch {
                    object: "boundary trivialization",
                });
            }
            if !presented.insert(trivialization.boundary_id.as_str()) {
                return Err(CharacterError::DuplicateBoundaryTrivialization {
                    id: trivialization.boundary_id.clone(),
                });
            }
        }
        for boundary in pair
            .boundaries
            .iter()
            .filter(|boundary| boundary.role == BoundaryRole::Terminal)
        {
            if !presented.contains(boundary.id.as_str()) {
                return Err(CharacterError::MissingTerminalTrivialization {
                    id: boundary.id.clone(),
                });
            }
        }
        boundary_trivializations.sort_by(|left, right| {
            left.boundary_id
                .cmp(&right.boundary_id)
                .then_with(|| left.representative_id.cmp(&right.representative_id))
        });

        let mut object = Self {
            pair,
            degree,
            coefficients,
            representative,
            relative_trivialization,
            boundary_trivializations,
            algebra_id: AlgebraId([0; 32]),
        };
        let canonical = object.canonical_bytes_without_id()?;
        object.algebra_id = stable_digest(&canonical);
        Ok(object)
    }

    /// Relative pair `(X, A)`.
    #[must_use]
    pub const fn pair(&self) -> &RelativePairSchema {
        &self.pair
    }

    /// Character degree.
    #[must_use]
    pub const fn degree(&self) -> CohomologicalDegree {
        self.degree
    }

    /// Coefficient system.
    #[must_use]
    pub const fn coefficients(&self) -> &CoefficientSystem {
        &self.coefficients
    }

    /// Representative family.
    #[must_use]
    pub const fn representative(&self) -> RepresentativeKind {
        self.representative
    }

    /// Mapping-cone trivialization on the entire relative subcomplex `A`.
    #[must_use]
    pub const fn relative_trivialization(&self) -> Option<&RelativeTrivialization> {
        self.relative_trivialization.as_ref()
    }

    /// Relative and terminal trivializations in canonical boundary-identity
    /// order.
    #[must_use]
    pub fn boundary_trivializations(&self) -> &[BoundaryTrivialization] {
        &self.boundary_trivializations
    }

    /// Deterministic schema identity.
    #[must_use]
    pub const fn algebra_id(&self) -> AlgebraId {
        self.algebra_id
    }

    /// Canonical, versioned schema bytes.
    #[must_use]
    pub fn canonical_bytes(&self) -> Result<Vec<u8>, CharacterError> {
        self.canonical_bytes_without_id()
    }

    /// Typed curvature map `curv: DiffChar^k -> Z^k_R`.
    pub fn curvature_map(&self) -> Result<TypedMap, CharacterError> {
        self.ensure_differential_character()?;
        let CoefficientSystem::RealModuloLattice(lattice) = &self.coefficients else {
            return Err(CharacterError::OperationRequiresDifferentialCharacter);
        };
        TypedMap::new(
            MapKind::Curvature,
            self.object(
                ObjectKind::DifferentialCharacters,
                self.degree,
                self.coefficients.clone(),
            ),
            self.object(
                ObjectKind::ClosedIntegralCurvatures,
                self.degree,
                CoefficientSystem::RealWithLatticePeriods(lattice.clone()),
            ),
        )
    }

    /// Typed characteristic map `c: DiffChar^k -> H^k(-; Lambda)`.
    pub fn characteristic_class_map(&self) -> Result<TypedMap, CharacterError> {
        self.ensure_differential_character()?;
        TypedMap::new(
            MapKind::CharacteristicClass,
            self.object(
                ObjectKind::DifferentialCharacters,
                self.degree,
                self.coefficients.clone(),
            ),
            self.object(
                ObjectKind::CohomologyClasses,
                self.degree,
                self.coefficients.characteristic_coefficients()?,
            ),
        )
    }

    /// Gauge-equivalence schema.  It states the quotient and degree of gauge
    /// representatives; it does not pretend the downstream orbit checker ran.
    pub fn gauge_equivalence(&self) -> Result<GaugeEquivalenceSchema, CharacterError> {
        self.ensure_differential_character()?;
        let gauge_degree = self
            .degree
            .predecessor()
            .ok_or(CharacterError::DegreeZeroHasNoGaugeParameter)?;
        Ok(GaugeEquivalenceSchema {
            gauge_parameters: self.object(
                ObjectKind::GaugeParameters,
                gauge_degree,
                self.coefficients.clone(),
            ),
            representatives: self.object(
                ObjectKind::CharacterRepresentatives,
                self.degree,
                self.coefficients.clone(),
            ),
            quotient: self.object(
                ObjectKind::DifferentialCharacters,
                self.degree,
                self.coefficients.clone(),
            ),
            cellular_nilpotence: NilpotenceLaw::CellularCoboundarySquaredZero,
            de_rham_nilpotence: NilpotenceLaw::DeRhamDifferentialSquaredZero,
        })
    }

    /// Standard curvature short exact-sequence schema
    /// `0 -> H^(k-1)(R/Lambda) -> DiffChar^k -> Z^k_Lambda(R) -> 0`.
    pub fn curvature_exact_sequence(&self) -> Result<ExactSequenceSchema, CharacterError> {
        self.ensure_differential_character()?;
        let CoefficientSystem::RealModuloLattice(lattice) = &self.coefficients else {
            return Err(CharacterError::StandardSequenceRequiresRealModuloLattice);
        };
        let previous = self
            .degree
            .predecessor()
            .ok_or(CharacterError::DegreeZeroHasNoExactSequence)?;
        let quotient = CoefficientSystem::RealModuloLattice(lattice.clone());
        let real_with_periods = CoefficientSystem::RealWithLatticePeriods(lattice.clone());
        let objects = vec![
            self.object(ObjectKind::Zero, previous, quotient.clone()),
            self.object(ObjectKind::FlatCohomology, previous, quotient),
            self.object(
                ObjectKind::DifferentialCharacters,
                self.degree,
                self.coefficients.clone(),
            ),
            self.object(
                ObjectKind::ClosedIntegralCurvatures,
                self.degree,
                real_with_periods.clone(),
            ),
            self.object(ObjectKind::Zero, self.degree, real_with_periods),
        ];
        ExactSequenceSchema::short_exact(
            ExactSequenceKind::Curvature,
            objects,
            [
                MapKind::Zero,
                MapKind::FlatInclusion,
                MapKind::Curvature,
                MapKind::Zero,
            ],
        )
    }

    /// Standard characteristic-class short exact-sequence schema
    /// `0 -> C^(k-1)(R)/C^(k-1)_Lambda -> DiffChar^k -> H^k(Lambda) -> 0`.
    pub fn characteristic_exact_sequence(&self) -> Result<ExactSequenceSchema, CharacterError> {
        self.ensure_differential_character()?;
        let CoefficientSystem::RealModuloLattice(lattice) = &self.coefficients else {
            return Err(CharacterError::StandardSequenceRequiresRealModuloLattice);
        };
        let previous = self
            .degree
            .predecessor()
            .ok_or(CharacterError::DegreeZeroHasNoExactSequence)?;
        let quotient = CoefficientSystem::RealModuloLattice(lattice.clone());
        let integral = CoefficientSystem::Integral(lattice.clone());
        let objects = vec![
            self.object(ObjectKind::Zero, previous, quotient.clone()),
            self.object(
                ObjectKind::RealCochainsModuloIntegralCocycles,
                previous,
                quotient,
            ),
            self.object(
                ObjectKind::DifferentialCharacters,
                self.degree,
                self.coefficients.clone(),
            ),
            self.object(ObjectKind::CohomologyClasses, self.degree, integral.clone()),
            self.object(ObjectKind::Zero, self.degree, integral),
        ];
        ExactSequenceSchema::short_exact(
            ExactSequenceKind::CharacteristicClass,
            objects,
            [
                MapKind::Zero,
                MapKind::TopologicallyTrivialInclusion,
                MapKind::CharacteristicClass,
                MapKind::Zero,
            ],
        )
    }

    /// Relative-pair long exact-sequence window
    /// `DiffChar^(k-1)(A) -> DiffChar^k(X,A) -> DiffChar^k(X) ->
    /// DiffChar^k(A) -> DiffChar^(k+1)(X,A)`.
    pub fn boundary_exact_sequence(&self) -> Result<ExactSequenceSchema, CharacterError> {
        self.ensure_differential_character()?;
        let previous = self
            .degree
            .predecessor()
            .ok_or(CharacterError::DegreeZeroHasNoExactSequence)?;
        let next = self
            .degree
            .successor()
            .filter(|degree| degree.get() <= self.maximum_character_degree())
            .ok_or(CharacterError::DegreeHasNoSuccessorInComplex {
                degree: self.degree.get(),
                maximum: self.maximum_character_degree(),
            })?;
        let coefficient = self.coefficients.clone();
        let objects = vec![
            self.object_with_support(
                ObjectKind::DifferentialCharacters,
                ObjectSupport::RelativeSubcomplex,
                previous,
                coefficient.clone(),
            ),
            self.object_with_support(
                ObjectKind::DifferentialCharacters,
                ObjectSupport::RelativePair,
                self.degree,
                coefficient.clone(),
            ),
            self.object_with_support(
                ObjectKind::DifferentialCharacters,
                ObjectSupport::AmbientComplex,
                self.degree,
                coefficient.clone(),
            ),
            self.object_with_support(
                ObjectKind::DifferentialCharacters,
                ObjectSupport::RelativeSubcomplex,
                self.degree,
                coefficient.clone(),
            ),
            self.object_with_support(
                ObjectKind::DifferentialCharacters,
                ObjectSupport::RelativePair,
                next,
                coefficient,
            ),
        ];
        ExactSequenceSchema::long_window(
            ExactSequenceKind::RelativeBoundary,
            objects,
            [
                MapKind::Connecting,
                MapKind::ForgetRelative,
                MapKind::BoundaryRestriction,
                MapKind::Connecting,
            ],
        )
    }

    /// Typed cup-product schema.  The output degree is `p + q`; coefficients,
    /// pair, support, and primal/dual lane must agree exactly.
    pub fn cup_product(
        &self,
        other: &RelativeDifferentialCharacter,
        coefficient_product: &CoefficientProductSchema,
    ) -> Result<BilinearMapSchema, CharacterError> {
        self.ensure_differential_character()?;
        other.ensure_differential_character()?;
        if self.pair.algebra_id != other.pair.algebra_id {
            return Err(CharacterError::PairMismatch);
        }
        if coefficient_product.left != self.coefficients {
            return Err(CharacterError::CoefficientProductMismatch { input: "left" });
        }
        if coefficient_product.right != other.coefficients {
            return Err(CharacterError::CoefficientProductMismatch { input: "right" });
        }
        if coefficient_product.output.sector() != CoefficientSector::RealModuloLattice {
            return Err(CharacterError::CoefficientProductMismatch { input: "output" });
        }
        let degree = self
            .degree
            .get()
            .checked_add(other.degree.get())
            .filter(|degree| *degree <= self.maximum_character_degree())
            .ok_or(CharacterError::ProductDegreeOutOfRange)?;
        Ok(BilinearMapSchema {
            kind: BilinearMapKind::CupProduct,
            left: self.object(
                ObjectKind::DifferentialCharacters,
                self.degree,
                self.coefficients.clone(),
            ),
            right: other.object(
                ObjectKind::DifferentialCharacters,
                other.degree,
                other.coefficients.clone(),
            ),
            output: self.object(
                ObjectKind::DifferentialCharacters,
                CohomologicalDegree::new(degree),
                coefficient_product.output.clone(),
            ),
            coefficient_rule: BilinearCoefficientRule::Declared(Box::new(
                coefficient_product.clone(),
            )),
        })
    }

    /// Typed holonomy pairing against relative `(k - 1)`-cycles.
    pub fn holonomy_pairing(&self) -> Result<BilinearMapSchema, CharacterError> {
        self.ensure_differential_character()?;
        let cycle_degree = self
            .degree
            .predecessor()
            .ok_or(CharacterError::DegreeZeroHasNoHolonomyCycle)?;
        let CoefficientSystem::RealModuloLattice(_) = &self.coefficients else {
            return Err(CharacterError::HolonomyRequiresRealModuloLattice);
        };
        Ok(BilinearMapSchema {
            kind: BilinearMapKind::HolonomyPairing,
            left: self.object(
                ObjectKind::DifferentialCharacters,
                self.degree,
                self.coefficients.clone(),
            ),
            right: self.object(
                ObjectKind::RelativeCycles,
                cycle_degree,
                CoefficientSystem::Integral(CoefficientLattice::dimensionless_integer_cycles()),
            ),
            output: self.object_with_support(
                ObjectKind::CoefficientValues,
                ObjectSupport::Point,
                CohomologicalDegree::new(0),
                self.coefficients.clone(),
            ),
            coefficient_rule: BilinearCoefficientRule::HolonomyEvaluation,
        })
    }

    /// Typed space of the declared representative family.  Bare integral,
    /// torsion, and real cochains remain supporting objects and are never
    /// mislabeled as differential characters.
    #[must_use]
    pub fn object_space(&self) -> ObjectSpace {
        let kind = match self.representative {
            RepresentativeKind::HopkinsSingerTriple => ObjectKind::DifferentialCharacters,
            RepresentativeKind::IntegralCocycle => ObjectKind::IntegralCocycles,
            RepresentativeKind::FlatTorsionCocycle => ObjectKind::TorsionCocycles,
            RepresentativeKind::RealCochain => ObjectKind::RealCochains,
        };
        self.object(kind, self.degree, self.coefficients.clone())
    }

    fn ensure_differential_character(&self) -> Result<(), CharacterError> {
        if self.representative == RepresentativeKind::HopkinsSingerTriple {
            Ok(())
        } else {
            Err(CharacterError::OperationRequiresDifferentialCharacter)
        }
    }

    fn maximum_character_degree(&self) -> u8 {
        self.pair.ambient.dimension + 1
    }

    fn object(
        &self,
        kind: ObjectKind,
        degree: CohomologicalDegree,
        coefficients: CoefficientSystem,
    ) -> ObjectSpace {
        self.object_with_support(kind, ObjectSupport::RelativePair, degree, coefficients)
    }

    fn object_with_support(
        &self,
        kind: ObjectKind,
        support: ObjectSupport,
        degree: CohomologicalDegree,
        coefficients: CoefficientSystem,
    ) -> ObjectSpace {
        ObjectSpace {
            pair: self.pair.algebra_id,
            lane: self.pair.ambient.lane,
            support,
            degree,
            coefficients,
            kind,
        }
    }

    fn canonical_bytes_without_id(&self) -> Result<Vec<u8>, CharacterError> {
        let mut encoder = CanonicalEncoder::new(
            b"frankensim.fs-feec.relative-diffchar",
            self.pair.budget.max_canonical_bytes,
        );
        encoder.u32(DIFFERENTIAL_CHARACTER_SCHEMA_VERSION);
        encoder.bytes(self.pair.algebra_id.as_bytes());
        encoder.u8(self.degree.get());
        self.coefficients.encode(&mut encoder);
        encoder.representative(self.representative);
        match &self.relative_trivialization {
            Some(trivialization) => {
                encoder.u8(1);
                trivialization.encode(&mut encoder);
            }
            None => encoder.u8(0),
        }
        encoder.u64(usize_to_u64(self.boundary_trivializations.len()));
        for trivialization in &self.boundary_trivializations {
            trivialization.encode(&mut encoder);
        }
        encoder.finish()
    }
}

/// Support carried by a typed algebra object.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ObjectSupport {
    /// Relative object on `(X, A)`.
    RelativePair,
    /// Absolute object on `X`.
    AmbientComplex,
    /// Absolute object on `A`.
    RelativeSubcomplex,
    /// Scalar output support.
    Point,
}

/// Semantic kind of a typed algebra object.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ObjectKind {
    /// Zero object with otherwise explicit type parameters.
    Zero,
    /// Relative or absolute differential-character group.
    DifferentialCharacters,
    /// Integral cocycles used as characteristic representatives.
    IntegralCocycles,
    /// Flat finite-torsion cocycles.
    TorsionCocycles,
    /// Ordinary real cochains.
    RealCochains,
    /// Hopkins--Singer or other representative space before quotienting.
    CharacterRepresentatives,
    /// Gauge cochains of degree `k - 1`.
    GaugeParameters,
    /// Closed real cochains with lattice periods.
    ClosedIntegralCurvatures,
    /// Flat `R/Lambda` cohomology.
    FlatCohomology,
    /// Integral or torsion cohomology classes.
    CohomologyClasses,
    /// Real `(k - 1)`-cochains modulo closed real cocycles whose class lies in
    /// the lattice image.  Nonclosed lattice-valued cochains are not quotiented
    /// away because they change curvature.
    RealCochainsModuloIntegralCocycles,
    /// Relative cycles paired with holonomy.
    RelativeCycles,
    /// Scalar/vector coefficient values.
    CoefficientValues,
}

/// Fully typed domain or codomain.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ObjectSpace {
    pair: AlgebraId,
    lane: ComplexLane,
    support: ObjectSupport,
    degree: CohomologicalDegree,
    coefficients: CoefficientSystem,
    kind: ObjectKind,
}

impl ObjectSpace {
    /// Relative-pair content identity.
    #[must_use]
    pub const fn pair_id(&self) -> AlgebraId {
        self.pair
    }

    /// Primal/dual lane.
    #[must_use]
    pub const fn lane(&self) -> ComplexLane {
        self.lane
    }

    /// Relative/absolute support.
    #[must_use]
    pub const fn support(&self) -> ObjectSupport {
        self.support
    }

    /// Cohomological degree.
    #[must_use]
    pub const fn degree(&self) -> CohomologicalDegree {
        self.degree
    }

    /// Coefficient group.
    #[must_use]
    pub const fn coefficients(&self) -> &CoefficientSystem {
        &self.coefficients
    }

    /// Semantic object kind.
    #[must_use]
    pub const fn kind(&self) -> ObjectKind {
        self.kind
    }
}

/// Named typed homomorphism.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum MapKind {
    /// Unique zero map.
    Zero,
    /// Inclusion of flat classes.
    FlatInclusion,
    /// Inclusion of topologically trivial representatives.
    TopologicallyTrivialInclusion,
    /// Curvature map.
    Curvature,
    /// Characteristic-class map.
    CharacteristicClass,
    /// Connecting homomorphism in the pair sequence.
    Connecting,
    /// Forget relative trivialization.
    ForgetRelative,
    /// Restrict an ambient character to `A`.
    BoundaryRestriction,
}

/// Typed map with explicit degree shift.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TypedMap {
    kind: MapKind,
    domain: ObjectSpace,
    codomain: ObjectSpace,
    degree_shift: i16,
}

impl TypedMap {
    fn new(
        kind: MapKind,
        domain: ObjectSpace,
        codomain: ObjectSpace,
    ) -> Result<Self, CharacterError> {
        if domain.pair != codomain.pair {
            return Err(CharacterError::PairMismatch);
        }
        if domain.lane != codomain.lane {
            return Err(CharacterError::LaneMismatch {
                expected: domain.lane,
                actual: codomain.lane,
                object: "typed map codomain",
            });
        }
        let degree_shift = i16::from(codomain.degree.get()) - i16::from(domain.degree.get());
        Ok(Self {
            kind,
            domain,
            codomain,
            degree_shift,
        })
    }

    /// Map semantics.
    #[must_use]
    pub const fn kind(&self) -> MapKind {
        self.kind
    }

    /// Typed domain.
    #[must_use]
    pub const fn domain(&self) -> &ObjectSpace {
        &self.domain
    }

    /// Typed codomain.
    #[must_use]
    pub const fn codomain(&self) -> &ObjectSpace {
        &self.codomain
    }

    /// Codomain degree minus domain degree.
    #[must_use]
    pub const fn degree_shift(&self) -> i16 {
        self.degree_shift
    }
}

/// The two nilpotence laws are nominally distinct.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum NilpotenceLaw {
    /// Cellular `delta ∘ delta = 0` on cochains.
    CellularCoboundarySquaredZero,
    /// De Rham `d ∘ d = 0` on forms/curvature representatives.
    DeRhamDifferentialSquaredZero,
}

/// Gauge quotient schema attached to one character family.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GaugeEquivalenceSchema {
    /// Degree-`k - 1` gauge parameters.
    gauge_parameters: ObjectSpace,
    /// Representative space before quotienting.
    representatives: ObjectSpace,
    /// Gauge-equivalence quotient.
    quotient: ObjectSpace,
    /// Cellular nilpotence obligation.
    cellular_nilpotence: NilpotenceLaw,
    /// de Rham nilpotence obligation.
    de_rham_nilpotence: NilpotenceLaw,
}

impl GaugeEquivalenceSchema {
    /// Degree-`k - 1` gauge parameters.
    #[must_use]
    pub const fn gauge_parameters(&self) -> &ObjectSpace {
        &self.gauge_parameters
    }

    /// Representative space before quotienting.
    #[must_use]
    pub const fn representatives(&self) -> &ObjectSpace {
        &self.representatives
    }

    /// Gauge-equivalence quotient.
    #[must_use]
    pub const fn quotient(&self) -> &ObjectSpace {
        &self.quotient
    }

    /// Cellular nilpotence obligation.
    #[must_use]
    pub const fn cellular_nilpotence(&self) -> NilpotenceLaw {
        self.cellular_nilpotence
    }

    /// de Rham nilpotence obligation.
    #[must_use]
    pub const fn de_rham_nilpotence(&self) -> NilpotenceLaw {
        self.de_rham_nilpotence
    }
}

/// Explicit coefficient-level bilinear map used by a cup product.
///
/// Equality of two coefficient groups does not manufacture a multiplication:
/// arbitrary lattices and physical units require a declared map
/// `left x right -> output`.  This schema names and versions that choice.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CoefficientProductSchema {
    id: String,
    revision: u32,
    left: CoefficientSystem,
    right: CoefficientSystem,
    output: CoefficientSystem,
    map_artifact: AlgebraId,
    budget: AlgebraBudget,
    algebra_id: AlgebraId,
}

impl CoefficientProductSchema {
    /// Construct a named coefficient product.  Algebraic associativity,
    /// graded commutativity, and normalization remain constructive-checker
    /// obligations; the schema never infers them from matching type names.
    pub fn new(
        id: impl Into<String>,
        revision: u32,
        left: CoefficientSystem,
        right: CoefficientSystem,
        output: CoefficientSystem,
        map_artifact: AlgebraId,
        budget: AlgebraBudget,
    ) -> Result<Self, CharacterError> {
        let id = nonempty(id.into(), "coefficient-product id")?;
        if revision == 0 {
            return Err(CharacterError::ZeroRevision {
                object: "coefficient product",
            });
        }
        let mut schema = Self {
            id,
            revision,
            left,
            right,
            output,
            map_artifact,
            budget,
            algebra_id: AlgebraId([0; 32]),
        };
        let canonical = schema.canonical_bytes_without_id()?;
        schema.algebra_id = stable_digest(&canonical);
        Ok(schema)
    }

    /// Human-readable product identity.
    #[must_use]
    pub fn id(&self) -> &str {
        &self.id
    }

    /// Semantic revision.
    #[must_use]
    pub const fn revision(&self) -> u32 {
        self.revision
    }

    /// Left coefficient group.
    #[must_use]
    pub const fn left(&self) -> &CoefficientSystem {
        &self.left
    }

    /// Right coefficient group.
    #[must_use]
    pub const fn right(&self) -> &CoefficientSystem {
        &self.right
    }

    /// Output coefficient group.
    #[must_use]
    pub const fn output(&self) -> &CoefficientSystem {
        &self.output
    }

    /// Immutable artifact containing the actual bilinear map/normalization.
    #[must_use]
    pub const fn map_artifact(&self) -> AlgebraId {
        self.map_artifact
    }

    /// Explicit construction budget.
    #[must_use]
    pub const fn budget(&self) -> AlgebraBudget {
        self.budget
    }

    /// Canonical, bounded coefficient-product bytes.
    pub fn canonical_bytes(&self) -> Result<Vec<u8>, CharacterError> {
        self.canonical_bytes_without_id()
    }

    /// Deterministic coefficient-map identity.
    #[must_use]
    pub const fn algebra_id(&self) -> AlgebraId {
        self.algebra_id
    }

    fn canonical_bytes_without_id(&self) -> Result<Vec<u8>, CharacterError> {
        let mut encoder = CanonicalEncoder::new(
            b"frankensim.fs-feec.coefficient-product",
            self.budget.max_canonical_bytes,
        );
        encoder.u32(DIFFERENTIAL_CHARACTER_SCHEMA_VERSION);
        encoder.string(&self.id);
        encoder.u32(self.revision);
        self.left.encode(&mut encoder);
        self.right.encode(&mut encoder);
        self.output.encode(&mut encoder);
        encoder.bytes(self.map_artifact.as_bytes());
        encoder.budget(self.budget);
        encoder.finish()
    }
}

/// Exact-sequence family.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ExactSequenceKind {
    /// Curvature short exact sequence.
    Curvature,
    /// Characteristic-class short exact sequence.
    CharacteristicClass,
    /// Relative-pair long exact sequence window.
    RelativeBoundary,
}

/// Finite assumptions under which a later checker may discharge an exactness
/// claim.  Merely constructing this enum never upgrades a claim to proven.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum FiniteExactnessAssumption {
    /// Boundary/coboundary matrices form a finite chain complex over the named
    /// coefficient system, with exact arithmetic and admitted quotient maps.
    ExactFiniteChainComplex,
    /// The relative inclusion is an admitted subcomplex and its quotient
    /// chain complex is exact at the stated window.
    AdmittedRelativeSubcomplex,
}

/// Proof state for one `image(previous) = kernel(next)` obligation.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum ExactnessStatus {
    /// Schema is well typed; constructive witness remains required.
    RequiresConstructiveWitness {
        /// Explicit finite assumption used by the future checker.
        assumption: FiniteExactnessAssumption,
    },
}

/// One exactness obligation at an interior sequence object.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct ExactnessClaim {
    /// Object index where exactness is claimed.
    pub at_object: usize,
    /// Incoming map whose image is compared.
    pub image_of_map: usize,
    /// Outgoing map whose kernel is compared.
    pub kernel_of_map: usize,
    /// Honest proof state.
    pub status: ExactnessStatus,
}

/// Typed exact-sequence schema.  Construction verifies map composability, but
/// does not fabricate algebraic witnesses.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ExactSequenceSchema {
    kind: ExactSequenceKind,
    objects: Vec<ObjectSpace>,
    maps: Vec<TypedMap>,
    exactness: Vec<ExactnessClaim>,
}

impl ExactSequenceSchema {
    fn short_exact(
        kind: ExactSequenceKind,
        objects: Vec<ObjectSpace>,
        map_kinds: [MapKind; 4],
    ) -> Result<Self, CharacterError> {
        Self::from_parts(
            kind,
            objects,
            map_kinds,
            FiniteExactnessAssumption::ExactFiniteChainComplex,
        )
    }

    fn long_window(
        kind: ExactSequenceKind,
        objects: Vec<ObjectSpace>,
        map_kinds: [MapKind; 4],
    ) -> Result<Self, CharacterError> {
        Self::from_parts(
            kind,
            objects,
            map_kinds,
            FiniteExactnessAssumption::AdmittedRelativeSubcomplex,
        )
    }

    fn from_parts(
        kind: ExactSequenceKind,
        objects: Vec<ObjectSpace>,
        map_kinds: [MapKind; 4],
        assumption: FiniteExactnessAssumption,
    ) -> Result<Self, CharacterError> {
        if objects.len() != 5 {
            return Err(CharacterError::SequenceArity);
        }
        let mut maps = Vec::with_capacity(4);
        for (index, map_kind) in map_kinds.into_iter().enumerate() {
            maps.push(TypedMap::new(
                map_kind,
                objects[index].clone(),
                objects[index + 1].clone(),
            )?);
        }
        let exactness = (1..=3)
            .map(|at_object| ExactnessClaim {
                at_object,
                image_of_map: at_object - 1,
                kernel_of_map: at_object,
                status: ExactnessStatus::RequiresConstructiveWitness { assumption },
            })
            .collect();
        Ok(Self {
            kind,
            objects,
            maps,
            exactness,
        })
    }

    /// Sequence family.
    #[must_use]
    pub const fn kind(&self) -> ExactSequenceKind {
        self.kind
    }

    /// Ordered objects.
    #[must_use]
    pub fn objects(&self) -> &[ObjectSpace] {
        &self.objects
    }

    /// Ordered maps between adjacent objects.
    #[must_use]
    pub fn maps(&self) -> &[TypedMap] {
        &self.maps
    }

    /// Explicit image/kernel obligations and their proof state.
    #[must_use]
    pub fn exactness_claims(&self) -> &[ExactnessClaim] {
        &self.exactness
    }

    /// Re-check structural composability.  This is intentionally weaker than
    /// an exactness proof.
    pub fn validate_composable(&self) -> Result<(), CharacterError> {
        if self.maps.len() + 1 != self.objects.len() {
            return Err(CharacterError::SequenceArity);
        }
        for (index, map) in self.maps.iter().enumerate() {
            if map.domain != self.objects[index] || map.codomain != self.objects[index + 1] {
                return Err(CharacterError::SequenceNotComposable { map: index });
            }
        }
        for claim in &self.exactness {
            if claim.at_object == 0
                || claim.at_object + 1 >= self.objects.len()
                || claim.image_of_map + 1 != claim.at_object
                || claim.kernel_of_map != claim.at_object
            {
                return Err(CharacterError::InvalidExactnessClaim);
            }
        }
        Ok(())
    }
}

/// Bilinear operation family.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum BilinearMapKind {
    /// Differential-cohomology cup product.
    CupProduct,
    /// Holonomy evaluation on a relative cycle.
    HolonomyPairing,
}

/// Provenance for the coefficient-level part of a bilinear map.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum BilinearCoefficientRule {
    /// Caller-declared coefficient product; no ring law is inferred.
    Declared(Box<CoefficientProductSchema>),
    /// Canonical evaluation of a lattice differential character on an
    /// integral cycle, with output in `R/Lambda`.
    HolonomyEvaluation,
}

/// Typed bilinear map schema.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct BilinearMapSchema {
    /// Operation semantics.
    kind: BilinearMapKind,
    /// Left input.
    left: ObjectSpace,
    /// Right input.
    right: ObjectSpace,
    /// Output.
    output: ObjectSpace,
    /// Explicit coefficient-level map.
    coefficient_rule: BilinearCoefficientRule,
}

impl BilinearMapSchema {
    /// Operation semantics.
    #[must_use]
    pub const fn kind(&self) -> BilinearMapKind {
        self.kind
    }

    /// Left typed input.
    #[must_use]
    pub const fn left(&self) -> &ObjectSpace {
        &self.left
    }

    /// Right typed input.
    #[must_use]
    pub const fn right(&self) -> &ObjectSpace {
        &self.right
    }

    /// Typed output.
    #[must_use]
    pub const fn output(&self) -> &ObjectSpace {
        &self.output
    }

    /// Explicit coefficient-level rule/provenance.
    #[must_use]
    pub const fn coefficient_rule(&self) -> &BilinearCoefficientRule {
        &self.coefficient_rule
    }
}

/// Refusals from the finite differential-character schema layer.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum CharacterError {
    /// Required identifier was empty.
    EmptyIdentifier { object: &'static str },
    /// Semantic revisions start at one.
    ZeroRevision { object: &'static str },
    /// One or more resource limits were zero.
    ZeroBudget,
    /// Cell-count vector does not cover exactly `0..=dimension`.
    CellCountArity { dimension: u8, actual: usize },
    /// Cell-count sum overflowed.
    CellCountOverflow,
    /// Requested cells exceed the explicit budget.
    CellBudgetExceeded { requested: u64, limit: u64 },
    /// Named boundary count exceeds the explicit budget.
    BoundaryBudgetExceeded,
    /// Canonical serialization exceeds the explicit budget.
    CanonicalBudgetExceeded { requested: u64, limit: u64 },
    /// Relative complex is higher-dimensional than its ambient complex.
    RelativeDimension { ambient: u8, relative: u8 },
    /// The degree representation cannot encode `dim(X) + 1`.
    CharacterDegreeRepresentationOverflow { dimension: u8 },
    /// Relative cell count exceeds the ambient count.
    RelativeCellCount {
        degree: u8,
        ambient: u64,
        relative: u64,
    },
    /// Primal/dual confusion.
    LaneMismatch {
        expected: ComplexLane,
        actual: ComplexLane,
        object: &'static str,
    },
    /// Named boundary exceeds the relative subcomplex dimension.
    BoundaryDimension { boundary: u8, relative: u8 },
    /// A terminal is not codimension one in the ambient complex.
    TerminalCodimension { ambient: u8, terminal: u8 },
    /// Boundary identities must be unique.
    DuplicateBoundary { id: String },
    /// Lattice rank is zero or too large.
    InvalidLatticeRank { rank: usize },
    /// Lattice scale is non-finite or non-positive.
    InvalidLatticeScale { index: usize },
    /// `Z/n` requires `n >= 2`.
    InvalidCyclicModulus { modulus: u64 },
    /// Object degree exceeds the representative-specific maximum.
    DegreeOutOfRange { degree: u8, maximum: u8 },
    /// Representative and coefficient sectors disagree.
    RepresentativeCoefficientMismatch {
        representative: RepresentativeKind,
        expected: CoefficientSector,
        actual: CoefficientSector,
    },
    /// Degree-zero character has no boundary `(k - 1)` trivialization.
    DegreeZeroHasNoBoundaryTrivialization,
    /// Degree-zero character has no gauge `(k - 1)` parameter.
    DegreeZeroHasNoGaugeParameter,
    /// Degree-zero character does not instantiate these short sequences.
    DegreeZeroHasNoExactSequence,
    /// Degree-zero character has no holonomy cycle.
    DegreeZeroHasNoHolonomyCycle,
    /// No in-range successor exists for this finite complex.
    DegreeHasNoSuccessorInComplex { degree: u8, maximum: u8 },
    /// Terminal names an unknown boundary.
    UnknownBoundary { id: String },
    /// Trivialization names an interface rather than relative/terminal data.
    BoundaryDoesNotAdmitTrivialization { id: String },
    /// Relative/terminal trivialization has the wrong degree.
    BoundaryTrivializationDegreeMismatch { expected: u8, actual: u8 },
    /// A terminal has no required trivialization.
    MissingTerminalTrivialization { id: String },
    /// Nonempty `A` has no whole-subcomplex mapping-cone trivialization.
    MissingRelativeSubcomplexTrivialization,
    /// Empty `A` was given a spurious mapping-cone trivialization.
    UnexpectedRelativeSubcomplexTrivialization,
    /// A boundary was trivialized more than once.
    DuplicateBoundaryTrivialization { id: String },
    /// Coefficient systems disagree.
    CoefficientMismatch { object: &'static str },
    /// Pair identities disagree.
    PairMismatch,
    /// Torsion coefficients have no nonzero real curvature lane.
    TorsionHasNoCurvature,
    /// An integral cocycle alone is not a differential-curvature object.
    IntegralCocycleHasNoCurvature,
    /// An ordinary real cochain alone is not a differential-curvature object.
    RealCochainHasNoCurvature,
    /// Plain real coefficients carry no integral characteristic class.
    RealHasNoCharacteristicClass,
    /// A curvature coefficient space is a codomain, not a character.
    CurvatureSpaceHasNoCharacteristicClass,
    /// Operation is defined only for a Hopkins--Singer differential character,
    /// not for a supporting raw cocycle/cochain schema.
    OperationRequiresDifferentialCharacter,
    /// Standard differential-character sequences require `R/Lambda`.
    StandardSequenceRequiresRealModuloLattice,
    /// Cup-product degree exceeds the differential-character degree ceiling.
    ProductDegreeOutOfRange,
    /// Declared coefficient product does not accept the supplied input.
    CoefficientProductMismatch { input: &'static str },
    /// Differential-character holonomy requires `R/Lambda` coefficients.
    HolonomyRequiresRealModuloLattice,
    /// Exact-sequence object/map count is malformed.
    SequenceArity,
    /// A map does not connect the adjacent declared objects.
    SequenceNotComposable { map: usize },
    /// Image/kernel indices are malformed.
    InvalidExactnessClaim,
}

impl fmt::Display for CharacterError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyIdentifier { object } => write!(f, "{object} must not be empty"),
            Self::ZeroRevision { object } => write!(f, "{object} revision must be nonzero"),
            Self::ZeroBudget => f.write_str("algebra budgets must be nonzero"),
            Self::CellCountArity { dimension, actual } => write!(
                f,
                "dimension {dimension} requires {} cell counts, got {actual}",
                usize::from(*dimension) + 1
            ),
            Self::CellCountOverflow => f.write_str("finite-complex cell count overflow"),
            Self::CellBudgetExceeded { requested, limit } => {
                write!(
                    f,
                    "cell budget exceeded: requested {requested}, limit {limit}"
                )
            }
            Self::BoundaryBudgetExceeded => f.write_str("boundary-component budget exceeded"),
            Self::CanonicalBudgetExceeded { requested, limit } => write!(
                f,
                "canonical serialization budget exceeded: requested {requested}, limit {limit}"
            ),
            Self::RelativeDimension { ambient, relative } => write!(
                f,
                "relative dimension {relative} exceeds ambient dimension {ambient}"
            ),
            Self::CharacterDegreeRepresentationOverflow { dimension } => write!(
                f,
                "ambient dimension {dimension} cannot represent the required dim(X) + 1 character degree"
            ),
            Self::RelativeCellCount {
                degree,
                ambient,
                relative,
            } => write!(
                f,
                "relative degree-{degree} cell count {relative} exceeds ambient {ambient}"
            ),
            Self::LaneMismatch {
                expected,
                actual,
                object,
            } => write!(f, "{object} is in {actual:?} lane; expected {expected:?}"),
            Self::BoundaryDimension { boundary, relative } => write!(
                f,
                "boundary dimension {boundary} exceeds relative dimension {relative}"
            ),
            Self::TerminalCodimension { ambient, terminal } => write!(
                f,
                "terminal dimension {terminal} is not codimension one in ambient dimension {ambient}"
            ),
            Self::DuplicateBoundary { id } => write!(f, "duplicate boundary '{id}'"),
            Self::InvalidLatticeRank { rank } => write!(f, "invalid lattice rank {rank}"),
            Self::InvalidLatticeScale { index } => {
                write!(
                    f,
                    "lattice generator scale {index} must be positive and finite"
                )
            }
            Self::InvalidCyclicModulus { modulus } => {
                write!(f, "cyclic modulus must be at least 2, got {modulus}")
            }
            Self::DegreeOutOfRange { degree, maximum } => write!(
                f,
                "object degree {degree} exceeds representative maximum {maximum}"
            ),
            Self::RepresentativeCoefficientMismatch {
                representative,
                expected,
                actual,
            } => write!(
                f,
                "representative {representative:?} requires {expected:?} coefficients, got {actual:?}"
            ),
            Self::DegreeZeroHasNoBoundaryTrivialization => {
                f.write_str("degree-zero character has no degree-minus-one boundary trivialization")
            }
            Self::DegreeZeroHasNoGaugeParameter => {
                f.write_str("degree-zero character has no degree-minus-one gauge parameter")
            }
            Self::DegreeZeroHasNoExactSequence => {
                f.write_str("degree-zero character does not instantiate this sequence window")
            }
            Self::DegreeZeroHasNoHolonomyCycle => {
                f.write_str("degree-zero character has no degree-minus-one holonomy cycle")
            }
            Self::DegreeHasNoSuccessorInComplex { degree, maximum } => write!(
                f,
                "degree {degree} has no successor within object maximum {maximum}"
            ),
            Self::UnknownBoundary { id } => write!(f, "unknown boundary '{id}'"),
            Self::BoundaryDoesNotAdmitTrivialization { id } => {
                write!(
                    f,
                    "interface boundary '{id}' does not admit relative trivialization"
                )
            }
            Self::BoundaryTrivializationDegreeMismatch { expected, actual } => write!(
                f,
                "boundary trivialization degree {actual}; expected {expected}"
            ),
            Self::MissingTerminalTrivialization { id } => {
                write!(f, "terminal '{id}' has no trivialization")
            }
            Self::MissingRelativeSubcomplexTrivialization => f.write_str(
                "mapping-cone object on nonempty A requires a whole-subcomplex trivialization",
            ),
            Self::UnexpectedRelativeSubcomplexTrivialization => f.write_str(
                "mapping-cone object on empty A must not carry a relative trivialization",
            ),
            Self::DuplicateBoundaryTrivialization { id } => {
                write!(f, "boundary '{id}' has duplicate trivializations")
            }
            Self::CoefficientMismatch { object } => {
                write!(f, "coefficient mismatch in {object}")
            }
            Self::PairMismatch => f.write_str("relative-pair identity mismatch"),
            Self::TorsionHasNoCurvature => {
                f.write_str("finite torsion coefficients have no real curvature map")
            }
            Self::IntegralCocycleHasNoCurvature => {
                f.write_str("an integral cocycle alone has no differential curvature map")
            }
            Self::RealCochainHasNoCurvature => {
                f.write_str("an ordinary real cochain alone has no differential curvature map")
            }
            Self::RealHasNoCharacteristicClass => {
                f.write_str("plain real coefficients have no integral characteristic class")
            }
            Self::CurvatureSpaceHasNoCharacteristicClass => {
                f.write_str("a curvature coefficient space has no characteristic-class map")
            }
            Self::OperationRequiresDifferentialCharacter => f.write_str(
                "operation requires a Hopkins-Singer differential character, not a raw cocycle",
            ),
            Self::StandardSequenceRequiresRealModuloLattice => {
                f.write_str("standard character sequence requires R/Lambda coefficients")
            }
            Self::ProductDegreeOutOfRange => {
                f.write_str("cup-product degree exceeds the differential-character ceiling")
            }
            Self::CoefficientProductMismatch { input } => {
                write!(f, "declared coefficient product rejects the {input} input")
            }
            Self::HolonomyRequiresRealModuloLattice => {
                f.write_str("holonomy pairing requires R/Lambda coefficients")
            }
            Self::SequenceArity => f.write_str("exact-sequence arity is malformed"),
            Self::SequenceNotComposable { map } => {
                write!(f, "exact-sequence map {map} is not composable")
            }
            Self::InvalidExactnessClaim => {
                f.write_str("exactness image/kernel indices are invalid")
            }
        }
    }
}

impl std::error::Error for CharacterError {}

fn nonempty(value: String, object: &'static str) -> Result<String, CharacterError> {
    if value.trim().is_empty() {
        Err(CharacterError::EmptyIdentifier { object })
    } else {
        Ok(value)
    }
}

fn usize_to_u64(value: usize) -> u64 {
    u64::try_from(value).unwrap_or(u64::MAX)
}

fn usize_to_u8(value: usize) -> u8 {
    u8::try_from(value).unwrap_or(u8::MAX)
}

struct CanonicalEncoder {
    output: Vec<u8>,
    limit: u64,
    logical_len: u64,
    overflow: Option<u64>,
}

impl CanonicalEncoder {
    fn new(domain: &[u8], limit: u64) -> Self {
        let mut encoder = Self {
            output: Vec::new(),
            limit,
            logical_len: 0,
            overflow: None,
        };
        encoder.bytes(domain);
        encoder
    }

    fn finish(self) -> Result<Vec<u8>, CharacterError> {
        match self.overflow {
            Some(requested) => Err(CharacterError::CanonicalBudgetExceeded {
                requested,
                limit: self.limit,
            }),
            None => Ok(self.output),
        }
    }

    fn bytes(&mut self, value: &[u8]) {
        self.u64(usize_to_u64(value.len()));
        self.push(value);
    }

    fn string(&mut self, value: &str) {
        self.bytes(value.as_bytes());
    }

    fn u8(&mut self, value: u8) {
        self.push(&[value]);
    }

    fn u16(&mut self, value: u16) {
        self.push(&value.to_le_bytes());
    }

    fn u32(&mut self, value: u32) {
        self.push(&value.to_le_bytes());
    }

    fn u64(&mut self, value: u64) {
        self.push(&value.to_le_bytes());
    }

    fn dims(&mut self, value: Dims) {
        for exponent in value.0 {
            self.u8(exponent.to_le_bytes()[0]);
        }
    }

    fn lane(&mut self, value: ComplexLane) {
        self.u8(match value {
            ComplexLane::Primal => 0,
            ComplexLane::Dual => 1,
        });
    }

    fn orientation(&mut self, value: Orientation) {
        self.u8(match value {
            Orientation::Positive => 0,
            Orientation::Negative => 1,
        });
    }

    fn boundary_orientation(&mut self, value: BoundaryOrientation) {
        self.u8(match value {
            BoundaryOrientation::Induced => 0,
            BoundaryOrientation::Opposite => 1,
        });
    }

    fn boundary_role(&mut self, value: BoundaryRole) {
        self.u8(match value {
            BoundaryRole::Relative => 0,
            BoundaryRole::Terminal => 1,
            BoundaryRole::Interface => 2,
        });
    }

    fn relative_model(&mut self, value: RelativeModel) {
        self.u8(match value {
            RelativeModel::MappingCone => 0,
        });
    }

    fn representative(&mut self, value: RepresentativeKind) {
        self.u8(match value {
            RepresentativeKind::HopkinsSingerTriple => 0,
            RepresentativeKind::IntegralCocycle => 1,
            RepresentativeKind::FlatTorsionCocycle => 2,
            RepresentativeKind::RealCochain => 3,
        });
    }

    fn budget(&mut self, value: AlgebraBudget) {
        self.u64(value.max_cells);
        self.u32(value.max_boundary_components);
        self.u64(value.max_canonical_bytes);
    }

    fn push(&mut self, value: &[u8]) {
        let requested = self
            .logical_len
            .checked_add(usize_to_u64(value.len()))
            .unwrap_or(u64::MAX);
        self.logical_len = requested;
        if self.overflow.is_some() {
            self.overflow = Some(requested);
            return;
        }
        if requested > self.limit {
            self.overflow = Some(requested);
            return;
        }
        self.output.extend_from_slice(value);
    }
}

fn stable_digest(bytes: &[u8]) -> AlgebraId {
    const PRIME: u64 = 0x0000_0100_0000_01b3;
    let mut lanes = [
        0xcbf2_9ce4_8422_2325,
        0x8422_2325_cbf2_9ce4,
        0x9e37_79b9_7f4a_7c15,
        0xd6e8_feb8_6659_fd93,
    ];
    for (index, byte) in bytes.iter().copied().enumerate() {
        let position = usize_to_u64(index);
        for (lane_index, lane) in lanes.iter_mut().enumerate() {
            *lane ^= u64::from(byte)
                .wrapping_add(position.rotate_left(u32::try_from(lane_index * 11).unwrap_or(0)));
            *lane = lane.wrapping_mul(PRIME);
            *lane ^= lane.rotate_right(u32::try_from(13 + lane_index * 7).unwrap_or(13));
        }
    }
    let mut digest = [0_u8; 32];
    for (index, lane) in lanes.into_iter().enumerate() {
        digest[index * 8..(index + 1) * 8].copy_from_slice(&lane.to_le_bytes());
    }
    AlgebraId(digest)
}
