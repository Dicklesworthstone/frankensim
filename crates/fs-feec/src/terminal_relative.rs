//! Typed terminal-relative cellular schema objects for winding topology (I13.2a).
//!
//! This module is a schema and exact-incidence boundary, not a homology
//! solver.  It keeps four physically different things nominally separate:
//! declared integral winding representatives, real current amplitudes, distributed real
//! current cochains, and geometric coil realizations.  Any bridge between
//! them must retain an explicit map artifact.

use core::fmt;

use std::collections::{BTreeMap, BTreeSet};

use fs_blake3::identity::{
    CanonicalEncoder, CanonicalError, CanonicalLimits, CanonicalSchema, ChildSpec, Field,
    FieldSpec, IdentityReceipt, NeverCancel, SemanticId, StrongIdentity, WireType,
};
use fs_couple::{
    ConservationRole, FieldMeasureSide, PortKind, PortOrientation, PortSchema, PortValueShape,
    PowerPairing, StableId,
};
use fs_qty::{Current, Dims};

/// Canonical schema version for terminal-relative pairs and winding representatives.
pub const TERMINAL_RELATIVE_SCHEMA_VERSION: u32 = 1;
/// Exact L6 MachineGraph identity domain accepted as presented data.
pub const PRESENTED_MACHINE_GRAPH_DOMAIN: &str = "org.frankensim.fs-ir.machine.graph.v1";
/// Exact L6 MachineGraph schema version accepted as presented data.
pub const PRESENTED_MACHINE_GRAPH_SCHEMA_VERSION: u32 = 1;
/// Largest admitted cell-complex dimension in this physical lane.
pub const MAX_TERMINAL_RELATIVE_DIMENSION: u8 = 3;
/// Maximum total cells in one admitted complex.
pub const MAX_TERMINAL_RELATIVE_CELLS: usize = 131_072;
/// Maximum nonzero boundary incidences in one admitted complex.
pub const MAX_TERMINAL_RELATIVE_INCIDENCES: usize = 1_048_576;
/// Maximum conductor components in one pair.
pub const MAX_CONDUCTOR_COMPONENTS: usize = 4_096;
/// Maximum physical terminals in one pair.
pub const MAX_PHYSICAL_TERMINALS: usize = 4_096;
/// Maximum canonical payload size before the typed identity frame is added.
pub const MAX_TERMINAL_RELATIVE_CANONICAL_BYTES: usize = 2 * 1_024 * 1_024;

const PAIR_IDENTITY_LIMITS: CanonicalLimits =
    CanonicalLimits::new(4 * 1_024 * 1_024, 2 * 1_024 * 1_024, 1, 1, 256);
const SIGNED_RELABEL_IDENTITY_LIMITS: CanonicalLimits =
    CanonicalLimits::new(4 * 1_024 * 1_024, 2 * 1_024 * 1_024, 3, 1, 256);
const PHYSICAL_RELABEL_IDENTITY_LIMITS: CanonicalLimits =
    CanonicalLimits::new(12 * 1_024 * 1_024, 2 * 1_024 * 1_024, 6, 1, 256);

/// Strong semantic identity of one admitted physical terminal-relative pair.
pub enum TerminalRelativePairSchemaV1 {}

impl CanonicalSchema for TerminalRelativePairSchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-feec.terminal-relative-pair.v1";
    const NAME: &'static str = "terminal-relative-pair";
    const VERSION: u32 = TERMINAL_RELATIVE_SCHEMA_VERSION;
    const CONTEXT: &'static str =
        "physical conductor, terminal, insulation, component, phase, and port pair";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required(
        "terminal-relative-payload",
        WireType::Bytes,
    )];
}

/// Nominal pair identity.  It cannot be confused with a representative ID.
pub type TerminalRelativePairId = SemanticId<TerminalRelativePairSchemaV1>;

static TERMINAL_RELATIVE_SIGNED_RELABEL_PAIR_CHILD: ChildSpec =
    ChildSpec::for_identity::<TerminalRelativePairId>();

/// Strong semantic identity of one admitted orientation-aware cell relabeling.
pub enum TerminalRelativeSignedRelabelSchemaV1 {}

impl CanonicalSchema for TerminalRelativeSignedRelabelSchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-feec.terminal-relative-signed-relabel.v1";
    const NAME: &'static str = "terminal-relative-signed-relabel";
    const VERSION: u32 = TERMINAL_RELATIVE_SCHEMA_VERSION;
    const CONTEXT: &'static str =
        "complete signed cell bijection preserving incidence and terminal-relative semantics";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::child_of("source-pair", &TERMINAL_RELATIVE_SIGNED_RELABEL_PAIR_CHILD),
        FieldSpec::child_of("target-pair", &TERMINAL_RELATIVE_SIGNED_RELABEL_PAIR_CHILD),
        FieldSpec::required("signed-cell-map", WireType::Bytes),
    ];
}

/// Nominal identity of one admitted terminal-relative signed relabeling.
pub type TerminalRelativeSignedRelabelId = SemanticId<TerminalRelativeSignedRelabelSchemaV1>;

/// Strong semantic identity of one admitted physical terminal-relative
/// relabeling with explicit semantic permutations and phase-current action.
pub enum TerminalRelativePhysicalRelabelSchemaV1 {}

impl CanonicalSchema for TerminalRelativePhysicalRelabelSchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-feec.terminal-relative-physical-relabel.v1";
    const NAME: &'static str = "terminal-relative-physical-relabel";
    const VERSION: u32 = TERMINAL_RELATIVE_SCHEMA_VERSION;
    const CONTEXT: &'static str = "complete signed cell bijection with explicit component, phase, terminal, and phase-current transport";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::child_of("source-pair", &TERMINAL_RELATIVE_SIGNED_RELABEL_PAIR_CHILD),
        FieldSpec::child_of("target-pair", &TERMINAL_RELATIVE_SIGNED_RELABEL_PAIR_CHILD),
        FieldSpec::required("signed-cell-map", WireType::Bytes),
        FieldSpec::required("component-map", WireType::Bytes),
        FieldSpec::required("phase-map", WireType::Bytes),
        FieldSpec::required("terminal-map", WireType::Bytes),
    ];
}

/// Nominal identity of one admitted physical terminal-relative relabeling.
pub type TerminalRelativePhysicalRelabelId = SemanticId<TerminalRelativePhysicalRelabelSchemaV1>;

/// Strong semantic identity of one declared integral winding representative.
pub enum IntegralWindingRepresentativeSchemaV1 {}

impl CanonicalSchema for IntegralWindingRepresentativeSchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-feec.integral-winding-representative.v1";
    const NAME: &'static str = "integral-winding-representative";
    const VERSION: u32 = TERMINAL_RELATIVE_SCHEMA_VERSION;
    const CONTEXT: &'static str =
        "declared integral relative one-cycle bound to a phase and pair; no quotient-class claim";
    const FIELDS: &'static [FieldSpec] = &[FieldSpec::required(
        "winding-representative-payload",
        WireType::Bytes,
    )];
}

/// Nominal identity of a declared integral winding representative.
pub type IntegralWindingRepresentativeId = SemanticId<IntegralWindingRepresentativeSchemaV1>;

macro_rules! typed_stable_id {
    ($name:ident, $role:literal) => {
        #[doc = concat!("Nominal stable identity for one ", $role, ".")]
        #[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name(StableId);

        impl $name {
            #[doc = concat!("Construct a canonical ", $role, " identity.")]
            pub fn new(value: impl Into<String>) -> Result<Self, TerminalRelativeError> {
                let value = value.into();
                StableId::new(value.clone())
                    .map(Self)
                    .map_err(|_| TerminalRelativeError::InvalidIdentity { role: $role, value })
            }

            /// Canonical text carried by this nominal identity.
            #[must_use]
            pub fn as_str(&self) -> &str {
                self.0.as_str()
            }
        }
    };
}

typed_stable_id!(ConductorComponentId, "conductor component");
typed_stable_id!(PhysicalTerminalId, "physical terminal");
typed_stable_id!(PhaseId, "electrical phase");
typed_stable_id!(TrivializationId, "terminal trivialization");
typed_stable_id!(PhysicalObjectId, "physical object");
typed_stable_id!(ConversionMapId, "physical conversion map");

/// Canonical reference to one oriented cell.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct CellRef {
    degree: u8,
    ordinal: u32,
}

impl CellRef {
    /// Construct a cell reference.  Extent validation occurs at complex
    /// admission so references remain cheap value types.
    #[must_use]
    pub const fn new(degree: u8, ordinal: u32) -> Self {
        Self { degree, ordinal }
    }

    /// Cell dimension/chain degree.
    #[must_use]
    pub const fn degree(self) -> u8 {
        self.degree
    }

    /// Canonical ordinal within this degree.
    #[must_use]
    pub const fn ordinal(self) -> u32 {
        self.ordinal
    }
}

/// Exact orientation coefficient of one cellular incidence.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(i8)]
pub enum IncidenceSign {
    /// Negative orientation coefficient.
    Negative = -1,
    /// Positive orientation coefficient.
    Positive = 1,
}

impl IncidenceSign {
    const fn as_i128(self) -> i128 {
        match self {
            Self::Negative => -1,
            Self::Positive => 1,
        }
    }

    const fn tag(self) -> u8 {
        match self {
            Self::Negative => 0,
            Self::Positive => 1,
        }
    }
}

/// One oriented source-to-target cell row in a complete signed relabeling.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SignedCellRelabelEntry {
    source: CellRef,
    target: CellRef,
    sign: IncidenceSign,
}

impl SignedCellRelabelEntry {
    /// Declare that the oriented source cell maps to `sign * target`.
    #[must_use]
    pub const fn new(source: CellRef, target: CellRef, sign: IncidenceSign) -> Self {
        Self {
            source,
            target,
            sign,
        }
    }

    /// Source-space cell.
    #[must_use]
    pub const fn source(self) -> CellRef {
        self.source
    }

    /// Target-space cell.
    #[must_use]
    pub const fn target(self) -> CellRef {
        self.target
    }

    /// Orientation action on the source basis cell.
    #[must_use]
    pub const fn sign(self) -> IncidenceSign {
        self.sign
    }
}

/// Explicit action on the positive-current coordinate of one mapped phase.
///
/// This is deliberately distinct from [`OrientationMapSign`], which describes
/// one terminal's port trivialization rather than a change of physical phase
/// convention.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum PhaseCurrentSign {
    /// Preserve the positive-current coordinate.
    Preserve,
    /// Negate the positive-current coordinate.
    Reverse,
}

impl PhaseCurrentSign {
    /// Apply this action to a finite real current amplitude.
    #[must_use]
    pub fn apply(self, value: Current) -> Current {
        match self {
            Self::Preserve => value,
            Self::Reverse => Current::new(-value.value()),
        }
    }

    /// Compose `next` after this action.
    #[must_use]
    pub const fn compose(self, next: Self) -> Self {
        match (self, next) {
            (Self::Preserve, Self::Preserve) | (Self::Reverse, Self::Reverse) => Self::Preserve,
            (Self::Preserve, Self::Reverse) | (Self::Reverse, Self::Preserve) => Self::Reverse,
        }
    }

    /// Every sign action is its own inverse.
    #[must_use]
    pub const fn inverse(self) -> Self {
        self
    }

    const fn tag(self) -> u8 {
        match self {
            Self::Preserve => 0,
            Self::Reverse => 1,
        }
    }

    const fn incidence_sign(self) -> IncidenceSign {
        match self {
            Self::Preserve => IncidenceSign::Positive,
            Self::Reverse => IncidenceSign::Negative,
        }
    }

    const fn map_role(self, role: TerminalRole) -> TerminalRole {
        match (self, role) {
            (Self::Preserve, role) => role,
            (Self::Reverse, TerminalRole::Driven) => TerminalRole::ReturnReference,
            (Self::Reverse, TerminalRole::ReturnReference) => TerminalRole::Driven,
        }
    }

    const fn map_orientation(self, orientation: TerminalOrientation) -> TerminalOrientation {
        match (self, orientation) {
            (Self::Preserve, orientation) => orientation,
            (Self::Reverse, TerminalOrientation::IntoConductor) => {
                TerminalOrientation::OutOfConductor
            }
            (Self::Reverse, TerminalOrientation::OutOfConductor) => {
                TerminalOrientation::IntoConductor
            }
        }
    }

    const fn map_trivialization(self, sign: OrientationMapSign) -> OrientationMapSign {
        match (self, sign) {
            (Self::Preserve, sign) => sign,
            (Self::Reverse, OrientationMapSign::Preserve) => OrientationMapSign::Reverse,
            (Self::Reverse, OrientationMapSign::Reverse) => OrientationMapSign::Preserve,
        }
    }
}

/// One explicitly declared source-to-target conductor-component row.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct ComponentRelabelEntry {
    source: ConductorComponentId,
    target: ConductorComponentId,
}

impl ComponentRelabelEntry {
    /// Declare one component image.
    #[must_use]
    pub const fn new(source: ConductorComponentId, target: ConductorComponentId) -> Self {
        Self { source, target }
    }

    /// Source component identity.
    #[must_use]
    pub const fn source(&self) -> &ConductorComponentId {
        &self.source
    }

    /// Target component identity.
    #[must_use]
    pub const fn target(&self) -> &ConductorComponentId {
        &self.target
    }
}

/// One explicitly declared source-to-target phase row and current action.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct PhaseRelabelEntry {
    source: PhaseId,
    target: PhaseId,
    current_sign: PhaseCurrentSign,
}

impl PhaseRelabelEntry {
    /// Declare one phase image and its physical positive-current action.
    #[must_use]
    pub const fn new(source: PhaseId, target: PhaseId, current_sign: PhaseCurrentSign) -> Self {
        Self {
            source,
            target,
            current_sign,
        }
    }

    /// Source phase identity.
    #[must_use]
    pub const fn source(&self) -> &PhaseId {
        &self.source
    }

    /// Target phase identity.
    #[must_use]
    pub const fn target(&self) -> &PhaseId {
        &self.target
    }

    /// Explicit positive-current action.
    #[must_use]
    pub const fn current_sign(&self) -> PhaseCurrentSign {
        self.current_sign
    }
}

/// One explicitly declared source-to-target physical-terminal row.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct TerminalRelabelEntry {
    source: PhysicalTerminalId,
    target: PhysicalTerminalId,
}

impl TerminalRelabelEntry {
    /// Declare one physical-terminal image.
    #[must_use]
    pub const fn new(source: PhysicalTerminalId, target: PhysicalTerminalId) -> Self {
        Self { source, target }
    }

    /// Source terminal identity.
    #[must_use]
    pub const fn source(&self) -> &PhysicalTerminalId {
        &self.source
    }

    /// Target terminal identity.
    #[must_use]
    pub const fn target(&self) -> &PhysicalTerminalId {
        &self.target
    }
}

/// Complete explicit semantic permutation carried by a physical relabeling.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TerminalRelativeSemanticPermutation {
    components: Vec<ComponentRelabelEntry>,
    phases: Vec<PhaseRelabelEntry>,
    terminals: Vec<TerminalRelabelEntry>,
}

impl TerminalRelativeSemanticPermutation {
    /// Declare the three total semantic maps. Admission and canonical sorting
    /// occur only with both endpoint pairs in
    /// [`TerminalRelativePhysicalRelabel::try_new`].
    #[must_use]
    pub const fn new(
        components: Vec<ComponentRelabelEntry>,
        phases: Vec<PhaseRelabelEntry>,
        terminals: Vec<TerminalRelabelEntry>,
    ) -> Self {
        Self {
            components,
            phases,
            terminals,
        }
    }

    /// Canonically source-sorted component rows after admission.
    #[must_use]
    pub fn components(&self) -> &[ComponentRelabelEntry] {
        &self.components
    }

    /// Canonically source-sorted phase rows after admission.
    #[must_use]
    pub fn phases(&self) -> &[PhaseRelabelEntry] {
        &self.phases
    }

    /// Canonically source-sorted terminal rows after admission.
    #[must_use]
    pub fn terminals(&self) -> &[TerminalRelabelEntry] {
        &self.terminals
    }
}

/// One nonzero entry in an exact integer boundary matrix.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct BoundaryIncidence {
    lower: CellRef,
    upper: CellRef,
    sign: IncidenceSign,
}

impl BoundaryIncidence {
    /// Declare `sign * lower` in the boundary of `upper`.
    #[must_use]
    pub const fn new(lower: CellRef, upper: CellRef, sign: IncidenceSign) -> Self {
        Self { lower, upper, sign }
    }

    /// Lower-dimensional face.
    #[must_use]
    pub const fn lower(self) -> CellRef {
        self.lower
    }

    /// Higher-dimensional cell.
    #[must_use]
    pub const fn upper(self) -> CellRef {
        self.upper
    }

    /// Exact orientation coefficient.
    #[must_use]
    pub const fn sign(self) -> IncidenceSign {
        self.sign
    }
}

/// Bounded oriented cell complex with exact admitted boundary matrices.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FiniteCellComplex {
    dimension: u8,
    cell_counts: Vec<u32>,
    incidences: Vec<BoundaryIncidence>,
}

impl FiniteCellComplex {
    /// Admit exact incidence after range, degree, uniqueness, budget, and
    /// `boundary * boundary = 0` checks.
    pub fn try_new(
        dimension: u8,
        cell_counts: Vec<u32>,
        mut incidences: Vec<BoundaryIncidence>,
    ) -> Result<Self, TerminalRelativeError> {
        if dimension > MAX_TERMINAL_RELATIVE_DIMENSION {
            return Err(TerminalRelativeError::DimensionTooLarge {
                actual: dimension,
                max: MAX_TERMINAL_RELATIVE_DIMENSION,
            });
        }
        if cell_counts.len() != usize::from(dimension) + 1 {
            return Err(TerminalRelativeError::CellCountArity {
                dimension,
                actual: cell_counts.len(),
            });
        }
        let total_cells = cell_counts.iter().try_fold(0_usize, |sum, count| {
            sum.checked_add(usize::try_from(*count).unwrap_or(usize::MAX))
        });
        let Some(total_cells) = total_cells else {
            return Err(TerminalRelativeError::CellBudgetExceeded {
                actual: usize::MAX,
                max: MAX_TERMINAL_RELATIVE_CELLS,
            });
        };
        if total_cells == 0 || total_cells > MAX_TERMINAL_RELATIVE_CELLS {
            return Err(TerminalRelativeError::CellBudgetExceeded {
                actual: total_cells,
                max: MAX_TERMINAL_RELATIVE_CELLS,
            });
        }
        if incidences.len() > MAX_TERMINAL_RELATIVE_INCIDENCES {
            return Err(TerminalRelativeError::IncidenceBudgetExceeded {
                actual: incidences.len(),
                max: MAX_TERMINAL_RELATIVE_INCIDENCES,
            });
        }

        for incidence in &incidences {
            validate_cell_ref(incidence.lower, &cell_counts)?;
            validate_cell_ref(incidence.upper, &cell_counts)?;
            if incidence.lower.degree.checked_add(1) != Some(incidence.upper.degree) {
                return Err(TerminalRelativeError::InvalidIncidenceDegree {
                    lower: incidence.lower,
                    upper: incidence.upper,
                });
            }
        }
        incidences.sort_unstable_by_key(|entry| {
            (entry.upper.degree, entry.upper.ordinal, entry.lower.ordinal)
        });
        for pair in incidences.windows(2) {
            if pair[0].upper == pair[1].upper && pair[0].lower == pair[1].lower {
                return Err(TerminalRelativeError::DuplicateIncidence {
                    lower: pair[0].lower,
                    upper: pair[0].upper,
                });
            }
        }

        let complex = Self {
            dimension,
            cell_counts,
            incidences,
        };
        complex.verify_boundary_squared()?;
        Ok(complex)
    }

    /// Top cell degree.
    #[must_use]
    pub const fn dimension(&self) -> u8 {
        self.dimension
    }

    /// Number of cells at each degree.
    #[must_use]
    pub fn cell_counts(&self) -> &[u32] {
        &self.cell_counts
    }

    /// Canonically sorted nonzero incidence entries.
    #[must_use]
    pub fn incidences(&self) -> &[BoundaryIncidence] {
        &self.incidences
    }

    fn contains(&self, cell: CellRef) -> bool {
        self.cell_counts
            .get(usize::from(cell.degree))
            .is_some_and(|count| cell.ordinal < *count)
    }

    fn verify_boundary_squared(&self) -> Result<(), TerminalRelativeError> {
        let mut by_upper = BTreeMap::<CellRef, Vec<BoundaryIncidence>>::new();
        for incidence in &self.incidences {
            by_upper
                .entry(incidence.upper)
                .or_default()
                .push(*incidence);
        }
        for degree in 2..=self.dimension {
            let upper_count = self.cell_counts[usize::from(degree)];
            for ordinal in 0..upper_count {
                let source = CellRef::new(degree, ordinal);
                let mut accumulated = BTreeMap::<CellRef, i128>::new();
                for upper_to_middle in by_upper.get(&source).into_iter().flatten() {
                    for middle_to_lower in
                        by_upper.get(&upper_to_middle.lower).into_iter().flatten()
                    {
                        let contribution =
                            upper_to_middle.sign.as_i128() * middle_to_lower.sign.as_i128();
                        *accumulated.entry(middle_to_lower.lower).or_default() += contribution;
                    }
                }
                if let Some((target, value)) =
                    accumulated.into_iter().find(|(_, value)| *value != 0)
                {
                    return Err(TerminalRelativeError::BoundarySquaredNonzero {
                        source,
                        target,
                        value,
                    });
                }
            }
        }
        Ok(())
    }
}

fn validate_cell_ref(cell: CellRef, counts: &[u32]) -> Result<(), TerminalRelativeError> {
    let Some(count) = counts.get(usize::from(cell.degree)) else {
        return Err(TerminalRelativeError::CellOutOfRange { cell });
    };
    if cell.ordinal >= *count {
        return Err(TerminalRelativeError::CellOutOfRange { cell });
    }
    Ok(())
}

fn validate_subcomplex_against(
    subcomplex: &CellularSubcomplex,
    ambient: &FiniteCellComplex,
) -> Result<(), TerminalRelativeError> {
    for cell in &subcomplex.cells {
        if !ambient.contains(*cell) {
            return Err(TerminalRelativeError::CellOutOfRange { cell: *cell });
        }
        for incidence in ambient
            .incidences
            .iter()
            .filter(|entry| entry.upper == *cell)
        {
            if !subcomplex.cells.contains(&incidence.lower) {
                return Err(TerminalRelativeError::NotASubcomplex {
                    id: subcomplex.id.as_str().to_owned(),
                    cell: *cell,
                    missing_boundary: incidence.lower,
                });
            }
        }
    }
    Ok(())
}

fn downward_closure(
    seeds: impl IntoIterator<Item = CellRef>,
    boundary_faces: &BTreeMap<CellRef, Vec<CellRef>>,
) -> BTreeSet<CellRef> {
    let mut closure = BTreeSet::new();
    let mut pending = Vec::new();
    for seed in seeds {
        if closure.insert(seed) {
            pending.push(seed);
        }
    }
    while let Some(cell) = pending.pop() {
        for face in boundary_faces.get(&cell).into_iter().flatten() {
            if closure.insert(*face) {
                pending.push(*face);
            }
        }
    }
    closure
}

/// A named, explicitly enumerated cellular subcomplex.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct CellularSubcomplex {
    id: StableId,
    cells: BTreeSet<CellRef>,
}

impl CellularSubcomplex {
    /// Admit a subcomplex.  Every declared cell must be in the ambient
    /// complex and all of its nonzero boundary faces must also be present.
    pub fn try_new(
        id: StableId,
        cells: impl IntoIterator<Item = CellRef>,
        ambient: &FiniteCellComplex,
    ) -> Result<Self, TerminalRelativeError> {
        let mut canonical = BTreeSet::new();
        for cell in cells {
            if !ambient.contains(cell) {
                return Err(TerminalRelativeError::CellOutOfRange { cell });
            }
            if !canonical.insert(cell) {
                return Err(TerminalRelativeError::DuplicateSubcomplexCell {
                    id: id.as_str().to_owned(),
                    cell,
                });
            }
        }
        for cell in &canonical {
            for incidence in ambient
                .incidences
                .iter()
                .filter(|entry| entry.upper == *cell)
            {
                if !canonical.contains(&incidence.lower) {
                    return Err(TerminalRelativeError::NotASubcomplex {
                        id: id.as_str().to_owned(),
                        cell: *cell,
                        missing_boundary: incidence.lower,
                    });
                }
            }
        }
        Ok(Self {
            id,
            cells: canonical,
        })
    }

    /// Stable subcomplex identity.
    #[must_use]
    pub const fn id(&self) -> &StableId {
        &self.id
    }

    /// Canonically ordered cells.
    #[must_use]
    pub const fn cells(&self) -> &BTreeSet<CellRef> {
        &self.cells
    }

    /// Whether the support contains no cells.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.cells.is_empty()
    }
}

/// One connected-by-declaration conductor component.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ConductorComponent {
    id: ConductorComponentId,
    support: CellularSubcomplex,
}

impl ConductorComponent {
    /// Bind a component identity to a nonempty admitted support.
    pub fn new(
        id: ConductorComponentId,
        support: CellularSubcomplex,
    ) -> Result<Self, TerminalRelativeError> {
        if support.is_empty() {
            return Err(TerminalRelativeError::EmptySupport {
                object: "conductor component",
                id: id.as_str().to_owned(),
            });
        }
        Ok(Self { id, support })
    }

    /// Component identity.
    #[must_use]
    pub const fn id(&self) -> &ConductorComponentId {
        &self.id
    }

    /// Exact cellular support.
    #[must_use]
    pub const fn support(&self) -> &CellularSubcomplex {
        &self.support
    }
}

/// Electrical role of a terminal within one phase.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum TerminalRole {
    /// Driven/source-side terminal.
    Driven,
    /// Explicit return and reference terminal.
    ReturnReference,
}

/// Physical positive-current orientation at a conductor terminal.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum TerminalOrientation {
    /// Positive current enters the conductor component.
    IntoConductor,
    /// Positive current leaves the conductor component.
    OutOfConductor,
}

/// Port coordinate selected by a physical terminal patch.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum TerminalPortCoordinate {
    /// Effort coordinate (voltage for the admitted electrical port).
    Effort,
    /// Flow coordinate (current for the admitted electrical port).
    Flow,
}

/// Sign of the explicit port-to-terminal trivialization map.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum OrientationMapSign {
    /// Port-positive and physical-terminal-positive coordinates agree.
    Preserve,
    /// The physical terminal coordinate is the negative port coordinate.
    Reverse,
}

/// Explicit voltage/current reference and orientation map for one terminal.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TerminalPortTrivialization {
    id: TrivializationId,
    port_id: StableId,
    sign: OrientationMapSign,
    voltage_reference: StableId,
    current_reference: StableId,
}

impl TerminalPortTrivialization {
    /// Declare a content-referenced terminal/port coordinate map.
    #[must_use]
    pub const fn new(
        id: TrivializationId,
        port_id: StableId,
        sign: OrientationMapSign,
        voltage_reference: StableId,
        current_reference: StableId,
    ) -> Self {
        Self {
            id,
            port_id,
            sign,
            voltage_reference,
            current_reference,
        }
    }

    /// Trivialization identity.
    #[must_use]
    pub const fn id(&self) -> &TrivializationId {
        &self.id
    }

    /// Exact port bound by this map.
    #[must_use]
    pub const fn port_id(&self) -> &StableId {
        &self.port_id
    }

    /// Orientation action.
    #[must_use]
    pub const fn sign(&self) -> OrientationMapSign {
        self.sign
    }

    /// Voltage-zero/reference artifact.
    #[must_use]
    pub const fn voltage_reference(&self) -> &StableId {
        &self.voltage_reference
    }

    /// Current-positive/reference artifact.
    #[must_use]
    pub const fn current_reference(&self) -> &StableId {
        &self.current_reference
    }
}

/// Presented dependency-neutral reference to one Machine-IR port and its
/// effort/flow terminals.
///
/// `fs-feec` deliberately does not depend on L6 `fs-ir`.  The L6 adapter must
/// validate these canonical references against an admitted MachineGraph; this
/// lower-layer schema retains a domain, version, digest, owner, port, and both
/// terminal keys without laundering them into verified authority.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PresentedMachinePortRef {
    authority_domain: StableId,
    schema_version: u32,
    graph_digest: [u8; 32],
    owner: StableId,
    port: StableId,
    effort_terminal: StableId,
    flow_terminal: StableId,
}

impl PresentedMachinePortRef {
    /// Retain presented external references.  This constructor performs local
    /// shape checks only; graph membership remains an L6 adapter obligation.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        authority_domain: StableId,
        schema_version: u32,
        graph_digest: [u8; 32],
        owner: StableId,
        port: StableId,
        effort_terminal: StableId,
        flow_terminal: StableId,
    ) -> Result<Self, TerminalRelativeError> {
        if authority_domain.as_str() != PRESENTED_MACHINE_GRAPH_DOMAIN
            || schema_version != PRESENTED_MACHINE_GRAPH_SCHEMA_VERSION
        {
            return Err(TerminalRelativeError::MachineGraphSchemaMismatch {
                expected_domain: PRESENTED_MACHINE_GRAPH_DOMAIN,
                expected_version: PRESENTED_MACHINE_GRAPH_SCHEMA_VERSION,
                actual_domain: authority_domain.as_str().to_owned(),
                actual_version: schema_version,
            });
        }
        if graph_digest == [0; 32] {
            return Err(TerminalRelativeError::ZeroMachineGraphDigest);
        }
        if effort_terminal == flow_terminal {
            return Err(TerminalRelativeError::MachinePortTerminalAlias {
                terminal: effort_terminal.as_str().to_owned(),
            });
        }
        Ok(Self {
            authority_domain,
            schema_version,
            graph_digest,
            owner,
            port,
            effort_terminal,
            flow_terminal,
        })
    }

    /// External authority domain; still only presented here.
    #[must_use]
    pub const fn authority_domain(&self) -> &StableId {
        &self.authority_domain
    }

    /// External schema version.
    #[must_use]
    pub const fn schema_version(&self) -> u32 {
        self.schema_version
    }

    /// Exact presented MachineGraph digest bytes.
    #[must_use]
    pub const fn graph_digest(&self) -> &[u8; 32] {
        &self.graph_digest
    }

    /// Owning Machine-IR subsystem key.
    #[must_use]
    pub const fn owner(&self) -> &StableId {
        &self.owner
    }

    /// Machine-IR coupling-port key.
    #[must_use]
    pub const fn port(&self) -> &StableId {
        &self.port
    }

    /// Machine-IR effort terminal key.
    #[must_use]
    pub const fn effort_terminal(&self) -> &StableId {
        &self.effort_terminal
    }

    /// Machine-IR flow terminal key.
    #[must_use]
    pub const fn flow_terminal(&self) -> &StableId {
        &self.flow_terminal
    }
}

/// One physical terminal with topology, phase, port, and Machine-IR bindings.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PhysicalTerminal {
    id: PhysicalTerminalId,
    support: CellularSubcomplex,
    component: ConductorComponentId,
    phase: PhaseId,
    role: TerminalRole,
    orientation: TerminalOrientation,
    coordinate: TerminalPortCoordinate,
    port: PortSchema,
    machine: PresentedMachinePortRef,
    trivialization: TerminalPortTrivialization,
}

impl PhysicalTerminal {
    /// Declare one terminal.  Cross-object checks occur when the complete pair
    /// is admitted.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        id: PhysicalTerminalId,
        support: CellularSubcomplex,
        component: ConductorComponentId,
        phase: PhaseId,
        role: TerminalRole,
        orientation: TerminalOrientation,
        coordinate: TerminalPortCoordinate,
        port: PortSchema,
        machine: PresentedMachinePortRef,
        trivialization: TerminalPortTrivialization,
    ) -> Result<Self, TerminalRelativeError> {
        if support.is_empty() {
            return Err(TerminalRelativeError::EmptySupport {
                object: "physical terminal",
                id: id.as_str().to_owned(),
            });
        }
        if port.kind() != PortKind::ElectricalVoltageCurrent {
            return Err(TerminalRelativeError::TerminalRequiresElectricalPort {
                terminal: id.as_str().to_owned(),
                actual: port.kind(),
            });
        }
        if coordinate != TerminalPortCoordinate::Flow {
            return Err(TerminalRelativeError::TerminalRequiresFlowCoordinate {
                terminal: id.as_str().to_owned(),
                actual: coordinate,
            });
        }
        if trivialization.port_id() != port.id() {
            return Err(TerminalRelativeError::TrivializationPortMismatch {
                terminal: id.as_str().to_owned(),
                expected: port.id().as_str().to_owned(),
                actual: trivialization.port_id().as_str().to_owned(),
            });
        }
        if machine.port() != port.id() {
            return Err(TerminalRelativeError::MachinePortSchemaMismatch {
                terminal: id.as_str().to_owned(),
                expected: port.id().as_str().to_owned(),
                actual: machine.port().as_str().to_owned(),
            });
        }
        if port.coordinates().orientation() != PortOrientation::OutwardFromOwner {
            return Err(TerminalRelativeError::UnsupportedPortOrientation {
                terminal: id.as_str().to_owned(),
                actual: port.coordinates().orientation(),
            });
        }
        let expected_orientation = match trivialization.sign() {
            OrientationMapSign::Preserve => TerminalOrientation::OutOfConductor,
            OrientationMapSign::Reverse => TerminalOrientation::IntoConductor,
        };
        if orientation != expected_orientation {
            return Err(TerminalRelativeError::TerminalOrientationMismatch {
                terminal: id.as_str().to_owned(),
                port_orientation: port.coordinates().orientation(),
                trivialization: trivialization.sign(),
                actual: orientation,
            });
        }
        Ok(Self {
            id,
            support,
            component,
            phase,
            role,
            orientation,
            coordinate,
            port,
            machine,
            trivialization,
        })
    }

    /// Physical terminal identity.
    #[must_use]
    pub const fn id(&self) -> &PhysicalTerminalId {
        &self.id
    }

    /// Exact cellular support.
    #[must_use]
    pub const fn support(&self) -> &CellularSubcomplex {
        &self.support
    }

    /// Owning conductor component.
    #[must_use]
    pub const fn component(&self) -> &ConductorComponentId {
        &self.component
    }

    /// Electrical phase identity.
    #[must_use]
    pub const fn phase(&self) -> &PhaseId {
        &self.phase
    }

    /// Driven or return/reference role.
    #[must_use]
    pub const fn role(&self) -> TerminalRole {
        self.role
    }

    /// Positive-current orientation.
    #[must_use]
    pub const fn orientation(&self) -> TerminalOrientation {
        self.orientation
    }

    /// Explicitly selected electrical flow/current coordinate.
    #[must_use]
    pub const fn coordinate(&self) -> TerminalPortCoordinate {
        self.coordinate
    }

    /// Complete neutral `PortSchema` declaration.
    #[must_use]
    pub const fn port(&self) -> &PortSchema {
        &self.port
    }

    /// Exact Machine-IR graph/terminal reference pair.
    #[must_use]
    pub const fn machine(&self) -> &PresentedMachinePortRef {
        &self.machine
    }

    /// Explicit terminal/port coordinate map.
    #[must_use]
    pub const fn trivialization(&self) -> &TerminalPortTrivialization {
        &self.trivialization
    }
}

fn phase_convention_mismatch(
    left: &PhysicalTerminal,
    right: &PhysicalTerminal,
) -> Option<&'static str> {
    let left_port = &left.port;
    let right_port = &right.port;
    if left_port.version() != right_port.version() {
        Some("port schema version")
    } else if left_port.kind() != right_port.kind() {
        Some("port kind")
    } else if left_port.effort_dimensions() != right_port.effort_dimensions() {
        Some("effort dimensions")
    } else if left_port.flow_dimensions() != right_port.flow_dimensions() {
        Some("flow dimensions")
    } else if left_port.shape() != right_port.shape() {
        Some("port value shape")
    } else if left_port.coordinates().basis() != right_port.coordinates().basis() {
        Some("coordinate basis")
    } else if left_port.coordinates().frame() != right_port.coordinates().frame() {
        Some("coordinate frame")
    } else if left_port.coordinates().orientation() != right_port.coordinates().orientation() {
        Some("coordinate orientation")
    } else if left_port.power_pairing() != right_port.power_pairing() {
        Some("power pairing")
    } else if left_port.timestamp().clock() != right_port.timestamp().clock() {
        Some("clock domain")
    } else if left_port.timestamp().tick() != right_port.timestamp().tick() {
        Some("clock tick")
    } else if left_port.conservation_roles() != right_port.conservation_roles() {
        Some("conservation roles")
    } else if left.trivialization.voltage_reference() != right.trivialization.voltage_reference() {
        Some("voltage reference")
    } else {
        None
    }
}

/// Coefficient sectors explicitly represented by the I13.2a schema.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum TerminalRelativeCoefficientDomain {
    /// Exact integral topology representative coefficients.
    Integers,
    /// Finite real physical amplitude/cochain coefficients.
    FiniteReal,
}

/// Exact incidence vocabulary admitted by this first physical schema lane.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum TerminalRelativeIncidenceDomain {
    /// Oriented chain-complex incidence with coefficients `-1` or `+1`.
    OrientedSignedUnit,
}

/// Authority status of the dependency-neutral Machine-IR references.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum MachineBindingStatus {
    /// Shape-checked and identity-bearing, but not L6 graph-membership proof.
    PresentedOnly,
}

/// Explicit no-claim states retained by a terminal-relative receipt.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum TerminalRelativeNoClaim {
    /// Relative homology classes have not been computed or quotient-identified.
    HomologyClass,
    /// Torsion and normal forms have not been computed.
    TorsionAndNormalForms,
    /// MachineGraph membership has not been verified in L6.
    MachineGraphMembership,
    /// Declared conversion artifacts have not been executed or verified.
    ConversionExecution,
    /// Geometry and manufacturing realizability have not been certified.
    GeometricRealizability,
    /// Field transfer and electromagnetic solve semantics are out of scope.
    PhysicalFieldTransfer,
    /// Relabeling/permutation/refinement equivalence requires explicit maps.
    RepresentationNaturality,
}

const TERMINAL_RELATIVE_NO_CLAIMS: [TerminalRelativeNoClaim; 7] = [
    TerminalRelativeNoClaim::HomologyClass,
    TerminalRelativeNoClaim::TorsionAndNormalForms,
    TerminalRelativeNoClaim::MachineGraphMembership,
    TerminalRelativeNoClaim::ConversionExecution,
    TerminalRelativeNoClaim::GeometricRealizability,
    TerminalRelativeNoClaim::PhysicalFieldTransfer,
    TerminalRelativeNoClaim::RepresentationNaturality,
];

/// Audit row for one terminal/PortSchema/Machine-IR binding.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TerminalBindingReceipt {
    terminal: PhysicalTerminalId,
    phase: PhaseId,
    role: TerminalRole,
    orientation: TerminalOrientation,
    coordinate: TerminalPortCoordinate,
    port_schema: PortSchema,
    machine: PresentedMachinePortRef,
    trivialization: TerminalPortTrivialization,
}

impl TerminalBindingReceipt {
    fn from_terminal(terminal: &PhysicalTerminal) -> Self {
        Self {
            terminal: terminal.id.clone(),
            phase: terminal.phase.clone(),
            role: terminal.role,
            orientation: terminal.orientation,
            coordinate: terminal.coordinate,
            port_schema: terminal.port.clone(),
            machine: terminal.machine.clone(),
            trivialization: terminal.trivialization.clone(),
        }
    }

    /// Physical terminal identity.
    #[must_use]
    pub const fn terminal(&self) -> &PhysicalTerminalId {
        &self.terminal
    }

    /// Electrical phase identity.
    #[must_use]
    pub const fn phase(&self) -> &PhaseId {
        &self.phase
    }

    /// Driven/return role.
    #[must_use]
    pub const fn role(&self) -> TerminalRole {
        self.role
    }

    /// Physical positive-current direction.
    #[must_use]
    pub const fn orientation(&self) -> TerminalOrientation {
        self.orientation
    }

    /// Selected effort/flow coordinate.
    #[must_use]
    pub const fn coordinate(&self) -> TerminalPortCoordinate {
        self.coordinate
    }

    /// Complete admitted PortSchema value.
    #[must_use]
    pub const fn port_schema(&self) -> &PortSchema {
        &self.port_schema
    }

    /// Complete presented Machine-IR binding, including domain and version.
    #[must_use]
    pub const fn machine(&self) -> &PresentedMachinePortRef {
        &self.machine
    }

    /// Complete terminal/port reference and orientation map.
    #[must_use]
    pub const fn trivialization(&self) -> &TerminalPortTrivialization {
        &self.trivialization
    }
}

/// I13.2a structural receipt for one admitted terminal-relative complex.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TerminalRelativeComplexReceipt {
    identity: IdentityReceipt<TerminalRelativePairId>,
    schema_version: u32,
    incidence_domain: TerminalRelativeIncidenceDomain,
    coefficient_domains: [TerminalRelativeCoefficientDomain; 2],
    current_dimensions: Dims,
    terminal_bindings: Vec<TerminalBindingReceipt>,
    conversion_boundaries: [DeclaredPhysicalMapKind; 2],
    machine_binding: MachineBindingStatus,
    no_claims: [TerminalRelativeNoClaim; 7],
}

impl TerminalRelativeComplexReceipt {
    fn new(
        identity: IdentityReceipt<TerminalRelativePairId>,
        terminals: &[PhysicalTerminal],
    ) -> Self {
        Self {
            identity,
            schema_version: TERMINAL_RELATIVE_SCHEMA_VERSION,
            incidence_domain: TerminalRelativeIncidenceDomain::OrientedSignedUnit,
            coefficient_domains: [
                TerminalRelativeCoefficientDomain::Integers,
                TerminalRelativeCoefficientDomain::FiniteReal,
            ],
            current_dimensions: Current::DIMS,
            terminal_bindings: terminals
                .iter()
                .map(TerminalBindingReceipt::from_terminal)
                .collect(),
            conversion_boundaries: [
                DeclaredPhysicalMapKind::WindingRealization,
                DeclaredPhysicalMapKind::CurrentRealization,
            ],
            machine_binding: MachineBindingStatus::PresentedOnly,
            no_claims: TERMINAL_RELATIVE_NO_CLAIMS,
        }
    }

    /// Strong identity plus canonical-preimage/schema audit roots.
    #[must_use]
    pub const fn identity_receipt(&self) -> IdentityReceipt<TerminalRelativePairId> {
        self.identity
    }

    /// Terminal-relative schema version.
    #[must_use]
    pub const fn schema_version(&self) -> u32 {
        self.schema_version
    }

    /// Exact admitted incidence coefficient vocabulary.
    #[must_use]
    pub const fn incidence_domain(&self) -> TerminalRelativeIncidenceDomain {
        self.incidence_domain
    }

    /// Nominally separate integral and real coefficient sectors.
    #[must_use]
    pub const fn coefficient_domains(&self) -> &[TerminalRelativeCoefficientDomain; 2] {
        &self.coefficient_domains
    }

    /// Physical dimensions required of distributed current coefficients.
    #[must_use]
    pub const fn current_dimensions(&self) -> Dims {
        self.current_dimensions
    }

    /// Canonically ordered terminal binding audit rows.
    #[must_use]
    pub fn terminal_bindings(&self) -> &[TerminalBindingReceipt] {
        &self.terminal_bindings
    }

    /// Only cross-sector map families admitted by this slice.
    #[must_use]
    pub const fn conversion_boundaries(&self) -> &[DeclaredPhysicalMapKind; 2] {
        &self.conversion_boundaries
    }

    /// Presented/verified authority state of Machine-IR references.
    #[must_use]
    pub const fn machine_binding(&self) -> MachineBindingStatus {
        self.machine_binding
    }

    /// Explicit claims this structural receipt does not make.
    #[must_use]
    pub const fn no_claims(&self) -> &[TerminalRelativeNoClaim; 7] {
        &self.no_claims
    }
}

/// Structurally admitted terminal-relative conductor pair.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TerminalRelativePair {
    complex: FiniteCellComplex,
    conductor: CellularSubcomplex,
    relative: CellularSubcomplex,
    insulation: CellularSubcomplex,
    components: Vec<ConductorComponent>,
    terminals: Vec<PhysicalTerminal>,
    phase_components: BTreeMap<PhaseId, ConductorComponentId>,
    receipt: TerminalRelativeComplexReceipt,
}

impl TerminalRelativePair {
    /// Admit the complete pair.  Components must partition the conductor;
    /// terminal supports must be disjoint subsets of their named components;
    /// terminal and insulation supports must be disjoint; each phase must
    /// declare a driven terminal and an explicit return/reference terminal;
    /// every component must be owned by exactly one phase in this first slice.
    #[allow(clippy::too_many_lines)]
    pub fn try_new(
        complex: FiniteCellComplex,
        conductor: CellularSubcomplex,
        relative: CellularSubcomplex,
        insulation: CellularSubcomplex,
        mut components: Vec<ConductorComponent>,
        mut terminals: Vec<PhysicalTerminal>,
    ) -> Result<Self, TerminalRelativeError> {
        validate_subcomplex_against(&conductor, &complex)?;
        validate_subcomplex_against(&relative, &complex)?;
        validate_subcomplex_against(&insulation, &complex)?;
        if conductor.is_empty() {
            return Err(TerminalRelativeError::EmptySupport {
                object: "conductor",
                id: conductor.id().as_str().to_owned(),
            });
        }
        if let Some(cell) = relative
            .cells
            .iter()
            .copied()
            .find(|cell| !conductor.cells.contains(cell))
        {
            return Err(TerminalRelativeError::RelativeOutsideConductor { cell });
        }
        if let Some(cell) = insulation
            .cells
            .iter()
            .copied()
            .find(|cell| !relative.cells.contains(cell))
        {
            return Err(TerminalRelativeError::InsulationOutsideRelativeSubcomplex { cell });
        }
        if components.is_empty() || components.len() > MAX_CONDUCTOR_COMPONENTS {
            return Err(TerminalRelativeError::ComponentBudgetExceeded {
                actual: components.len(),
                max: MAX_CONDUCTOR_COMPONENTS,
            });
        }
        if terminals.len() < 2 || terminals.len() > MAX_PHYSICAL_TERMINALS {
            return Err(TerminalRelativeError::TerminalBudgetExceeded {
                actual: terminals.len(),
                max: MAX_PHYSICAL_TERMINALS,
            });
        }
        components.sort_by(|left, right| left.id.cmp(&right.id));
        terminals.sort_by(|left, right| left.id.cmp(&right.id));

        let mut boundary_faces = BTreeMap::<CellRef, Vec<CellRef>>::new();
        for incidence in &complex.incidences {
            boundary_faces
                .entry(incidence.upper)
                .or_default()
                .push(incidence.lower);
        }

        let mut component_ids = BTreeSet::new();
        let mut covered = BTreeSet::new();
        let mut component_union = BTreeSet::new();
        for component in &components {
            validate_subcomplex_against(&component.support, &complex)?;
            if !component_ids.insert(component.id.clone()) {
                return Err(TerminalRelativeError::DuplicateIdentity {
                    role: "conductor component",
                    id: component.id.as_str().to_owned(),
                });
            }
            let top_cells: Vec<_> = component
                .support
                .cells
                .iter()
                .copied()
                .filter(|cell| cell.degree == complex.dimension)
                .collect();
            if top_cells.is_empty() {
                return Err(TerminalRelativeError::ComponentHasNoTopCell {
                    component: component.id.as_str().to_owned(),
                    top_degree: complex.dimension,
                });
            }
            let expected_support = downward_closure(top_cells, &boundary_faces);
            if expected_support != component.support.cells {
                let cell = expected_support
                    .symmetric_difference(&component.support.cells)
                    .next()
                    .copied()
                    .expect("unequal finite sets have a witness");
                return Err(TerminalRelativeError::ComponentSupportNotTopClosure {
                    component: component.id.as_str().to_owned(),
                    cell,
                });
            }
            for cell in component.support.cells() {
                if !conductor.cells.contains(cell) {
                    return Err(TerminalRelativeError::ComponentOutsideConductor {
                        component: component.id.as_str().to_owned(),
                        cell: *cell,
                    });
                }
                component_union.insert(*cell);
                if cell.degree == complex.dimension && !covered.insert(*cell) {
                    return Err(TerminalRelativeError::OverlappingComponents { cell: *cell });
                }
            }
        }
        if component_union != conductor.cells {
            let missing = conductor
                .cells
                .difference(&component_union)
                .next()
                .copied()
                .or_else(|| component_union.difference(&conductor.cells).next().copied())
                .expect("unequal finite sets have a witness");
            return Err(TerminalRelativeError::ComponentPartitionMismatch { cell: missing });
        }
        let conductor_top_cells: BTreeSet<_> = conductor
            .cells
            .iter()
            .copied()
            .filter(|cell| cell.degree == complex.dimension)
            .collect();
        if covered != conductor_top_cells {
            let missing = conductor_top_cells
                .difference(&covered)
                .next()
                .copied()
                .or_else(|| covered.difference(&conductor_top_cells).next().copied())
                .expect("unequal finite sets have a witness");
            return Err(TerminalRelativeError::ComponentPartitionMismatch { cell: missing });
        }

        let component_map: BTreeMap<_, _> = components
            .iter()
            .map(|component| (&component.id, component))
            .collect();
        let mut terminal_ids = BTreeSet::new();
        let mut port_ids = BTreeSet::new();
        let mut machine_terminals = BTreeSet::new();
        let mut trivialization_ids = BTreeSet::new();
        let mut terminal_cells = BTreeSet::new();
        let mut phases = BTreeMap::<PhaseId, Vec<&PhysicalTerminal>>::new();

        for terminal in &terminals {
            validate_subcomplex_against(&terminal.support, &complex)?;
            if !terminal_ids.insert(terminal.id.clone()) {
                return Err(TerminalRelativeError::DuplicateIdentity {
                    role: "physical terminal",
                    id: terminal.id.as_str().to_owned(),
                });
            }
            if !port_ids.insert(terminal.port.id().as_str().to_owned()) {
                return Err(TerminalRelativeError::DuplicateIdentity {
                    role: "port",
                    id: terminal.port.id().as_str().to_owned(),
                });
            }
            let machine_key = (
                *terminal.machine.graph_digest(),
                terminal.machine.port().as_str().to_owned(),
            );
            if !machine_terminals.insert(machine_key) {
                return Err(TerminalRelativeError::DuplicateIdentity {
                    role: "presented Machine-IR port binding",
                    id: terminal.machine.port().as_str().to_owned(),
                });
            }
            if !trivialization_ids.insert(terminal.trivialization.id.clone()) {
                return Err(TerminalRelativeError::DuplicateIdentity {
                    role: "terminal trivialization",
                    id: terminal.trivialization.id.as_str().to_owned(),
                });
            }
            let Some(component) = component_map.get(&terminal.component) else {
                return Err(TerminalRelativeError::UnknownTerminalComponent {
                    terminal: terminal.id.as_str().to_owned(),
                    component: terminal.component.as_str().to_owned(),
                });
            };
            let Some(terminal_degree) = complex.dimension.checked_sub(1) else {
                return Err(TerminalRelativeError::TerminalCodimension {
                    terminal: terminal.id.as_str().to_owned(),
                    ambient_dimension: complex.dimension,
                });
            };
            if !terminal
                .support
                .cells
                .iter()
                .any(|cell| cell.degree == terminal_degree)
                || terminal
                    .support
                    .cells
                    .iter()
                    .any(|cell| cell.degree > terminal_degree)
            {
                return Err(TerminalRelativeError::TerminalCodimension {
                    terminal: terminal.id.as_str().to_owned(),
                    ambient_dimension: complex.dimension,
                });
            }
            let patch_cells: Vec<_> = terminal
                .support
                .cells
                .iter()
                .copied()
                .filter(|cell| cell.degree == terminal_degree)
                .collect();
            let expected_support = downward_closure(patch_cells, &boundary_faces);
            if expected_support != terminal.support.cells {
                let cell = expected_support
                    .symmetric_difference(&terminal.support.cells)
                    .next()
                    .copied()
                    .expect("unequal finite sets have a witness");
                return Err(TerminalRelativeError::TerminalSupportNotPatchClosure {
                    terminal: terminal.id.as_str().to_owned(),
                    cell,
                });
            }
            for cell in terminal.support.cells() {
                if !component.support.cells.contains(cell) {
                    return Err(TerminalRelativeError::TerminalOutsideComponent {
                        terminal: terminal.id.as_str().to_owned(),
                        component: component.id.as_str().to_owned(),
                        cell: *cell,
                    });
                }
                if insulation.cells.contains(cell) {
                    return Err(TerminalRelativeError::TerminalInsulationOverlap {
                        terminal: terminal.id.as_str().to_owned(),
                        cell: *cell,
                    });
                }
                if !relative.cells.contains(cell) {
                    return Err(TerminalRelativeError::TerminalOutsideRelativeSubcomplex {
                        terminal: terminal.id.as_str().to_owned(),
                        cell: *cell,
                    });
                }
                if cell.degree == terminal_degree {
                    let incident_component_cells = complex
                        .incidences
                        .iter()
                        .filter(|incidence| {
                            incidence.lower == *cell
                                && incidence.upper.degree == complex.dimension
                                && component.support.cells.contains(&incidence.upper)
                        })
                        .count();
                    if incident_component_cells != 1 {
                        return Err(TerminalRelativeError::TerminalNotOnComponentBoundary {
                            terminal: terminal.id.as_str().to_owned(),
                            cell: *cell,
                        });
                    }
                }
                if !terminal_cells.insert(*cell) {
                    return Err(TerminalRelativeError::OverlappingTerminals { cell: *cell });
                }
            }
            phases
                .entry(terminal.phase.clone())
                .or_default()
                .push(terminal);
        }
        let mut phase_components = BTreeMap::new();
        let mut component_phases = BTreeMap::<ConductorComponentId, PhaseId>::new();
        for (phase, phase_terminals) in &phases {
            if phase_terminals.len() != 2 {
                return Err(TerminalRelativeError::PhaseTerminalCardinality {
                    phase: phase.as_str().to_owned(),
                    actual: phase_terminals.len(),
                });
            }
            let driven = phase_terminals
                .iter()
                .filter(|terminal| terminal.role == TerminalRole::Driven)
                .count();
            let return_reference = phase_terminals
                .iter()
                .filter(|terminal| terminal.role == TerminalRole::ReturnReference)
                .count();
            if driven == 0 {
                return Err(TerminalRelativeError::MissingPhaseRole {
                    phase: phase.as_str().to_owned(),
                    role: TerminalRole::Driven,
                });
            }
            if return_reference == 0 {
                return Err(TerminalRelativeError::MissingPhaseRole {
                    phase: phase.as_str().to_owned(),
                    role: TerminalRole::ReturnReference,
                });
            }
            if driven != 1 || return_reference != 1 {
                return Err(TerminalRelativeError::DuplicatePhaseRole {
                    phase: phase.as_str().to_owned(),
                });
            }
            let driven_terminal = phase_terminals
                .iter()
                .find(|terminal| terminal.role == TerminalRole::Driven)
                .ok_or_else(|| TerminalRelativeError::MissingPhaseRole {
                    phase: phase.as_str().to_owned(),
                    role: TerminalRole::Driven,
                })?;
            let return_terminal = phase_terminals
                .iter()
                .find(|terminal| terminal.role == TerminalRole::ReturnReference)
                .ok_or_else(|| TerminalRelativeError::MissingPhaseRole {
                    phase: phase.as_str().to_owned(),
                    role: TerminalRole::ReturnReference,
                })?;
            let component = &driven_terminal.component;
            if return_terminal.component != *component {
                return Err(TerminalRelativeError::PhaseComponentMismatch {
                    phase: phase.as_str().to_owned(),
                    driven_component: component.as_str().to_owned(),
                    return_component: return_terminal.component.as_str().to_owned(),
                });
            }
            if let Some(existing) = component_phases.insert(component.clone(), phase.clone()) {
                return Err(TerminalRelativeError::ComponentPhaseConflict {
                    component: component.as_str().to_owned(),
                    first_phase: existing.as_str().to_owned(),
                    second_phase: phase.as_str().to_owned(),
                });
            }
            phase_components.insert(phase.clone(), component.clone());
            let into = phase_terminals
                .iter()
                .filter(|terminal| terminal.orientation == TerminalOrientation::IntoConductor)
                .count();
            let out = phase_terminals
                .iter()
                .filter(|terminal| terminal.orientation == TerminalOrientation::OutOfConductor)
                .count();
            if into != 1 || out != 1 {
                return Err(TerminalRelativeError::PhaseOrientationDoesNotClose {
                    phase: phase.as_str().to_owned(),
                });
            }
            if let Some(field) = phase_convention_mismatch(phase_terminals[0], phase_terminals[1]) {
                return Err(TerminalRelativeError::PhaseConventionMismatch {
                    phase: phase.as_str().to_owned(),
                    field,
                });
            }
        }
        for component in &components {
            if !component_phases.contains_key(&component.id) {
                return Err(TerminalRelativeError::UnboundConductorComponent {
                    component: component.id.as_str().to_owned(),
                });
            }
        }
        drop(phases);

        let payload = canonical_pair_payload(
            &complex,
            &conductor,
            &relative,
            &insulation,
            &components,
            &terminals,
        )?;
        let receipt =
            CanonicalEncoder::<TerminalRelativePairId, _>::new(PAIR_IDENTITY_LIMITS, NeverCancel)?
                .bytes(Field::new(0, "terminal-relative-payload"), &payload)?
                .finish()?;
        let receipt = TerminalRelativeComplexReceipt::new(receipt, &terminals);
        Ok(Self {
            complex,
            conductor,
            relative,
            insulation,
            components,
            terminals,
            phase_components,
            receipt,
        })
    }

    /// Strong semantic identity of the complete pair.
    #[must_use]
    pub const fn identity(&self) -> TerminalRelativePairId {
        self.receipt.identity.id()
    }

    /// Exact typed canonical-frame byte count absorbed for the identity.
    #[must_use]
    pub const fn canonical_bytes(&self) -> u64 {
        self.receipt.identity.canonical_bytes()
    }

    /// Complete I13.2a structural and no-claim receipt.
    #[must_use]
    pub const fn complex_receipt(&self) -> &TerminalRelativeComplexReceipt {
        &self.receipt
    }

    /// Exact admitted incidence complex.
    #[must_use]
    pub const fn complex(&self) -> &FiniteCellComplex {
        &self.complex
    }

    /// Conductor subcomplex.
    #[must_use]
    pub const fn conductor(&self) -> &CellularSubcomplex {
        &self.conductor
    }

    /// Explicit quotient subcomplex `A`; it is never inferred from material
    /// or insulation labels.
    #[must_use]
    pub const fn relative(&self) -> &CellularSubcomplex {
        &self.relative
    }

    /// Insulation subcomplex.
    #[must_use]
    pub const fn insulation(&self) -> &CellularSubcomplex {
        &self.insulation
    }

    /// Canonically ordered conductor components.
    #[must_use]
    pub fn components(&self) -> &[ConductorComponent] {
        &self.components
    }

    /// Canonically ordered physical terminals.
    #[must_use]
    pub fn terminals(&self) -> &[PhysicalTerminal] {
        &self.terminals
    }

    /// Conductor component explicitly owned by a phase's two terminal rows.
    #[must_use]
    pub fn phase_component(&self, phase: &PhaseId) -> Option<&ConductorComponentId> {
        self.phase_components.get(phase)
    }

    /// Canonical quotient basis: conductor cells not killed by the explicitly
    /// declared relative subcomplex.
    #[must_use]
    pub fn relative_basis(&self, degree: u8) -> Vec<CellRef> {
        self.conductor
            .cells
            .iter()
            .copied()
            .filter(|cell| cell.degree == degree && !self.relative.cells.contains(cell))
            .collect()
    }

    /// Canonical quotient basis restricted to the component owned by `phase`.
    pub fn phase_relative_basis(
        &self,
        phase: &PhaseId,
        degree: u8,
    ) -> Result<Vec<CellRef>, TerminalRelativeError> {
        let Some(component_id) = self.phase_components.get(phase) else {
            return Err(TerminalRelativeError::UnknownPhase {
                phase: phase.as_str().to_owned(),
            });
        };
        let Some(component) = self
            .components
            .iter()
            .find(|component| &component.id == component_id)
        else {
            return Err(TerminalRelativeError::PhaseComponentBindingLost {
                phase: phase.as_str().to_owned(),
                component: component_id.as_str().to_owned(),
            });
        };
        Ok(component
            .support
            .cells
            .iter()
            .copied()
            .filter(|cell| cell.degree == degree && !self.relative.cells.contains(cell))
            .collect())
    }

    /// Whether a phase has an admitted component binding.
    #[must_use]
    pub fn contains_phase(&self, phase: &PhaseId) -> bool {
        self.phase_components.contains_key(phase)
    }

    /// Apply the exact relative boundary map to an integral chain.
    pub fn boundary(
        &self,
        chain: &IntegralRelativeChain,
    ) -> Result<IntegralRelativeChain, TerminalRelativeError> {
        if chain.pair != self.identity() {
            return Err(TerminalRelativeError::PairIdentityMismatch);
        }
        let Some(target_degree) = chain.degree.checked_sub(1) else {
            return Err(TerminalRelativeError::NoBoundaryPredecessor);
        };
        let source_basis = self.phase_relative_basis(&chain.phase, chain.degree)?;
        let target_basis = self.phase_relative_basis(&chain.phase, target_degree)?;
        if source_basis.len() != chain.coefficients.len() {
            return Err(TerminalRelativeError::CoefficientArity {
                expected: source_basis.len(),
                actual: chain.coefficients.len(),
            });
        }
        let target_indices: BTreeMap<_, _> = target_basis
            .iter()
            .copied()
            .enumerate()
            .map(|(index, cell)| (cell, index))
            .collect();
        let mut accumulated = vec![0_i128; target_basis.len()];
        for (source, coefficient) in source_basis.iter().zip(&chain.coefficients) {
            for incidence in self
                .complex
                .incidences
                .iter()
                .filter(|entry| entry.upper == *source)
            {
                if let Some(target) = target_indices.get(&incidence.lower) {
                    accumulated[*target] += i128::from(*coefficient) * incidence.sign.as_i128();
                }
            }
        }
        let coefficients = accumulated
            .into_iter()
            .map(|value| {
                i64::try_from(value).map_err(|_| TerminalRelativeError::CoefficientOverflow)
            })
            .collect::<Result<Vec<_>, _>>()?;
        IntegralRelativeChain::try_new(self, chain.phase.clone(), target_degree, coefficients)
    }

    /// Apply the transpose incidence as a real relative coboundary map.
    pub fn coboundary(
        &self,
        cochain: &RealRelativeCochain,
    ) -> Result<RealRelativeCochain, TerminalRelativeError> {
        if cochain.pair != self.identity() {
            return Err(TerminalRelativeError::PairIdentityMismatch);
        }
        let Some(target_degree) = cochain.degree.checked_add(1) else {
            return Err(TerminalRelativeError::NoCoboundarySuccessor);
        };
        if target_degree > self.complex.dimension {
            return Err(TerminalRelativeError::NoCoboundarySuccessor);
        }
        let source_basis = self.phase_relative_basis(&cochain.phase, cochain.degree)?;
        let target_basis = self.phase_relative_basis(&cochain.phase, target_degree)?;
        if source_basis.len() != cochain.values.len() {
            return Err(TerminalRelativeError::CoefficientArity {
                expected: source_basis.len(),
                actual: cochain.values.len(),
            });
        }
        let source_indices: BTreeMap<_, _> = source_basis
            .iter()
            .copied()
            .enumerate()
            .map(|(index, cell)| (cell, index))
            .collect();
        let mut values = vec![0.0_f64; target_basis.len()];
        for (target_index, target) in target_basis.iter().enumerate() {
            for incidence in self
                .complex
                .incidences
                .iter()
                .filter(|entry| entry.upper == *target)
            {
                if let Some(source_index) = source_indices.get(&incidence.lower) {
                    let sign = match incidence.sign {
                        IncidenceSign::Negative => -1.0,
                        IncidenceSign::Positive => 1.0,
                    };
                    values[target_index] += sign * cochain.values[*source_index];
                }
            }
            if !values[target_index].is_finite() {
                return Err(TerminalRelativeError::NonFiniteRealCoefficient {
                    index: target_index,
                });
            }
        }
        RealRelativeCochain::try_new(
            self,
            cochain.phase.clone(),
            target_degree,
            cochain.units,
            values,
        )
    }

    /// Apply the transpose incidence as an exact integral relative
    /// coboundary map. No real-to-integer coercion is available.
    pub fn integral_coboundary(
        &self,
        cochain: &IntegralRelativeCochain,
    ) -> Result<IntegralRelativeCochain, TerminalRelativeError> {
        if cochain.pair != self.identity() {
            return Err(TerminalRelativeError::PairIdentityMismatch);
        }
        let Some(target_degree) = cochain.degree.checked_add(1) else {
            return Err(TerminalRelativeError::NoCoboundarySuccessor);
        };
        if target_degree > self.complex.dimension {
            return Err(TerminalRelativeError::NoCoboundarySuccessor);
        }
        let source_basis = self.phase_relative_basis(&cochain.phase, cochain.degree)?;
        let target_basis = self.phase_relative_basis(&cochain.phase, target_degree)?;
        if source_basis.len() != cochain.coefficients.len() {
            return Err(TerminalRelativeError::CoefficientArity {
                expected: source_basis.len(),
                actual: cochain.coefficients.len(),
            });
        }
        let source_indices: BTreeMap<_, _> = source_basis
            .iter()
            .copied()
            .enumerate()
            .map(|(index, cell)| (cell, index))
            .collect();
        let mut accumulated = vec![0_i128; target_basis.len()];
        for (target_index, target) in target_basis.iter().enumerate() {
            for incidence in self
                .complex
                .incidences
                .iter()
                .filter(|entry| entry.upper == *target)
            {
                if let Some(source_index) = source_indices.get(&incidence.lower) {
                    let contribution =
                        i128::from(cochain.coefficients[*source_index]) * incidence.sign.as_i128();
                    accumulated[target_index] = accumulated[target_index]
                        .checked_add(contribution)
                        .ok_or(TerminalRelativeError::CoefficientOverflow)?;
                }
            }
        }
        let coefficients = accumulated
            .into_iter()
            .map(|value| {
                i64::try_from(value).map_err(|_| TerminalRelativeError::CoefficientOverflow)
            })
            .collect::<Result<Vec<_>, _>>()?;
        IntegralRelativeCochain::try_new(self, cochain.phase.clone(), target_degree, coefficients)
    }

    /// Exact evaluation pairing between an integral cochain and chain.
    pub fn integral_pairing(
        &self,
        cochain: &IntegralRelativeCochain,
        chain: &IntegralRelativeChain,
    ) -> Result<i128, TerminalRelativeError> {
        if cochain.pair != self.identity() || chain.pair != self.identity() {
            return Err(TerminalRelativeError::PairIdentityMismatch);
        }
        if cochain.phase != chain.phase {
            return Err(TerminalRelativeError::PairingPhaseMismatch {
                cochain: cochain.phase.as_str().to_owned(),
                chain: chain.phase.as_str().to_owned(),
            });
        }
        if cochain.degree != chain.degree {
            return Err(TerminalRelativeError::PairingDegreeMismatch {
                cochain: cochain.degree,
                chain: chain.degree,
            });
        }
        let expected = self.phase_relative_basis(&chain.phase, chain.degree)?.len();
        if cochain.coefficients.len() != expected {
            return Err(TerminalRelativeError::CoefficientArity {
                expected,
                actual: cochain.coefficients.len(),
            });
        }
        if chain.coefficients.len() != expected {
            return Err(TerminalRelativeError::CoefficientArity {
                expected,
                actual: chain.coefficients.len(),
            });
        }
        cochain
            .coefficients
            .iter()
            .zip(&chain.coefficients)
            .try_fold(0_i128, |sum, (dual, primal)| {
                let product = i128::from(*dual) * i128::from(*primal);
                sum.checked_add(product)
                    .ok_or(TerminalRelativeError::PairingOverflow)
            })
    }
}

/// Admitted complete orientation-aware relabeling between terminal-relative
/// pairs.
///
/// The receipt proves only this explicit cell bijection and the structural
/// semantics checked by [`Self::try_new`].  It does not discharge the general
/// representation-naturality no-claim retained by either pair.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TerminalRelativeSignedRelabel {
    source_pair: TerminalRelativePairId,
    target_pair: TerminalRelativePairId,
    entries: Vec<SignedCellRelabelEntry>,
    identity_receipt: IdentityReceipt<TerminalRelativeSignedRelabelId>,
}

impl TerminalRelativeSignedRelabel {
    /// Admit a complete same-degree cell bijection which commutes exactly with
    /// signed incidence and preserves all terminal-relative semantic labels.
    ///
    /// Caller declaration order is discarded before validation and identity
    /// encoding.  No phase-current or other physical-coordinate compensation
    /// is inferred from cell-orientation signs.
    #[allow(clippy::too_many_lines)]
    pub fn try_new(
        source: &TerminalRelativePair,
        target: &TerminalRelativePair,
        mut entries: Vec<SignedCellRelabelEntry>,
    ) -> Result<Self, TerminalRelativeSignedRelabelError> {
        if source.complex.dimension != target.complex.dimension {
            return Err(
                TerminalRelativeSignedRelabelError::ComplexDimensionMismatch {
                    source: source.complex.dimension,
                    target: target.complex.dimension,
                },
            );
        }
        for (degree, (source_count, target_count)) in source
            .complex
            .cell_counts
            .iter()
            .zip(&target.complex.cell_counts)
            .enumerate()
        {
            if source_count != target_count {
                return Err(TerminalRelativeSignedRelabelError::CellCountMismatch {
                    degree: u8::try_from(degree)
                        .expect("admitted cell-complex degree always fits in u8"),
                    source: *source_count,
                    target: *target_count,
                });
            }
        }

        let expected_entries = source
            .complex
            .cell_counts
            .iter()
            .fold(0_usize, |total, count| {
                total.saturating_add(usize::try_from(*count).unwrap_or(usize::MAX))
            });
        if entries.len() != expected_entries {
            return Err(TerminalRelativeSignedRelabelError::EntryCountMismatch {
                expected: expected_entries,
                actual: entries.len(),
            });
        }

        entries.sort_unstable_by_key(|entry| entry.source);
        let mut target_cells = BTreeSet::new();
        for entry in &entries {
            if !source.complex.contains(entry.source) {
                return Err(TerminalRelativeSignedRelabelError::SourceCellOutOfRange {
                    cell: entry.source,
                });
            }
            if !target.complex.contains(entry.target) {
                return Err(TerminalRelativeSignedRelabelError::TargetCellOutOfRange {
                    cell: entry.target,
                });
            }
            if entry.source.degree != entry.target.degree {
                return Err(TerminalRelativeSignedRelabelError::CellDegreeMismatch {
                    source: entry.source,
                    target: entry.target,
                });
            }
            if !target_cells.insert(entry.target) {
                return Err(TerminalRelativeSignedRelabelError::DuplicateTargetCell {
                    cell: entry.target,
                });
            }
        }
        if let Some(duplicate) = entries
            .windows(2)
            .find(|pair| pair[0].source == pair[1].source)
            .map(|pair| pair[0].source)
        {
            return Err(TerminalRelativeSignedRelabelError::DuplicateSourceCell {
                cell: duplicate,
            });
        }

        let entry_map: BTreeMap<_, _> = entries
            .iter()
            .copied()
            .map(|entry| (entry.source, entry))
            .collect();
        verify_signed_incidence_commutation(source, target, &entry_map)?;
        verify_mapped_subcomplex(
            "conductor",
            None,
            &source.conductor,
            &target.conductor,
            &entry_map,
        )?;
        verify_mapped_subcomplex(
            "relative subcomplex",
            None,
            &source.relative,
            &target.relative,
            &entry_map,
        )?;
        verify_mapped_subcomplex(
            "insulation",
            None,
            &source.insulation,
            &target.insulation,
            &entry_map,
        )?;
        verify_component_semantics(source, target, &entry_map)?;
        verify_phase_semantics(source, target)?;
        verify_terminal_semantics(source, target, &entry_map)?;

        let payload = canonical_signed_relabel_payload(&entries)?;
        let identity_receipt = CanonicalEncoder::<TerminalRelativeSignedRelabelId, _>::new(
            SIGNED_RELABEL_IDENTITY_LIMITS,
            NeverCancel,
        )?
        .child(Field::new(0, "source-pair"), source.identity())?
        .child(Field::new(1, "target-pair"), target.identity())?
        .bytes(Field::new(2, "signed-cell-map"), &payload)?
        .finish()?;
        Ok(Self {
            source_pair: source.identity(),
            target_pair: target.identity(),
            entries,
            identity_receipt,
        })
    }

    /// Strong identity of the exact directed signed relabeling.
    #[must_use]
    pub const fn identity(&self) -> TerminalRelativeSignedRelabelId {
        self.identity_receipt.id()
    }

    /// Strong identity plus canonical-preimage/schema audit roots.
    #[must_use]
    pub const fn identity_receipt(&self) -> IdentityReceipt<TerminalRelativeSignedRelabelId> {
        self.identity_receipt
    }

    /// Strong identity of the admitted source pair.
    #[must_use]
    pub const fn source_pair_id(&self) -> TerminalRelativePairId {
        self.source_pair
    }

    /// Strong identity of the admitted target pair.
    #[must_use]
    pub const fn target_pair_id(&self) -> TerminalRelativePairId {
        self.target_pair
    }

    /// Canonically source-sorted complete signed cell map.
    #[must_use]
    pub fn entries(&self) -> &[SignedCellRelabelEntry] {
        &self.entries
    }

    /// Target cell and orientation sign for one source cell.
    #[must_use]
    pub fn image(&self, source: CellRef) -> Option<(CellRef, IncidenceSign)> {
        self.entries
            .binary_search_by_key(&source, |entry| entry.source)
            .ok()
            .map(|index| {
                let entry = self.entries[index];
                (entry.target, entry.sign)
            })
    }

    /// Push an integral chain through only the admitted signed cell
    /// reindexing.
    pub fn transport_integral_chain(
        &self,
        source: &TerminalRelativePair,
        target: &TerminalRelativePair,
        chain: &IntegralRelativeChain,
    ) -> Result<IntegralRelativeChain, TerminalRelativeSignedRelabelError> {
        self.verify_pair_bindings(source, target)?;
        if chain.pair != self.source_pair {
            return Err(TerminalRelativeSignedRelabelError::PairIdentityMismatch {
                role: "integral chain source",
                expected: self.source_pair,
                actual: chain.pair,
            });
        }
        let coefficients = self.transport_integral_coefficients(
            source,
            target,
            &chain.phase,
            chain.degree,
            &chain.coefficients,
        )?;
        IntegralRelativeChain::try_new(target, chain.phase.clone(), chain.degree, coefficients)
            .map_err(Into::into)
    }

    /// Push an integral cochain through the same signed basis reindexing.
    ///
    /// Signed permutations are orthogonal, so this action preserves raw
    /// chain/cochain evaluation pairing without introducing any phase-current
    /// compensation.
    pub fn transport_integral_cochain(
        &self,
        source: &TerminalRelativePair,
        target: &TerminalRelativePair,
        cochain: &IntegralRelativeCochain,
    ) -> Result<IntegralRelativeCochain, TerminalRelativeSignedRelabelError> {
        self.verify_pair_bindings(source, target)?;
        if cochain.pair != self.source_pair {
            return Err(TerminalRelativeSignedRelabelError::PairIdentityMismatch {
                role: "integral cochain source",
                expected: self.source_pair,
                actual: cochain.pair,
            });
        }
        let coefficients = self.transport_integral_coefficients(
            source,
            target,
            &cochain.phase,
            cochain.degree,
            &cochain.coefficients,
        )?;
        IntegralRelativeCochain::try_new(
            target,
            cochain.phase.clone(),
            cochain.degree,
            coefficients,
        )
        .map_err(Into::into)
    }

    /// Transport a winding representative with the same cell-only chain map
    /// and re-admit its exact relative-cycle receipt on the target pair.
    pub fn transport_winding_representative(
        &self,
        source: &TerminalRelativePair,
        target: &TerminalRelativePair,
        representative: &IntegralWindingRepresentative,
    ) -> Result<IntegralWindingRepresentative, TerminalRelativeSignedRelabelError> {
        let chain = self.transport_integral_chain(source, target, &representative.chain)?;
        IntegralWindingRepresentative::try_new(target, chain.phase.clone(), chain.coefficients)
            .map_err(Into::into)
    }

    /// Admit the exact inverse signed permutation.
    pub fn inverse(
        &self,
        source: &TerminalRelativePair,
        target: &TerminalRelativePair,
    ) -> Result<Self, TerminalRelativeSignedRelabelError> {
        self.verify_pair_bindings(source, target)?;
        let entries = self
            .entries
            .iter()
            .map(|entry| SignedCellRelabelEntry::new(entry.target, entry.source, entry.sign))
            .collect();
        Self::try_new(target, source, entries)
    }

    /// Compose `next` after `self`, admitting `next ∘ self` as a fresh
    /// canonical signed relabeling.
    pub fn compose(
        &self,
        next: &Self,
        source: &TerminalRelativePair,
        intermediate: &TerminalRelativePair,
        target: &TerminalRelativePair,
    ) -> Result<Self, TerminalRelativeSignedRelabelError> {
        self.verify_pair_bindings(source, intermediate)?;
        next.verify_pair_bindings(intermediate, target)?;
        let mut entries = Vec::with_capacity(self.entries.len());
        for entry in &self.entries {
            let Some((target_cell, next_sign)) = next.image(entry.target) else {
                return Err(TerminalRelativeSignedRelabelError::MissingSourceCell {
                    cell: entry.target,
                });
            };
            entries.push(SignedCellRelabelEntry::new(
                entry.source,
                target_cell,
                multiply_incidence_sign(entry.sign, next_sign),
            ));
        }
        Self::try_new(source, target, entries)
    }

    fn verify_pair_bindings(
        &self,
        source: &TerminalRelativePair,
        target: &TerminalRelativePair,
    ) -> Result<(), TerminalRelativeSignedRelabelError> {
        if source.identity() != self.source_pair {
            return Err(TerminalRelativeSignedRelabelError::PairIdentityMismatch {
                role: "source pair",
                expected: self.source_pair,
                actual: source.identity(),
            });
        }
        if target.identity() != self.target_pair {
            return Err(TerminalRelativeSignedRelabelError::PairIdentityMismatch {
                role: "target pair",
                expected: self.target_pair,
                actual: target.identity(),
            });
        }
        Ok(())
    }

    fn transport_integral_coefficients(
        &self,
        source: &TerminalRelativePair,
        target: &TerminalRelativePair,
        phase: &PhaseId,
        degree: u8,
        coefficients: &[i64],
    ) -> Result<Vec<i64>, TerminalRelativeSignedRelabelError> {
        let source_basis = source.phase_relative_basis(phase, degree)?;
        if source_basis.len() != coefficients.len() {
            return Err(TerminalRelativeError::CoefficientArity {
                expected: source_basis.len(),
                actual: coefficients.len(),
            }
            .into());
        }
        let target_basis = target.phase_relative_basis(phase, degree)?;
        let target_indices: BTreeMap<_, _> = target_basis
            .iter()
            .copied()
            .enumerate()
            .map(|(index, cell)| (cell, index))
            .collect();
        let mut transported = vec![0_i64; target_basis.len()];
        for (source_cell, coefficient) in source_basis.iter().zip(coefficients) {
            let Some((target_cell, sign)) = self.image(*source_cell) else {
                return Err(TerminalRelativeSignedRelabelError::MissingSourceCell {
                    cell: *source_cell,
                });
            };
            let Some(target_index) = target_indices.get(&target_cell) else {
                return Err(TerminalRelativeSignedRelabelError::MappedBasisCellMissing {
                    phase: phase.as_str().to_owned(),
                    degree,
                    source: *source_cell,
                    target: target_cell,
                });
            };
            transported[*target_index] = match sign {
                IncidenceSign::Positive => *coefficient,
                IncidenceSign::Negative => coefficient.checked_neg().ok_or(
                    TerminalRelativeSignedRelabelError::CoefficientOverflow { cell: *source_cell },
                )?,
            };
        }
        Ok(transported)
    }
}

/// Admitted physical relabeling with explicit cell, component, phase,
/// terminal, and phase-current actions.
///
/// This receipt proves only the declared finite bijections and the commuting
/// squares checked by [`Self::try_new`]. It is not refinement, remesh,
/// topology-event, homology-class, field-transfer, or MachineGraph authority.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TerminalRelativePhysicalRelabel {
    source_pair: TerminalRelativePairId,
    target_pair: TerminalRelativePairId,
    cells: Vec<SignedCellRelabelEntry>,
    semantics: TerminalRelativeSemanticPermutation,
    identity_receipt: IdentityReceipt<TerminalRelativePhysicalRelabelId>,
}

impl TerminalRelativePhysicalRelabel {
    /// Admit complete physical relabeling data between two terminal-relative
    /// pairs.
    ///
    /// Every map is supplied explicitly. Declaration order is discarded only
    /// after total source/target bijection checks; no semantic row or current
    /// action is inferred from canonical ordering or terminal metadata.
    #[allow(clippy::too_many_lines)]
    pub fn try_new(
        source: &TerminalRelativePair,
        target: &TerminalRelativePair,
        cells: Vec<SignedCellRelabelEntry>,
        semantics: TerminalRelativeSemanticPermutation,
    ) -> Result<Self, TerminalRelativePhysicalRelabelError> {
        let (cells, cell_map) = admit_physical_cell_map(source, target, cells)
            .map_err(TerminalRelativePhysicalRelabelError::CellRelabel)?;
        let semantics = canonicalize_semantic_permutation(source, target, semantics)?;
        verify_physical_component_map(source, target, &cell_map, &semantics)?;
        verify_physical_phase_map(source, target, &semantics)?;
        verify_physical_terminal_map(source, target, &cell_map, &semantics)?;

        let cell_payload = canonical_signed_relabel_payload(&cells)
            .map_err(TerminalRelativePhysicalRelabelError::CellRelabel)?;
        let component_payload = canonical_component_relabel_payload(&semantics.components)?;
        let phase_payload = canonical_phase_relabel_payload(&semantics.phases)?;
        let terminal_payload = canonical_terminal_relabel_payload(&semantics.terminals)?;
        let identity_receipt = CanonicalEncoder::<TerminalRelativePhysicalRelabelId, _>::new(
            PHYSICAL_RELABEL_IDENTITY_LIMITS,
            NeverCancel,
        )?
        .child(Field::new(0, "source-pair"), source.identity())?
        .child(Field::new(1, "target-pair"), target.identity())?
        .bytes(Field::new(2, "signed-cell-map"), &cell_payload)?
        .bytes(Field::new(3, "component-map"), &component_payload)?
        .bytes(Field::new(4, "phase-map"), &phase_payload)?
        .bytes(Field::new(5, "terminal-map"), &terminal_payload)?
        .finish()?;

        Ok(Self {
            source_pair: source.identity(),
            target_pair: target.identity(),
            cells,
            semantics,
            identity_receipt,
        })
    }

    /// Construct the explicit identity relabeling of one pair.
    pub fn identity_on(
        pair: &TerminalRelativePair,
    ) -> Result<Self, TerminalRelativePhysicalRelabelError> {
        let mut cells = Vec::new();
        for (degree, count) in pair.complex.cell_counts.iter().copied().enumerate() {
            let degree =
                u8::try_from(degree).expect("admitted terminal-relative dimension always fits u8");
            for ordinal in 0..count {
                let cell = CellRef::new(degree, ordinal);
                cells.push(SignedCellRelabelEntry::new(
                    cell,
                    cell,
                    IncidenceSign::Positive,
                ));
            }
        }
        let components = pair
            .components
            .iter()
            .map(|component| ComponentRelabelEntry::new(component.id.clone(), component.id.clone()))
            .collect();
        let phases = pair
            .phase_components
            .keys()
            .map(|phase| {
                PhaseRelabelEntry::new(phase.clone(), phase.clone(), PhaseCurrentSign::Preserve)
            })
            .collect();
        let terminals = pair
            .terminals
            .iter()
            .map(|terminal| TerminalRelabelEntry::new(terminal.id.clone(), terminal.id.clone()))
            .collect();
        Self::try_new(
            pair,
            pair,
            cells,
            TerminalRelativeSemanticPermutation::new(components, phases, terminals),
        )
    }

    /// Strong identity of this exact directed physical relabeling.
    #[must_use]
    pub const fn identity(&self) -> TerminalRelativePhysicalRelabelId {
        self.identity_receipt.id()
    }

    /// Strong identity plus canonical-preimage/schema audit roots.
    #[must_use]
    pub const fn identity_receipt(&self) -> IdentityReceipt<TerminalRelativePhysicalRelabelId> {
        self.identity_receipt
    }

    /// Strong identity of the admitted source pair.
    #[must_use]
    pub const fn source_pair_id(&self) -> TerminalRelativePairId {
        self.source_pair
    }

    /// Strong identity of the admitted target pair.
    #[must_use]
    pub const fn target_pair_id(&self) -> TerminalRelativePairId {
        self.target_pair
    }

    /// Canonically source-sorted complete signed cell map.
    #[must_use]
    pub fn cell_entries(&self) -> &[SignedCellRelabelEntry] {
        &self.cells
    }

    /// Canonically source-sorted semantic maps.
    #[must_use]
    pub const fn semantic_permutation(&self) -> &TerminalRelativeSemanticPermutation {
        &self.semantics
    }

    /// Target cell and orientation sign for one source cell.
    #[must_use]
    pub fn cell_image(&self, source: CellRef) -> Option<(CellRef, IncidenceSign)> {
        self.cells
            .binary_search_by_key(&source, |entry| entry.source)
            .ok()
            .map(|index| {
                let entry = self.cells[index];
                (entry.target, entry.sign)
            })
    }

    /// Target component for one source component.
    #[must_use]
    pub fn component_image(&self, source: &ConductorComponentId) -> Option<&ConductorComponentId> {
        self.semantics
            .components
            .binary_search_by(|entry| entry.source.cmp(source))
            .ok()
            .map(|index| &self.semantics.components[index].target)
    }

    /// Target phase and explicit current action for one source phase.
    #[must_use]
    pub fn phase_image(&self, source: &PhaseId) -> Option<(&PhaseId, PhaseCurrentSign)> {
        self.semantics
            .phases
            .binary_search_by(|entry| entry.source.cmp(source))
            .ok()
            .map(|index| {
                let entry = &self.semantics.phases[index];
                (&entry.target, entry.current_sign)
            })
    }

    /// Target terminal for one source terminal.
    #[must_use]
    pub fn terminal_image(&self, source: &PhysicalTerminalId) -> Option<&PhysicalTerminalId> {
        self.semantics
            .terminals
            .binary_search_by(|entry| entry.source.cmp(source))
            .ok()
            .map(|index| &self.semantics.terminals[index].target)
    }

    /// Transport an arbitrary integral chain by the cell action and phase
    /// permutation only. The physical current sign is intentionally not
    /// applied to a generic algebraic chain.
    pub fn transport_integral_chain(
        &self,
        source: &TerminalRelativePair,
        target: &TerminalRelativePair,
        chain: &IntegralRelativeChain,
    ) -> Result<IntegralRelativeChain, TerminalRelativePhysicalRelabelError> {
        self.verify_pair_bindings(source, target)?;
        self.verify_value_pair("integral chain source", chain.pair)?;
        let (target_phase, _) = self.required_phase_image(&chain.phase)?;
        let target_phase = target_phase.clone();
        let coefficients = self.transport_integral_coefficients(
            source,
            target,
            &chain.phase,
            &target_phase,
            chain.degree,
            &chain.coefficients,
            PhaseCurrentSign::Preserve,
        )?;
        IntegralRelativeChain::try_new(target, target_phase, chain.degree, coefficients)
            .map_err(Into::into)
    }

    /// Transport an arbitrary integral cochain by the cell action and phase
    /// permutation only. The physical current sign is intentionally not
    /// applied to a generic algebraic cochain.
    pub fn transport_integral_cochain(
        &self,
        source: &TerminalRelativePair,
        target: &TerminalRelativePair,
        cochain: &IntegralRelativeCochain,
    ) -> Result<IntegralRelativeCochain, TerminalRelativePhysicalRelabelError> {
        self.verify_pair_bindings(source, target)?;
        self.verify_value_pair("integral cochain source", cochain.pair)?;
        let (target_phase, _) = self.required_phase_image(&cochain.phase)?;
        let target_phase = target_phase.clone();
        let coefficients = self.transport_integral_coefficients(
            source,
            target,
            &cochain.phase,
            &target_phase,
            cochain.degree,
            &cochain.coefficients,
            PhaseCurrentSign::Preserve,
        )?;
        IntegralRelativeCochain::try_new(target, target_phase, cochain.degree, coefficients)
            .map_err(Into::into)
    }

    /// Transport a winding representative through the combined cell and
    /// phase-current sign in one exact coefficient operation.
    ///
    /// Combining the two signs before touching a coefficient ensures that two
    /// reversals preserve `i64::MIN` instead of spuriously overflowing through
    /// an intermediate negation.
    pub fn transport_winding_representative(
        &self,
        source: &TerminalRelativePair,
        target: &TerminalRelativePair,
        representative: &IntegralWindingRepresentative,
    ) -> Result<IntegralWindingRepresentative, TerminalRelativePhysicalRelabelError> {
        self.verify_pair_bindings(source, target)?;
        self.verify_value_pair("winding representative source", representative.chain.pair)?;
        let (target_phase, current_sign) =
            self.required_phase_image(&representative.chain.phase)?;
        let target_phase = target_phase.clone();
        let coefficients = self.transport_integral_coefficients(
            source,
            target,
            &representative.chain.phase,
            &target_phase,
            representative.chain.degree,
            &representative.chain.coefficients,
            current_sign,
        )?;
        IntegralWindingRepresentative::try_new(target, target_phase, coefficients)
            .map_err(Into::into)
    }

    /// Transport one finite scalar current amplitude with the phase's explicit
    /// current action. The caller must declare the target physical-object ID;
    /// this API never silently reuses or invents object identity.
    pub fn transport_current_amplitude(
        &self,
        source: &TerminalRelativePair,
        target: &TerminalRelativePair,
        amplitude: &RealCurrentAmplitude,
        target_id: PhysicalObjectId,
    ) -> Result<RealCurrentAmplitude, TerminalRelativePhysicalRelabelError> {
        self.verify_pair_bindings(source, target)?;
        self.verify_value_pair("current amplitude source", amplitude.pair)?;
        let (target_phase, current_sign) = self.required_phase_image(&amplitude.phase)?;
        RealCurrentAmplitude::try_new(
            target_id,
            target,
            target_phase.clone(),
            current_sign.apply(amplitude.value),
        )
        .map_err(Into::into)
    }

    /// Transport a distributed physical current through the signed cell action
    /// while requiring fresh target-side constraint receipt identifiers.
    ///
    /// A distributed current is already a physical cochain, so the phase's
    /// scalar current-coordinate sign is not applied a second time. Freshness
    /// here is deliberately nominal and fail-closed: neither target receipt
    /// may reuse either source receipt. This method does not verify the external
    /// receipt authorities or transport a current-realization map; callers must
    /// admit those target artifacts separately.
    #[allow(clippy::too_many_arguments)]
    pub fn transport_distributed_current(
        &self,
        source: &TerminalRelativePair,
        target: &TerminalRelativePair,
        current: &DistributedCurrent,
        target_id: PhysicalObjectId,
        target_divergence_receipt: StableId,
        target_terminal_constraint_receipt: StableId,
    ) -> Result<DistributedCurrent, TerminalRelativePhysicalRelabelError> {
        self.verify_pair_bindings(source, target)?;
        self.verify_value_pair("distributed current source", current.cochain.pair)?;
        let (target_phase, _) = self.required_phase_image(&current.cochain.phase)?;
        let target_phase = target_phase.clone();

        for (role, receipt) in [
            ("divergence", &target_divergence_receipt),
            ("terminal constraint", &target_terminal_constraint_receipt),
        ] {
            if receipt == current.divergence_receipt()
                || receipt == current.terminal_constraint_receipt()
            {
                return Err(
                    TerminalRelativePhysicalRelabelError::ConstraintReceiptNotFresh {
                        role,
                        receipt: receipt.as_str().to_owned(),
                    },
                );
            }
        }

        let values = self.transport_real_coefficients(
            source,
            target,
            &current.cochain.phase,
            &target_phase,
            current.cochain.degree,
            &current.cochain.values,
        )?;
        let cochain = RealRelativeCochain::try_new(
            target,
            target_phase,
            current.cochain.degree,
            current.cochain.units,
            values,
        )?;
        DistributedCurrent::new(
            target_id,
            cochain,
            target_divergence_receipt,
            target_terminal_constraint_receipt,
        )
        .map_err(Into::into)
    }

    /// Re-declare a geometric coil on the mapped target phase/component using
    /// fresh target-side realization artifacts.
    ///
    /// This operation intentionally does not claim to transport geometry. It
    /// maps only the admitted nominal phase/component bindings; the caller must
    /// supply new connectivity and manufacturing artifacts, neither of which
    /// may reuse either source artifact. A target winding-realization map must
    /// also be declared separately.
    #[allow(clippy::too_many_arguments)]
    pub fn redeclare_geometric_coil(
        &self,
        source: &TerminalRelativePair,
        target: &TerminalRelativePair,
        coil: &GeometricCoil,
        target_id: PhysicalObjectId,
        target_connectivity_artifact: StableId,
        target_manufacturing_artifact: StableId,
    ) -> Result<GeometricCoil, TerminalRelativePhysicalRelabelError> {
        self.verify_pair_bindings(source, target)?;
        self.verify_value_pair("geometric coil source", coil.pair)?;
        let (target_phase, _) = self.required_phase_image(&coil.phase)?;
        let target_phase = target_phase.clone();
        let target_component = self.component_image(&coil.component).ok_or_else(|| {
            TerminalRelativePhysicalRelabelError::MissingSemanticImage {
                role: "conductor component",
                id: coil.component.as_str().to_owned(),
            }
        })?;
        let target_component = target_component.clone();

        for (role, artifact) in [
            ("connectivity", &target_connectivity_artifact),
            ("manufacturing", &target_manufacturing_artifact),
        ] {
            if artifact == coil.connectivity_artifact() || artifact == coil.manufacturing_artifact()
            {
                return Err(
                    TerminalRelativePhysicalRelabelError::GeometricCoilArtifactNotFresh {
                        role,
                        artifact: artifact.as_str().to_owned(),
                    },
                );
            }
        }

        GeometricCoil::try_new(
            target_id,
            target,
            target_phase,
            target_component,
            target_connectivity_artifact,
            target_manufacturing_artifact,
        )
        .map_err(Into::into)
    }

    /// Admit the exact inverse physical relabeling.
    pub fn inverse(
        &self,
        source: &TerminalRelativePair,
        target: &TerminalRelativePair,
    ) -> Result<Self, TerminalRelativePhysicalRelabelError> {
        self.verify_pair_bindings(source, target)?;
        let cells = self
            .cells
            .iter()
            .map(|entry| SignedCellRelabelEntry::new(entry.target, entry.source, entry.sign))
            .collect();
        let components = self
            .semantics
            .components
            .iter()
            .map(|entry| ComponentRelabelEntry::new(entry.target.clone(), entry.source.clone()))
            .collect();
        let phases = self
            .semantics
            .phases
            .iter()
            .map(|entry| {
                PhaseRelabelEntry::new(
                    entry.target.clone(),
                    entry.source.clone(),
                    entry.current_sign.inverse(),
                )
            })
            .collect();
        let terminals = self
            .semantics
            .terminals
            .iter()
            .map(|entry| TerminalRelabelEntry::new(entry.target.clone(), entry.source.clone()))
            .collect();
        Self::try_new(
            target,
            source,
            cells,
            TerminalRelativeSemanticPermutation::new(components, phases, terminals),
        )
    }

    /// Compose `next` after `self`, including semantic permutations and the
    /// exact composed phase-current action.
    pub fn compose(
        &self,
        next: &Self,
        source: &TerminalRelativePair,
        intermediate: &TerminalRelativePair,
        target: &TerminalRelativePair,
    ) -> Result<Self, TerminalRelativePhysicalRelabelError> {
        self.verify_pair_bindings(source, intermediate)?;
        next.verify_pair_bindings(intermediate, target)?;

        let mut cells = Vec::with_capacity(self.cells.len());
        for entry in &self.cells {
            let Some((target_cell, next_sign)) = next.cell_image(entry.target) else {
                return Err(TerminalRelativePhysicalRelabelError::MissingSemanticImage {
                    role: "cell",
                    id: format!("{}:{}", entry.target.degree, entry.target.ordinal),
                });
            };
            cells.push(SignedCellRelabelEntry::new(
                entry.source,
                target_cell,
                multiply_incidence_sign(entry.sign, next_sign),
            ));
        }

        let mut components = Vec::with_capacity(self.semantics.components.len());
        for entry in &self.semantics.components {
            let Some(target_component) = next.component_image(&entry.target) else {
                return Err(TerminalRelativePhysicalRelabelError::MissingSemanticImage {
                    role: "conductor component",
                    id: entry.target.as_str().to_owned(),
                });
            };
            components.push(ComponentRelabelEntry::new(
                entry.source.clone(),
                target_component.clone(),
            ));
        }

        let mut phases = Vec::with_capacity(self.semantics.phases.len());
        for entry in &self.semantics.phases {
            let Some((target_phase, next_sign)) = next.phase_image(&entry.target) else {
                return Err(TerminalRelativePhysicalRelabelError::MissingSemanticImage {
                    role: "phase",
                    id: entry.target.as_str().to_owned(),
                });
            };
            phases.push(PhaseRelabelEntry::new(
                entry.source.clone(),
                target_phase.clone(),
                entry.current_sign.compose(next_sign),
            ));
        }

        let mut terminals = Vec::with_capacity(self.semantics.terminals.len());
        for entry in &self.semantics.terminals {
            let Some(target_terminal) = next.terminal_image(&entry.target) else {
                return Err(TerminalRelativePhysicalRelabelError::MissingSemanticImage {
                    role: "physical terminal",
                    id: entry.target.as_str().to_owned(),
                });
            };
            terminals.push(TerminalRelabelEntry::new(
                entry.source.clone(),
                target_terminal.clone(),
            ));
        }

        Self::try_new(
            source,
            target,
            cells,
            TerminalRelativeSemanticPermutation::new(components, phases, terminals),
        )
    }

    fn verify_pair_bindings(
        &self,
        source: &TerminalRelativePair,
        target: &TerminalRelativePair,
    ) -> Result<(), TerminalRelativePhysicalRelabelError> {
        if source.identity() != self.source_pair {
            return Err(TerminalRelativePhysicalRelabelError::PairIdentityMismatch {
                role: "source pair",
                expected: self.source_pair,
                actual: source.identity(),
            });
        }
        if target.identity() != self.target_pair {
            return Err(TerminalRelativePhysicalRelabelError::PairIdentityMismatch {
                role: "target pair",
                expected: self.target_pair,
                actual: target.identity(),
            });
        }
        Ok(())
    }

    fn verify_value_pair(
        &self,
        role: &'static str,
        actual: TerminalRelativePairId,
    ) -> Result<(), TerminalRelativePhysicalRelabelError> {
        if actual != self.source_pair {
            return Err(TerminalRelativePhysicalRelabelError::PairIdentityMismatch {
                role,
                expected: self.source_pair,
                actual,
            });
        }
        Ok(())
    }

    fn required_phase_image(
        &self,
        source: &PhaseId,
    ) -> Result<(&PhaseId, PhaseCurrentSign), TerminalRelativePhysicalRelabelError> {
        self.phase_image(source).ok_or_else(|| {
            TerminalRelativePhysicalRelabelError::MissingSemanticImage {
                role: "phase",
                id: source.as_str().to_owned(),
            }
        })
    }

    #[allow(clippy::too_many_arguments)]
    fn transport_integral_coefficients(
        &self,
        source: &TerminalRelativePair,
        target: &TerminalRelativePair,
        source_phase: &PhaseId,
        target_phase: &PhaseId,
        degree: u8,
        coefficients: &[i64],
        coefficient_action: PhaseCurrentSign,
    ) -> Result<Vec<i64>, TerminalRelativePhysicalRelabelError> {
        let source_basis = source.phase_relative_basis(source_phase, degree)?;
        if source_basis.len() != coefficients.len() {
            return Err(TerminalRelativeError::CoefficientArity {
                expected: source_basis.len(),
                actual: coefficients.len(),
            }
            .into());
        }
        let target_basis = target.phase_relative_basis(target_phase, degree)?;
        let target_indices: BTreeMap<_, _> = target_basis
            .iter()
            .copied()
            .enumerate()
            .map(|(index, cell)| (cell, index))
            .collect();
        let mut transported = vec![0_i64; target_basis.len()];
        for (source_cell, coefficient) in source_basis.iter().zip(coefficients) {
            let Some((target_cell, cell_sign)) = self.cell_image(*source_cell) else {
                return Err(TerminalRelativePhysicalRelabelError::MissingSemanticImage {
                    role: "cell",
                    id: format!("{}:{}", source_cell.degree, source_cell.ordinal),
                });
            };
            let Some(target_index) = target_indices.get(&target_cell) else {
                return Err(
                    TerminalRelativePhysicalRelabelError::MappedBasisCellMissing {
                        source_phase: source_phase.as_str().to_owned(),
                        target_phase: target_phase.as_str().to_owned(),
                        degree,
                        source: *source_cell,
                        target: target_cell,
                    },
                );
            };
            let effective_sign =
                multiply_incidence_sign(cell_sign, coefficient_action.incidence_sign());
            transported[*target_index] = match effective_sign {
                IncidenceSign::Positive => *coefficient,
                IncidenceSign::Negative => coefficient.checked_neg().ok_or(
                    TerminalRelativePhysicalRelabelError::CoefficientOverflow {
                        cell: *source_cell,
                    },
                )?,
            };
        }
        Ok(transported)
    }

    #[allow(clippy::too_many_arguments)]
    fn transport_real_coefficients(
        &self,
        source: &TerminalRelativePair,
        target: &TerminalRelativePair,
        source_phase: &PhaseId,
        target_phase: &PhaseId,
        degree: u8,
        values: &[f64],
    ) -> Result<Vec<f64>, TerminalRelativePhysicalRelabelError> {
        let source_basis = source.phase_relative_basis(source_phase, degree)?;
        if source_basis.len() != values.len() {
            return Err(TerminalRelativeError::CoefficientArity {
                expected: source_basis.len(),
                actual: values.len(),
            }
            .into());
        }
        let target_basis = target.phase_relative_basis(target_phase, degree)?;
        let target_indices: BTreeMap<_, _> = target_basis
            .iter()
            .copied()
            .enumerate()
            .map(|(index, cell)| (cell, index))
            .collect();
        let mut transported = vec![0.0_f64; target_basis.len()];
        for (source_cell, value) in source_basis.iter().zip(values) {
            let Some((target_cell, cell_sign)) = self.cell_image(*source_cell) else {
                return Err(TerminalRelativePhysicalRelabelError::MissingSemanticImage {
                    role: "cell",
                    id: format!("{}:{}", source_cell.degree, source_cell.ordinal),
                });
            };
            let Some(target_index) = target_indices.get(&target_cell) else {
                return Err(
                    TerminalRelativePhysicalRelabelError::MappedBasisCellMissing {
                        source_phase: source_phase.as_str().to_owned(),
                        target_phase: target_phase.as_str().to_owned(),
                        degree,
                        source: *source_cell,
                        target: target_cell,
                    },
                );
            };
            transported[*target_index] = match cell_sign {
                IncidenceSign::Positive => *value,
                IncidenceSign::Negative => -*value,
            };
        }
        Ok(transported)
    }
}

const fn multiply_incidence_sign(left: IncidenceSign, right: IncidenceSign) -> IncidenceSign {
    match (left, right) {
        (IncidenceSign::Negative, IncidenceSign::Positive)
        | (IncidenceSign::Positive, IncidenceSign::Negative) => IncidenceSign::Negative,
        (IncidenceSign::Negative, IncidenceSign::Negative)
        | (IncidenceSign::Positive, IncidenceSign::Positive) => IncidenceSign::Positive,
    }
}

fn admit_physical_cell_map(
    source: &TerminalRelativePair,
    target: &TerminalRelativePair,
    mut entries: Vec<SignedCellRelabelEntry>,
) -> Result<
    (
        Vec<SignedCellRelabelEntry>,
        BTreeMap<CellRef, SignedCellRelabelEntry>,
    ),
    TerminalRelativeSignedRelabelError,
> {
    if source.complex.dimension != target.complex.dimension {
        return Err(
            TerminalRelativeSignedRelabelError::ComplexDimensionMismatch {
                source: source.complex.dimension,
                target: target.complex.dimension,
            },
        );
    }
    for (degree, (source_count, target_count)) in source
        .complex
        .cell_counts
        .iter()
        .zip(&target.complex.cell_counts)
        .enumerate()
    {
        if source_count != target_count {
            return Err(TerminalRelativeSignedRelabelError::CellCountMismatch {
                degree: u8::try_from(degree)
                    .expect("admitted cell-complex degree always fits in u8"),
                source: *source_count,
                target: *target_count,
            });
        }
    }
    let expected_entries = source
        .complex
        .cell_counts
        .iter()
        .fold(0_usize, |total, count| {
            total.saturating_add(usize::try_from(*count).unwrap_or(usize::MAX))
        });
    if entries.len() != expected_entries {
        return Err(TerminalRelativeSignedRelabelError::EntryCountMismatch {
            expected: expected_entries,
            actual: entries.len(),
        });
    }

    entries.sort_unstable_by_key(|entry| entry.source);
    for entry in &entries {
        if !source.complex.contains(entry.source) {
            return Err(TerminalRelativeSignedRelabelError::SourceCellOutOfRange {
                cell: entry.source,
            });
        }
        if !target.complex.contains(entry.target) {
            return Err(TerminalRelativeSignedRelabelError::TargetCellOutOfRange {
                cell: entry.target,
            });
        }
        if entry.source.degree != entry.target.degree {
            return Err(TerminalRelativeSignedRelabelError::CellDegreeMismatch {
                source: entry.source,
                target: entry.target,
            });
        }
    }
    if let Some(duplicate) = entries
        .windows(2)
        .find(|pair| pair[0].source == pair[1].source)
        .map(|pair| pair[0].source)
    {
        return Err(TerminalRelativeSignedRelabelError::DuplicateSourceCell { cell: duplicate });
    }
    let mut target_cells = BTreeSet::new();
    for entry in &entries {
        if !target_cells.insert(entry.target) {
            return Err(TerminalRelativeSignedRelabelError::DuplicateTargetCell {
                cell: entry.target,
            });
        }
    }

    let entry_map: BTreeMap<_, _> = entries
        .iter()
        .copied()
        .map(|entry| (entry.source, entry))
        .collect();
    verify_signed_incidence_commutation(source, target, &entry_map)?;
    verify_mapped_subcomplex(
        "conductor",
        None,
        &source.conductor,
        &target.conductor,
        &entry_map,
    )?;
    verify_mapped_subcomplex(
        "relative subcomplex",
        None,
        &source.relative,
        &target.relative,
        &entry_map,
    )?;
    verify_mapped_subcomplex(
        "insulation",
        None,
        &source.insulation,
        &target.insulation,
        &entry_map,
    )?;
    Ok((entries, entry_map))
}

fn incidence_coordinate(incidence: &BoundaryIncidence) -> (CellRef, CellRef) {
    (incidence.upper, incidence.lower)
}

fn verify_signed_incidence_commutation(
    source: &TerminalRelativePair,
    target: &TerminalRelativePair,
    entries: &BTreeMap<CellRef, SignedCellRelabelEntry>,
) -> Result<(), TerminalRelativeSignedRelabelError> {
    let mut mapped = Vec::with_capacity(source.complex.incidences.len());
    for incidence in &source.complex.incidences {
        let lower = entries.get(&incidence.lower).ok_or(
            TerminalRelativeSignedRelabelError::MissingSourceCell {
                cell: incidence.lower,
            },
        )?;
        let upper = entries.get(&incidence.upper).ok_or(
            TerminalRelativeSignedRelabelError::MissingSourceCell {
                cell: incidence.upper,
            },
        )?;
        let sign = multiply_incidence_sign(
            multiply_incidence_sign(incidence.sign, lower.sign),
            upper.sign,
        );
        mapped.push(BoundaryIncidence::new(lower.target, upper.target, sign));
    }
    mapped.sort_unstable_by_key(|incidence| {
        (
            incidence.upper.degree,
            incidence.upper.ordinal,
            incidence.lower.ordinal,
        )
    });

    let expected = &mapped;
    let actual = &target.complex.incidences;
    let mut expected_index = 0;
    let mut actual_index = 0;
    while expected_index < expected.len() || actual_index < actual.len() {
        match (expected.get(expected_index), actual.get(actual_index)) {
            (Some(expected_entry), Some(actual_entry)) => {
                let expected_coordinate = incidence_coordinate(expected_entry);
                let actual_coordinate = incidence_coordinate(actual_entry);
                match expected_coordinate.cmp(&actual_coordinate) {
                    core::cmp::Ordering::Less => {
                        return Err(
                            TerminalRelativeSignedRelabelError::MappedIncidenceMismatch {
                                lower: expected_entry.lower,
                                upper: expected_entry.upper,
                                expected: Some(expected_entry.sign),
                                actual: None,
                            },
                        );
                    }
                    core::cmp::Ordering::Greater => {
                        return Err(
                            TerminalRelativeSignedRelabelError::MappedIncidenceMismatch {
                                lower: actual_entry.lower,
                                upper: actual_entry.upper,
                                expected: None,
                                actual: Some(actual_entry.sign),
                            },
                        );
                    }
                    core::cmp::Ordering::Equal => {
                        if expected_entry.sign != actual_entry.sign {
                            return Err(
                                TerminalRelativeSignedRelabelError::MappedIncidenceMismatch {
                                    lower: actual_entry.lower,
                                    upper: actual_entry.upper,
                                    expected: Some(expected_entry.sign),
                                    actual: Some(actual_entry.sign),
                                },
                            );
                        }
                        expected_index += 1;
                        actual_index += 1;
                    }
                }
            }
            (Some(expected_entry), None) => {
                return Err(
                    TerminalRelativeSignedRelabelError::MappedIncidenceMismatch {
                        lower: expected_entry.lower,
                        upper: expected_entry.upper,
                        expected: Some(expected_entry.sign),
                        actual: None,
                    },
                );
            }
            (None, Some(actual_entry)) => {
                return Err(
                    TerminalRelativeSignedRelabelError::MappedIncidenceMismatch {
                        lower: actual_entry.lower,
                        upper: actual_entry.upper,
                        expected: None,
                        actual: Some(actual_entry.sign),
                    },
                );
            }
            (None, None) => break,
        }
    }
    Ok(())
}

fn verify_mapped_subcomplex(
    role: &'static str,
    owner: Option<&str>,
    source: &CellularSubcomplex,
    target: &CellularSubcomplex,
    entries: &BTreeMap<CellRef, SignedCellRelabelEntry>,
) -> Result<(), TerminalRelativeSignedRelabelError> {
    if source.id != target.id {
        return Err(
            TerminalRelativeSignedRelabelError::SubcomplexIdentityMismatch {
                role,
                owner: owner.map(str::to_owned),
                source: source.id.as_str().to_owned(),
                target: target.id.as_str().to_owned(),
            },
        );
    }
    let mut mapped = BTreeSet::new();
    for source_cell in &source.cells {
        let entry = entries
            .get(source_cell)
            .ok_or(TerminalRelativeSignedRelabelError::MissingSourceCell { cell: *source_cell })?;
        mapped.insert(entry.target);
    }
    if mapped != target.cells {
        let cell = mapped
            .symmetric_difference(&target.cells)
            .next()
            .copied()
            .expect("unequal finite support sets have a witness");
        return Err(TerminalRelativeSignedRelabelError::MappedSupportMismatch {
            role,
            owner: owner.map(str::to_owned),
            cell,
            expected_mapped: mapped.contains(&cell),
            actual_target: target.cells.contains(&cell),
        });
    }
    Ok(())
}

fn verify_component_semantics(
    source: &TerminalRelativePair,
    target: &TerminalRelativePair,
    entries: &BTreeMap<CellRef, SignedCellRelabelEntry>,
) -> Result<(), TerminalRelativeSignedRelabelError> {
    let source_components: BTreeMap<_, _> = source
        .components
        .iter()
        .map(|component| (component.id.as_str(), component))
        .collect();
    let target_components: BTreeMap<_, _> = target
        .components
        .iter()
        .map(|component| (component.id.as_str(), component))
        .collect();
    for source_component in &source.components {
        let id = source_component.id.as_str();
        let Some(target_component) = target_components.get(id).copied() else {
            return Err(
                TerminalRelativeSignedRelabelError::SemanticIdentitySetMismatch {
                    role: "conductor component",
                    id: id.to_owned(),
                    source_present: true,
                    target_present: false,
                },
            );
        };
        verify_mapped_subcomplex(
            "conductor component support",
            Some(id),
            &source_component.support,
            &target_component.support,
            entries,
        )?;
    }
    if let Some(id) = target_components
        .keys()
        .copied()
        .find(|id| !source_components.contains_key(id))
    {
        return Err(
            TerminalRelativeSignedRelabelError::SemanticIdentitySetMismatch {
                role: "conductor component",
                id: id.to_owned(),
                source_present: false,
                target_present: true,
            },
        );
    }
    Ok(())
}

fn verify_phase_semantics(
    source: &TerminalRelativePair,
    target: &TerminalRelativePair,
) -> Result<(), TerminalRelativeSignedRelabelError> {
    for (phase, source_component) in &source.phase_components {
        let target_component = target.phase_components.get(phase);
        if target_component != Some(source_component) {
            return Err(
                TerminalRelativeSignedRelabelError::PhaseComponentBindingMismatch {
                    phase: phase.as_str().to_owned(),
                    source_component: Some(source_component.as_str().to_owned()),
                    target_component: target_component
                        .map(|component| component.as_str().to_owned()),
                },
            );
        }
    }
    if let Some((phase, target_component)) = target
        .phase_components
        .iter()
        .find(|(phase, _)| !source.phase_components.contains_key(*phase))
    {
        return Err(
            TerminalRelativeSignedRelabelError::PhaseComponentBindingMismatch {
                phase: phase.as_str().to_owned(),
                source_component: None,
                target_component: Some(target_component.as_str().to_owned()),
            },
        );
    }
    Ok(())
}

fn verify_terminal_semantics(
    source: &TerminalRelativePair,
    target: &TerminalRelativePair,
    entries: &BTreeMap<CellRef, SignedCellRelabelEntry>,
) -> Result<(), TerminalRelativeSignedRelabelError> {
    let source_terminals: BTreeMap<_, _> = source
        .terminals
        .iter()
        .map(|terminal| (terminal.id.as_str(), terminal))
        .collect();
    let target_terminals: BTreeMap<_, _> = target
        .terminals
        .iter()
        .map(|terminal| (terminal.id.as_str(), terminal))
        .collect();
    for source_terminal in &source.terminals {
        let id = source_terminal.id.as_str();
        let Some(target_terminal) = target_terminals.get(id).copied() else {
            return Err(
                TerminalRelativeSignedRelabelError::SemanticIdentitySetMismatch {
                    role: "physical terminal",
                    id: id.to_owned(),
                    source_present: true,
                    target_present: false,
                },
            );
        };
        verify_mapped_subcomplex(
            "physical terminal support",
            Some(id),
            &source_terminal.support,
            &target_terminal.support,
            entries,
        )?;
        let mismatched_field = if source_terminal.component != target_terminal.component {
            Some("component")
        } else if source_terminal.phase != target_terminal.phase {
            Some("phase")
        } else if source_terminal.role != target_terminal.role {
            Some("role")
        } else if source_terminal.orientation != target_terminal.orientation {
            Some("orientation")
        } else if source_terminal.coordinate != target_terminal.coordinate {
            Some("coordinate")
        } else if source_terminal.port != target_terminal.port {
            Some("port schema")
        } else if source_terminal.machine != target_terminal.machine {
            Some("presented Machine-IR binding")
        } else if source_terminal.trivialization != target_terminal.trivialization {
            Some("port trivialization")
        } else {
            None
        };
        if let Some(field) = mismatched_field {
            return Err(
                TerminalRelativeSignedRelabelError::TerminalMetadataMismatch {
                    terminal: id.to_owned(),
                    field,
                },
            );
        }
    }
    if let Some(id) = target_terminals
        .keys()
        .copied()
        .find(|id| !source_terminals.contains_key(id))
    {
        return Err(
            TerminalRelativeSignedRelabelError::SemanticIdentitySetMismatch {
                role: "physical terminal",
                id: id.to_owned(),
                source_present: false,
                target_present: true,
            },
        );
    }
    Ok(())
}

fn canonicalize_semantic_permutation(
    source: &TerminalRelativePair,
    target: &TerminalRelativePair,
    mut semantics: TerminalRelativeSemanticPermutation,
) -> Result<TerminalRelativeSemanticPermutation, TerminalRelativePhysicalRelabelError> {
    if source.components.len() != target.components.len() {
        return Err(
            TerminalRelativePhysicalRelabelError::SemanticCardinalityMismatch {
                role: "conductor component",
                source: source.components.len(),
                target: target.components.len(),
            },
        );
    }
    if semantics.components.len() != source.components.len() {
        return Err(
            TerminalRelativePhysicalRelabelError::SemanticEntryCountMismatch {
                role: "conductor component",
                expected: source.components.len(),
                actual: semantics.components.len(),
            },
        );
    }
    semantics
        .components
        .sort_unstable_by(|left, right| left.source.cmp(&right.source));
    let source_components: BTreeSet<_> = source
        .components
        .iter()
        .map(|component| component.id.clone())
        .collect();
    let target_components: BTreeSet<_> = target
        .components
        .iter()
        .map(|component| component.id.clone())
        .collect();
    for entry in &semantics.components {
        if !source_components.contains(&entry.source) {
            return Err(
                TerminalRelativePhysicalRelabelError::UnknownSemanticSource {
                    role: "conductor component",
                    id: entry.source.as_str().to_owned(),
                },
            );
        }
        if !target_components.contains(&entry.target) {
            return Err(
                TerminalRelativePhysicalRelabelError::UnknownSemanticTarget {
                    role: "conductor component",
                    id: entry.target.as_str().to_owned(),
                },
            );
        }
    }
    if let Some(duplicate) = semantics
        .components
        .windows(2)
        .find(|pair| pair[0].source == pair[1].source)
        .map(|pair| pair[0].source.as_str().to_owned())
    {
        return Err(
            TerminalRelativePhysicalRelabelError::DuplicateSemanticSource {
                role: "conductor component",
                id: duplicate,
            },
        );
    }
    let mut component_targets = BTreeSet::new();
    for entry in &semantics.components {
        if !component_targets.insert(entry.target.clone()) {
            return Err(
                TerminalRelativePhysicalRelabelError::DuplicateSemanticTarget {
                    role: "conductor component",
                    id: entry.target.as_str().to_owned(),
                },
            );
        }
    }

    if source.phase_components.len() != target.phase_components.len() {
        return Err(
            TerminalRelativePhysicalRelabelError::SemanticCardinalityMismatch {
                role: "phase",
                source: source.phase_components.len(),
                target: target.phase_components.len(),
            },
        );
    }
    if semantics.phases.len() != source.phase_components.len() {
        return Err(
            TerminalRelativePhysicalRelabelError::SemanticEntryCountMismatch {
                role: "phase",
                expected: source.phase_components.len(),
                actual: semantics.phases.len(),
            },
        );
    }
    semantics
        .phases
        .sort_unstable_by(|left, right| left.source.cmp(&right.source));
    let source_phases: BTreeSet<_> = source.phase_components.keys().cloned().collect();
    let target_phases: BTreeSet<_> = target.phase_components.keys().cloned().collect();
    for entry in &semantics.phases {
        if !source_phases.contains(&entry.source) {
            return Err(
                TerminalRelativePhysicalRelabelError::UnknownSemanticSource {
                    role: "phase",
                    id: entry.source.as_str().to_owned(),
                },
            );
        }
        if !target_phases.contains(&entry.target) {
            return Err(
                TerminalRelativePhysicalRelabelError::UnknownSemanticTarget {
                    role: "phase",
                    id: entry.target.as_str().to_owned(),
                },
            );
        }
    }
    if let Some(duplicate) = semantics
        .phases
        .windows(2)
        .find(|pair| pair[0].source == pair[1].source)
        .map(|pair| pair[0].source.as_str().to_owned())
    {
        return Err(
            TerminalRelativePhysicalRelabelError::DuplicateSemanticSource {
                role: "phase",
                id: duplicate,
            },
        );
    }
    let mut phase_targets = BTreeSet::new();
    for entry in &semantics.phases {
        if !phase_targets.insert(entry.target.clone()) {
            return Err(
                TerminalRelativePhysicalRelabelError::DuplicateSemanticTarget {
                    role: "phase",
                    id: entry.target.as_str().to_owned(),
                },
            );
        }
    }

    if source.terminals.len() != target.terminals.len() {
        return Err(
            TerminalRelativePhysicalRelabelError::SemanticCardinalityMismatch {
                role: "physical terminal",
                source: source.terminals.len(),
                target: target.terminals.len(),
            },
        );
    }
    if semantics.terminals.len() != source.terminals.len() {
        return Err(
            TerminalRelativePhysicalRelabelError::SemanticEntryCountMismatch {
                role: "physical terminal",
                expected: source.terminals.len(),
                actual: semantics.terminals.len(),
            },
        );
    }
    semantics
        .terminals
        .sort_unstable_by(|left, right| left.source.cmp(&right.source));
    let source_terminals: BTreeSet<_> = source
        .terminals
        .iter()
        .map(|terminal| terminal.id.clone())
        .collect();
    let target_terminals: BTreeSet<_> = target
        .terminals
        .iter()
        .map(|terminal| terminal.id.clone())
        .collect();
    for entry in &semantics.terminals {
        if !source_terminals.contains(&entry.source) {
            return Err(
                TerminalRelativePhysicalRelabelError::UnknownSemanticSource {
                    role: "physical terminal",
                    id: entry.source.as_str().to_owned(),
                },
            );
        }
        if !target_terminals.contains(&entry.target) {
            return Err(
                TerminalRelativePhysicalRelabelError::UnknownSemanticTarget {
                    role: "physical terminal",
                    id: entry.target.as_str().to_owned(),
                },
            );
        }
    }
    if let Some(duplicate) = semantics
        .terminals
        .windows(2)
        .find(|pair| pair[0].source == pair[1].source)
        .map(|pair| pair[0].source.as_str().to_owned())
    {
        return Err(
            TerminalRelativePhysicalRelabelError::DuplicateSemanticSource {
                role: "physical terminal",
                id: duplicate,
            },
        );
    }
    let mut terminal_targets = BTreeSet::new();
    for entry in &semantics.terminals {
        if !terminal_targets.insert(entry.target.clone()) {
            return Err(
                TerminalRelativePhysicalRelabelError::DuplicateSemanticTarget {
                    role: "physical terminal",
                    id: entry.target.as_str().to_owned(),
                },
            );
        }
    }

    Ok(semantics)
}

fn verify_physical_component_map(
    source: &TerminalRelativePair,
    target: &TerminalRelativePair,
    cells: &BTreeMap<CellRef, SignedCellRelabelEntry>,
    semantics: &TerminalRelativeSemanticPermutation,
) -> Result<(), TerminalRelativePhysicalRelabelError> {
    let source_components: BTreeMap<_, _> = source
        .components
        .iter()
        .map(|component| (component.id.clone(), component))
        .collect();
    let target_components: BTreeMap<_, _> = target
        .components
        .iter()
        .map(|component| (component.id.clone(), component))
        .collect();
    for entry in &semantics.components {
        let source_component = source_components
            .get(&entry.source)
            .expect("canonical semantic admission retains every source component");
        let target_component = target_components
            .get(&entry.target)
            .expect("canonical semantic admission retains every target component");
        verify_physical_mapped_support(
            "conductor component support",
            entry.source.as_str(),
            entry.target.as_str(),
            &source_component.support,
            &target_component.support,
            cells,
        )?;
    }
    Ok(())
}

fn verify_physical_phase_map(
    source: &TerminalRelativePair,
    target: &TerminalRelativePair,
    semantics: &TerminalRelativeSemanticPermutation,
) -> Result<(), TerminalRelativePhysicalRelabelError> {
    let component_images: BTreeMap<_, _> = semantics
        .components
        .iter()
        .map(|entry| (entry.source.clone(), entry.target.clone()))
        .collect();
    for entry in &semantics.phases {
        let source_component = source
            .phase_components
            .get(&entry.source)
            .expect("canonical semantic admission retains every source phase");
        let target_component = target
            .phase_components
            .get(&entry.target)
            .expect("canonical semantic admission retains every target phase");
        let expected_target = component_images
            .get(source_component)
            .expect("total component permutation has every phase component");
        if target_component != expected_target {
            return Err(
                TerminalRelativePhysicalRelabelError::PhaseComponentSquareMismatch {
                    source_phase: entry.source.as_str().to_owned(),
                    target_phase: entry.target.as_str().to_owned(),
                    expected_target_component: expected_target.as_str().to_owned(),
                    actual_target_component: target_component.as_str().to_owned(),
                },
            );
        }
    }
    Ok(())
}

fn verify_physical_terminal_map(
    source: &TerminalRelativePair,
    target: &TerminalRelativePair,
    cells: &BTreeMap<CellRef, SignedCellRelabelEntry>,
    semantics: &TerminalRelativeSemanticPermutation,
) -> Result<(), TerminalRelativePhysicalRelabelError> {
    let source_terminals: BTreeMap<_, _> = source
        .terminals
        .iter()
        .map(|terminal| (terminal.id.clone(), terminal))
        .collect();
    let target_terminals: BTreeMap<_, _> = target
        .terminals
        .iter()
        .map(|terminal| (terminal.id.clone(), terminal))
        .collect();
    let component_images: BTreeMap<_, _> = semantics
        .components
        .iter()
        .map(|entry| (entry.source.clone(), entry.target.clone()))
        .collect();
    let phase_images: BTreeMap<_, _> = semantics
        .phases
        .iter()
        .map(|entry| {
            (
                entry.source.clone(),
                (entry.target.clone(), entry.current_sign),
            )
        })
        .collect();

    for entry in &semantics.terminals {
        let source_terminal = source_terminals
            .get(&entry.source)
            .expect("canonical semantic admission retains every source terminal");
        let target_terminal = target_terminals
            .get(&entry.target)
            .expect("canonical semantic admission retains every target terminal");
        verify_physical_mapped_support(
            "physical terminal support",
            entry.source.as_str(),
            entry.target.as_str(),
            &source_terminal.support,
            &target_terminal.support,
            cells,
        )?;

        let (expected_phase, current_sign) = phase_images
            .get(&source_terminal.phase)
            .expect("total phase permutation has every terminal phase");
        if &target_terminal.phase != expected_phase {
            return Err(
                TerminalRelativePhysicalRelabelError::TerminalPhaseSquareMismatch {
                    source_terminal: entry.source.as_str().to_owned(),
                    target_terminal: entry.target.as_str().to_owned(),
                    expected_target_phase: expected_phase.as_str().to_owned(),
                    actual_target_phase: target_terminal.phase.as_str().to_owned(),
                },
            );
        }
        let expected_component = component_images
            .get(&source_terminal.component)
            .expect("total component permutation has every terminal component");
        if &target_terminal.component != expected_component {
            return Err(
                TerminalRelativePhysicalRelabelError::TerminalComponentSquareMismatch {
                    source_terminal: entry.source.as_str().to_owned(),
                    target_terminal: entry.target.as_str().to_owned(),
                    expected_target_component: expected_component.as_str().to_owned(),
                    actual_target_component: target_terminal.component.as_str().to_owned(),
                },
            );
        }

        let mismatch = if target_terminal.role != current_sign.map_role(source_terminal.role) {
            Some("terminal role/current-sign parity")
        } else if target_terminal.orientation
            != current_sign.map_orientation(source_terminal.orientation)
        {
            Some("terminal orientation/current-sign parity")
        } else if target_terminal.trivialization.sign()
            != current_sign.map_trivialization(source_terminal.trivialization.sign())
        {
            Some("terminal trivialization/current-sign parity")
        } else {
            physical_terminal_convention_mismatch(source_terminal, target_terminal)
        };
        if let Some(field) = mismatch {
            return Err(
                TerminalRelativePhysicalRelabelError::TerminalConventionMismatch {
                    source_terminal: entry.source.as_str().to_owned(),
                    target_terminal: entry.target.as_str().to_owned(),
                    field,
                },
            );
        }
    }
    Ok(())
}

fn verify_physical_mapped_support(
    role: &'static str,
    source_owner: &str,
    target_owner: &str,
    source: &CellularSubcomplex,
    target: &CellularSubcomplex,
    cells: &BTreeMap<CellRef, SignedCellRelabelEntry>,
) -> Result<(), TerminalRelativePhysicalRelabelError> {
    let mut mapped = BTreeSet::new();
    for source_cell in &source.cells {
        let entry = cells.get(source_cell).ok_or_else(|| {
            TerminalRelativePhysicalRelabelError::CellRelabel(
                TerminalRelativeSignedRelabelError::MissingSourceCell { cell: *source_cell },
            )
        })?;
        mapped.insert(entry.target);
    }
    if mapped != target.cells {
        let cell = mapped
            .symmetric_difference(&target.cells)
            .next()
            .copied()
            .expect("unequal finite support sets have a witness");
        return Err(
            TerminalRelativePhysicalRelabelError::MappedSemanticSupportMismatch {
                role,
                source_owner: source_owner.to_owned(),
                target_owner: target_owner.to_owned(),
                cell,
                expected_mapped: mapped.contains(&cell),
                actual_target: target.cells.contains(&cell),
            },
        );
    }
    Ok(())
}

fn physical_terminal_convention_mismatch(
    source: &PhysicalTerminal,
    target: &PhysicalTerminal,
) -> Option<&'static str> {
    let source_port = &source.port;
    let target_port = &target.port;
    if source.coordinate != target.coordinate {
        Some("terminal coordinate")
    } else if source_port.version() != target_port.version() {
        Some("port schema version")
    } else if source_port.kind() != target_port.kind() {
        Some("port kind")
    } else if source_port.effort_dimensions() != target_port.effort_dimensions() {
        Some("effort dimensions")
    } else if source_port.flow_dimensions() != target_port.flow_dimensions() {
        Some("flow dimensions")
    } else if source_port.shape() != target_port.shape() {
        Some("port value shape")
    } else if source_port.coordinates().basis() != target_port.coordinates().basis() {
        Some("coordinate basis")
    } else if source_port.coordinates().frame() != target_port.coordinates().frame() {
        Some("coordinate frame")
    } else if source_port.coordinates().orientation() != target_port.coordinates().orientation() {
        Some("coordinate orientation")
    } else if source_port.power_pairing() != target_port.power_pairing() {
        Some("power pairing")
    } else if source_port.timestamp().clock() != target_port.timestamp().clock() {
        Some("clock domain")
    } else if source_port.timestamp().tick() != target_port.timestamp().tick() {
        Some("clock tick")
    } else if source_port.conservation_roles() != target_port.conservation_roles() {
        Some("conservation roles")
    } else if source.machine.authority_domain() != target.machine.authority_domain() {
        Some("MachineGraph authority domain")
    } else if source.machine.schema_version() != target.machine.schema_version() {
        Some("MachineGraph schema version")
    } else if source.machine.graph_digest() != target.machine.graph_digest() {
        Some("MachineGraph digest")
    } else if source.machine.owner() != target.machine.owner() {
        Some("MachineGraph owner")
    } else if source.trivialization.voltage_reference() != target.trivialization.voltage_reference()
    {
        Some("voltage reference")
    } else if source.trivialization.current_reference() != target.trivialization.current_reference()
    {
        Some("current reference")
    } else {
        None
    }
}

/// Integral chain on the canonical terminal-relative quotient basis.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IntegralRelativeChain {
    pair: TerminalRelativePairId,
    phase: PhaseId,
    degree: u8,
    coefficients: Vec<i64>,
}

impl IntegralRelativeChain {
    /// Construct an integral chain.  No real-to-integer conversion exists.
    pub fn try_new(
        pair: &TerminalRelativePair,
        phase: PhaseId,
        degree: u8,
        coefficients: Vec<i64>,
    ) -> Result<Self, TerminalRelativeError> {
        if !pair.contains_phase(&phase) {
            return Err(TerminalRelativeError::UnknownPhase {
                phase: phase.as_str().to_owned(),
            });
        }
        if degree > pair.complex.dimension {
            return Err(TerminalRelativeError::DegreeOutOfRange {
                degree,
                dimension: pair.complex.dimension,
            });
        }
        let expected = pair.phase_relative_basis(&phase, degree)?.len();
        if coefficients.len() != expected {
            return Err(TerminalRelativeError::CoefficientArity {
                expected,
                actual: coefficients.len(),
            });
        }
        Ok(Self {
            pair: pair.identity(),
            phase,
            degree,
            coefficients,
        })
    }

    /// Pair identity.
    #[must_use]
    pub const fn pair_id(&self) -> TerminalRelativePairId {
        self.pair
    }

    /// Phase identity.
    #[must_use]
    pub const fn phase(&self) -> &PhaseId {
        &self.phase
    }

    /// Chain degree.
    #[must_use]
    pub const fn degree(&self) -> u8 {
        self.degree
    }

    /// Exact integral coefficients.
    #[must_use]
    pub fn coefficients(&self) -> &[i64] {
        &self.coefficients
    }
}

/// Integral cochain on the phase-local relative quotient basis.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IntegralRelativeCochain {
    pair: TerminalRelativePairId,
    phase: PhaseId,
    degree: u8,
    coefficients: Vec<i64>,
}

impl IntegralRelativeCochain {
    /// Construct an exact integral cochain. It is neither a real physical
    /// field nor a cohomology-class witness.
    pub fn try_new(
        pair: &TerminalRelativePair,
        phase: PhaseId,
        degree: u8,
        coefficients: Vec<i64>,
    ) -> Result<Self, TerminalRelativeError> {
        if !pair.contains_phase(&phase) {
            return Err(TerminalRelativeError::UnknownPhase {
                phase: phase.as_str().to_owned(),
            });
        }
        if degree > pair.complex.dimension {
            return Err(TerminalRelativeError::DegreeOutOfRange {
                degree,
                dimension: pair.complex.dimension,
            });
        }
        let expected = pair.phase_relative_basis(&phase, degree)?.len();
        if coefficients.len() != expected {
            return Err(TerminalRelativeError::CoefficientArity {
                expected,
                actual: coefficients.len(),
            });
        }
        Ok(Self {
            pair: pair.identity(),
            phase,
            degree,
            coefficients,
        })
    }

    /// Pair identity.
    #[must_use]
    pub const fn pair_id(&self) -> TerminalRelativePairId {
        self.pair
    }

    /// Phase identity.
    #[must_use]
    pub const fn phase(&self) -> &PhaseId {
        &self.phase
    }

    /// Cochain degree.
    #[must_use]
    pub const fn degree(&self) -> u8 {
        self.degree
    }

    /// Exact integral coefficients in dual-basis order.
    #[must_use]
    pub fn coefficients(&self) -> &[i64] {
        &self.coefficients
    }
}

/// Real cochain with explicit physical dimensions.
#[derive(Clone, Debug, PartialEq)]
pub struct RealRelativeCochain {
    pair: TerminalRelativePairId,
    phase: PhaseId,
    degree: u8,
    units: Dims,
    values: Vec<f64>,
}

impl RealRelativeCochain {
    /// Construct a real cochain.  It is not an integral topology class.
    pub fn try_new(
        pair: &TerminalRelativePair,
        phase: PhaseId,
        degree: u8,
        units: Dims,
        values: Vec<f64>,
    ) -> Result<Self, TerminalRelativeError> {
        if !pair.contains_phase(&phase) {
            return Err(TerminalRelativeError::UnknownPhase {
                phase: phase.as_str().to_owned(),
            });
        }
        if degree > pair.complex.dimension {
            return Err(TerminalRelativeError::DegreeOutOfRange {
                degree,
                dimension: pair.complex.dimension,
            });
        }
        let expected = pair.phase_relative_basis(&phase, degree)?.len();
        if values.len() != expected {
            return Err(TerminalRelativeError::CoefficientArity {
                expected,
                actual: values.len(),
            });
        }
        if let Some(index) = values.iter().position(|value| !value.is_finite()) {
            return Err(TerminalRelativeError::NonFiniteRealCoefficient { index });
        }
        Ok(Self {
            pair: pair.identity(),
            phase,
            degree,
            units,
            values,
        })
    }

    /// Pair identity.
    #[must_use]
    pub const fn pair_id(&self) -> TerminalRelativePairId {
        self.pair
    }

    /// Phase identity.
    #[must_use]
    pub const fn phase(&self) -> &PhaseId {
        &self.phase
    }

    /// Cochain degree.
    #[must_use]
    pub const fn degree(&self) -> u8 {
        self.degree
    }

    /// Physical dimensions of every coefficient.
    #[must_use]
    pub const fn units(&self) -> Dims {
        self.units
    }

    /// Finite real coefficients.
    #[must_use]
    pub fn values(&self) -> &[f64] {
        &self.values
    }
}

/// Integral terminal-relative one-cycle representing winding topology only.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct IntegralWindingRepresentative {
    chain: IntegralRelativeChain,
    identity_receipt: IdentityReceipt<IntegralWindingRepresentativeId>,
}

impl IntegralWindingRepresentative {
    /// Admit an integral relative one-cycle.  A nonzero relative boundary
    /// refuses; no rounding or tolerance is involved.
    pub fn try_new(
        pair: &TerminalRelativePair,
        phase: PhaseId,
        coefficients: Vec<i64>,
    ) -> Result<Self, TerminalRelativeError> {
        let chain = IntegralRelativeChain::try_new(pair, phase, 1, coefficients)?;
        let boundary = pair.boundary(&chain)?;
        if let Some((index, coefficient)) = boundary
            .coefficients
            .iter()
            .copied()
            .enumerate()
            .find(|(_, coefficient)| *coefficient != 0)
        {
            return Err(TerminalRelativeError::NotARelativeCycle { index, coefficient });
        }
        let payload = canonical_winding_payload(&chain)?;
        let receipt = CanonicalEncoder::<IntegralWindingRepresentativeId, _>::new(
            PAIR_IDENTITY_LIMITS,
            NeverCancel,
        )?
        .bytes(Field::new(0, "winding-representative-payload"), &payload)?
        .finish()?;
        Ok(Self {
            chain,
            identity_receipt: receipt,
        })
    }

    /// Strong representative identity.  It is not a homology-class witness.
    #[must_use]
    pub const fn identity(&self) -> IntegralWindingRepresentativeId {
        self.identity_receipt.id()
    }

    /// Strong identity plus canonical-preimage/schema audit roots.
    #[must_use]
    pub const fn identity_receipt(&self) -> IdentityReceipt<IntegralWindingRepresentativeId> {
        self.identity_receipt
    }

    /// Exact integral relative cycle.
    #[must_use]
    pub const fn chain(&self) -> &IntegralRelativeChain {
        &self.chain
    }

    /// Nominal reference for a declared physical conversion map.
    #[must_use]
    pub fn object_ref(&self) -> PhysicalObjectRef {
        PhysicalObjectRef {
            identity: PhysicalObjectIdentity::IntegralWindingRepresentative(self.identity()),
            pair: self.chain.pair,
            phase: self.chain.phase.clone(),
            kind: PhysicalObjectKind::IntegralWindingRepresentative,
        }
    }
}

/// Real current amplitude; distinct from an integral winding representative.
#[derive(Clone, PartialEq)]
pub struct RealCurrentAmplitude {
    id: PhysicalObjectId,
    pair: TerminalRelativePairId,
    phase: PhaseId,
    value: Current,
}

impl fmt::Debug for RealCurrentAmplitude {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("RealCurrentAmplitude")
            .field("id", &self.id)
            .field("pair", &self.pair)
            .field("phase", &self.phase)
            .field("amperes", &self.value.value())
            .finish()
    }
}

impl RealCurrentAmplitude {
    /// Bind a finite real current amplitude to one pair/phase.
    pub fn try_new(
        id: PhysicalObjectId,
        pair: &TerminalRelativePair,
        phase: PhaseId,
        value: Current,
    ) -> Result<Self, TerminalRelativeError> {
        if !pair.contains_phase(&phase) {
            return Err(TerminalRelativeError::UnknownPhase {
                phase: phase.as_str().to_owned(),
            });
        }
        if !value.value().is_finite() {
            return Err(TerminalRelativeError::NonFiniteCurrentAmplitude);
        }
        Ok(Self {
            id,
            pair: pair.identity(),
            phase,
            value,
        })
    }

    /// Current value in coherent amperes.
    #[must_use]
    pub const fn value(&self) -> Current {
        self.value
    }

    /// Declared physical object identity.
    #[must_use]
    pub const fn id(&self) -> &PhysicalObjectId {
        &self.id
    }

    /// Terminal-relative pair identity.
    #[must_use]
    pub const fn pair_id(&self) -> TerminalRelativePairId {
        self.pair
    }

    /// Electrical phase identity.
    #[must_use]
    pub const fn phase(&self) -> &PhaseId {
        &self.phase
    }

    /// Nominal object reference.
    #[must_use]
    pub fn object_ref(&self) -> PhysicalObjectRef {
        PhysicalObjectRef {
            identity: PhysicalObjectIdentity::Declared(self.id.clone()),
            pair: self.pair,
            phase: self.phase.clone(),
            kind: PhysicalObjectKind::RealCurrentAmplitude,
        }
    }
}

/// Distributed physical current cochain with explicit constraint receipts.
#[derive(Clone, Debug, PartialEq)]
pub struct DistributedCurrent {
    id: PhysicalObjectId,
    cochain: RealRelativeCochain,
    divergence_receipt: StableId,
    terminal_constraint_receipt: StableId,
}

impl DistributedCurrent {
    /// Bind a current-dimensioned real cochain and the exact external receipts
    /// that claim divergence and terminal closure.
    pub fn new(
        id: PhysicalObjectId,
        cochain: RealRelativeCochain,
        divergence_receipt: StableId,
        terminal_constraint_receipt: StableId,
    ) -> Result<Self, TerminalRelativeError> {
        if cochain.units != Current::DIMS {
            return Err(TerminalRelativeError::DistributedCurrentUnits {
                actual: cochain.units,
            });
        }
        if divergence_receipt == terminal_constraint_receipt {
            return Err(TerminalRelativeError::DuplicateIdentity {
                role: "distributed-current constraint receipt",
                id: divergence_receipt.as_str().to_owned(),
            });
        }
        Ok(Self {
            id,
            cochain,
            divergence_receipt,
            terminal_constraint_receipt,
        })
    }

    /// Underlying real cochain.
    #[must_use]
    pub const fn cochain(&self) -> &RealRelativeCochain {
        &self.cochain
    }

    /// Presented receipt for the claimed divergence constraint.
    #[must_use]
    pub const fn divergence_receipt(&self) -> &StableId {
        &self.divergence_receipt
    }

    /// Presented receipt for the claimed terminal closure constraint.
    #[must_use]
    pub const fn terminal_constraint_receipt(&self) -> &StableId {
        &self.terminal_constraint_receipt
    }

    /// Nominal object reference.
    #[must_use]
    pub fn object_ref(&self) -> PhysicalObjectRef {
        PhysicalObjectRef {
            identity: PhysicalObjectIdentity::Declared(self.id.clone()),
            pair: self.cochain.pair,
            phase: self.cochain.phase.clone(),
            kind: PhysicalObjectKind::DistributedCurrent,
        }
    }
}

/// Geometric/manufacturing coil realization reference.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GeometricCoil {
    id: PhysicalObjectId,
    pair: TerminalRelativePairId,
    phase: PhaseId,
    component: ConductorComponentId,
    connectivity_artifact: StableId,
    manufacturing_artifact: StableId,
}

impl GeometricCoil {
    /// Declare a geometric coil realization without claiming that it realizes
    /// any integral class unless a separate conversion map is admitted.
    pub fn try_new(
        id: PhysicalObjectId,
        pair: &TerminalRelativePair,
        phase: PhaseId,
        component: ConductorComponentId,
        connectivity_artifact: StableId,
        manufacturing_artifact: StableId,
    ) -> Result<Self, TerminalRelativeError> {
        let Some(expected_component) = pair.phase_component(&phase) else {
            return Err(TerminalRelativeError::UnknownPhase {
                phase: phase.as_str().to_owned(),
            });
        };
        if !pair.components.iter().any(|entry| entry.id == component) {
            return Err(TerminalRelativeError::UnknownCoilComponent {
                component: component.as_str().to_owned(),
            });
        }
        if expected_component != &component {
            return Err(TerminalRelativeError::CoilPhaseComponentMismatch {
                phase: phase.as_str().to_owned(),
                expected_component: expected_component.as_str().to_owned(),
                actual_component: component.as_str().to_owned(),
            });
        }
        if connectivity_artifact == manufacturing_artifact {
            return Err(TerminalRelativeError::DuplicateIdentity {
                role: "coil realization artifact",
                id: connectivity_artifact.as_str().to_owned(),
            });
        }
        Ok(Self {
            id,
            pair: pair.identity(),
            phase,
            component,
            connectivity_artifact,
            manufacturing_artifact,
        })
    }

    /// Declared conductor component realized by this geometry.
    #[must_use]
    pub const fn component(&self) -> &ConductorComponentId {
        &self.component
    }

    /// Presented connectivity artifact.
    #[must_use]
    pub const fn connectivity_artifact(&self) -> &StableId {
        &self.connectivity_artifact
    }

    /// Presented manufacturing artifact.
    #[must_use]
    pub const fn manufacturing_artifact(&self) -> &StableId {
        &self.manufacturing_artifact
    }

    /// Nominal object reference.
    #[must_use]
    pub fn object_ref(&self) -> PhysicalObjectRef {
        PhysicalObjectRef {
            identity: PhysicalObjectIdentity::Declared(self.id.clone()),
            pair: self.pair,
            phase: self.phase.clone(),
            kind: PhysicalObjectKind::GeometricCoil,
        }
    }
}

/// Nominal physical object sector.  Equal numbers across sectors do not cast.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum PhysicalObjectKind {
    /// Declared integral one-cycle representative; no quotient-class claim.
    IntegralWindingRepresentative,
    /// Real scalar current amplitude.
    RealCurrentAmplitude,
    /// Distributed real current cochain.
    DistributedCurrent,
    /// Geometric/manufacturing coil realization.
    GeometricCoil,
}

/// Identity carried by one nominal conversion endpoint.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum PhysicalObjectIdentity {
    /// Strong coefficient-bearing identity of an integral representative.
    IntegralWindingRepresentative(IntegralWindingRepresentativeId),
    /// Stable declared identity of a real or geometric object.
    Declared(PhysicalObjectId),
}

/// Typed reference retained by a declared conversion map.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PhysicalObjectRef {
    identity: PhysicalObjectIdentity,
    pair: TerminalRelativePairId,
    phase: PhaseId,
    kind: PhysicalObjectKind,
}

impl PhysicalObjectRef {
    /// Exact nominal/strong endpoint identity.
    #[must_use]
    pub const fn identity(&self) -> &PhysicalObjectIdentity {
        &self.identity
    }

    /// Nominal object kind.
    #[must_use]
    pub const fn kind(&self) -> PhysicalObjectKind {
        self.kind
    }

    /// Pair identity.
    #[must_use]
    pub const fn pair_id(&self) -> TerminalRelativePairId {
        self.pair
    }

    /// Phase identity.
    #[must_use]
    pub const fn phase(&self) -> &PhaseId {
        &self.phase
    }
}

/// Explicit map between two otherwise non-convertible physical object sectors.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum DeclaredPhysicalMapKind {
    /// Bind an integral representative to a geometric conductor realization.
    WindingRealization,
    /// Bind one real scalar current normalization to a distributed cochain.
    CurrentRealization,
}

/// Explicit, typed relationship between otherwise non-convertible sectors.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DeclaredPhysicalMap {
    id: ConversionMapId,
    kind: DeclaredPhysicalMapKind,
    source: PhysicalObjectRef,
    target: PhysicalObjectRef,
    map_artifact: StableId,
}

impl DeclaredPhysicalMap {
    /// Bind an exact map artifact.  This records a conversion boundary; it
    /// does not execute the map or prove its physical correctness.
    pub fn try_new(
        id: ConversionMapId,
        kind: DeclaredPhysicalMapKind,
        source: PhysicalObjectRef,
        target: PhysicalObjectRef,
        map_artifact: StableId,
    ) -> Result<Self, TerminalRelativeError> {
        if source.pair != target.pair {
            return Err(TerminalRelativeError::PairIdentityMismatch);
        }
        if source.phase != target.phase {
            return Err(TerminalRelativeError::ConversionPhaseMismatch {
                source: source.phase.as_str().to_owned(),
                target: target.phase.as_str().to_owned(),
            });
        }
        if source.kind == target.kind {
            return Err(TerminalRelativeError::ConversionKindUnchanged { kind: source.kind });
        }
        let expected = match kind {
            DeclaredPhysicalMapKind::WindingRealization => (
                PhysicalObjectKind::IntegralWindingRepresentative,
                PhysicalObjectKind::GeometricCoil,
            ),
            DeclaredPhysicalMapKind::CurrentRealization => (
                PhysicalObjectKind::RealCurrentAmplitude,
                PhysicalObjectKind::DistributedCurrent,
            ),
        };
        if (source.kind, target.kind) != expected {
            return Err(TerminalRelativeError::ConversionKindMismatch {
                map: kind,
                source: source.kind,
                target: target.kind,
            });
        }
        if source.identity == target.identity {
            return Err(TerminalRelativeError::DuplicateIdentity {
                role: "conversion endpoint",
                id: format!("{:?}", source.identity),
            });
        }
        Ok(Self {
            id,
            kind,
            source,
            target,
            map_artifact,
        })
    }

    /// Conversion-map identity.
    #[must_use]
    pub const fn id(&self) -> &ConversionMapId {
        &self.id
    }

    /// Nominal map family; algebraic maps are not conflated with realization
    /// bindings.
    #[must_use]
    pub const fn kind(&self) -> DeclaredPhysicalMapKind {
        self.kind
    }

    /// Source physical object.
    #[must_use]
    pub const fn source(&self) -> &PhysicalObjectRef {
        &self.source
    }

    /// Target physical object.
    #[must_use]
    pub const fn target(&self) -> &PhysicalObjectRef {
        &self.target
    }

    /// Immutable executable/proof artifact reference.
    #[must_use]
    pub const fn map_artifact(&self) -> &StableId {
        &self.map_artifact
    }
}

fn canonical_pair_payload(
    complex: &FiniteCellComplex,
    conductor: &CellularSubcomplex,
    relative: &CellularSubcomplex,
    insulation: &CellularSubcomplex,
    components: &[ConductorComponent],
    terminals: &[PhysicalTerminal],
) -> Result<Vec<u8>, TerminalRelativeError> {
    let mut writer = CanonicalWriter::new();
    writer.u32(TERMINAL_RELATIVE_SCHEMA_VERSION)?;
    writer.u8(complex.dimension)?;
    writer.len(complex.cell_counts.len())?;
    for count in &complex.cell_counts {
        writer.u32(*count)?;
    }
    writer.len(complex.incidences.len())?;
    for incidence in &complex.incidences {
        writer.cell(incidence.lower)?;
        writer.cell(incidence.upper)?;
        writer.u8(incidence.sign.tag())?;
    }
    writer.subcomplex(conductor)?;
    writer.subcomplex(relative)?;
    writer.subcomplex(insulation)?;
    writer.len(components.len())?;
    for component in components {
        writer.text(component.id.as_str())?;
        writer.subcomplex(&component.support)?;
    }
    writer.len(terminals.len())?;
    for terminal in terminals {
        writer.text(terminal.id.as_str())?;
        writer.subcomplex(&terminal.support)?;
        writer.text(terminal.component.as_str())?;
        writer.text(terminal.phase.as_str())?;
        writer.u8(match terminal.role {
            TerminalRole::Driven => 0,
            TerminalRole::ReturnReference => 1,
        })?;
        writer.u8(match terminal.orientation {
            TerminalOrientation::IntoConductor => 0,
            TerminalOrientation::OutOfConductor => 1,
        })?;
        writer.u8(match terminal.coordinate {
            TerminalPortCoordinate::Effort => 0,
            TerminalPortCoordinate::Flow => 1,
        })?;
        writer.port_schema(&terminal.port)?;
        writer.text(terminal.machine.authority_domain().as_str())?;
        writer.u32(terminal.machine.schema_version())?;
        writer.bytes(terminal.machine.graph_digest())?;
        writer.text(terminal.machine.owner().as_str())?;
        writer.text(terminal.machine.port().as_str())?;
        writer.text(terminal.machine.effort_terminal().as_str())?;
        writer.text(terminal.machine.flow_terminal().as_str())?;
        writer.text(terminal.trivialization.id.as_str())?;
        writer.text(terminal.trivialization.port_id.as_str())?;
        writer.u8(match terminal.trivialization.sign {
            OrientationMapSign::Preserve => 0,
            OrientationMapSign::Reverse => 1,
        })?;
        writer.text(terminal.trivialization.voltage_reference.as_str())?;
        writer.text(terminal.trivialization.current_reference.as_str())?;
    }
    Ok(writer.finish())
}

fn canonical_winding_payload(
    chain: &IntegralRelativeChain,
) -> Result<Vec<u8>, TerminalRelativeError> {
    let mut writer = CanonicalWriter::new();
    writer.u32(TERMINAL_RELATIVE_SCHEMA_VERSION)?;
    writer.bytes(chain.pair.as_bytes())?;
    writer.text(chain.phase.as_str())?;
    writer.u8(chain.degree)?;
    writer.len(chain.coefficients.len())?;
    for coefficient in &chain.coefficients {
        writer.i64(*coefficient)?;
    }
    Ok(writer.finish())
}

fn canonical_signed_relabel_payload(
    entries: &[SignedCellRelabelEntry],
) -> Result<Vec<u8>, TerminalRelativeSignedRelabelError> {
    let mut writer = CanonicalWriter::new();
    writer.u32(TERMINAL_RELATIVE_SCHEMA_VERSION)?;
    writer.len(entries.len())?;
    for entry in entries {
        writer.cell(entry.source)?;
        writer.cell(entry.target)?;
        writer.u8(entry.sign.tag())?;
    }
    Ok(writer.finish())
}

fn canonical_component_relabel_payload(
    entries: &[ComponentRelabelEntry],
) -> Result<Vec<u8>, TerminalRelativePhysicalRelabelError> {
    let mut writer = CanonicalWriter::new();
    writer.u32(TERMINAL_RELATIVE_SCHEMA_VERSION)?;
    writer.len(entries.len())?;
    for entry in entries {
        writer.text(entry.source.as_str())?;
        writer.text(entry.target.as_str())?;
    }
    Ok(writer.finish())
}

fn canonical_phase_relabel_payload(
    entries: &[PhaseRelabelEntry],
) -> Result<Vec<u8>, TerminalRelativePhysicalRelabelError> {
    let mut writer = CanonicalWriter::new();
    writer.u32(TERMINAL_RELATIVE_SCHEMA_VERSION)?;
    writer.len(entries.len())?;
    for entry in entries {
        writer.text(entry.source.as_str())?;
        writer.text(entry.target.as_str())?;
        writer.u8(entry.current_sign.tag())?;
    }
    Ok(writer.finish())
}

fn canonical_terminal_relabel_payload(
    entries: &[TerminalRelabelEntry],
) -> Result<Vec<u8>, TerminalRelativePhysicalRelabelError> {
    let mut writer = CanonicalWriter::new();
    writer.u32(TERMINAL_RELATIVE_SCHEMA_VERSION)?;
    writer.len(entries.len())?;
    for entry in entries {
        writer.text(entry.source.as_str())?;
        writer.text(entry.target.as_str())?;
    }
    Ok(writer.finish())
}

struct CanonicalWriter {
    bytes: Vec<u8>,
}

impl CanonicalWriter {
    fn new() -> Self {
        Self { bytes: Vec::new() }
    }

    fn finish(self) -> Vec<u8> {
        self.bytes
    }

    fn reserve(&self, additional: usize) -> Result<(), TerminalRelativeError> {
        let requested = self.bytes.len().checked_add(additional).ok_or(
            TerminalRelativeError::CanonicalBytesExceeded {
                requested: usize::MAX,
                max: MAX_TERMINAL_RELATIVE_CANONICAL_BYTES,
            },
        )?;
        if requested > MAX_TERMINAL_RELATIVE_CANONICAL_BYTES {
            return Err(TerminalRelativeError::CanonicalBytesExceeded {
                requested,
                max: MAX_TERMINAL_RELATIVE_CANONICAL_BYTES,
            });
        }
        Ok(())
    }

    fn raw(&mut self, bytes: &[u8]) -> Result<(), TerminalRelativeError> {
        self.reserve(bytes.len())?;
        self.bytes.extend_from_slice(bytes);
        Ok(())
    }

    fn u8(&mut self, value: u8) -> Result<(), TerminalRelativeError> {
        self.raw(&[value])
    }

    fn u32(&mut self, value: u32) -> Result<(), TerminalRelativeError> {
        self.raw(&value.to_le_bytes())
    }

    fn u64(&mut self, value: u64) -> Result<(), TerminalRelativeError> {
        self.raw(&value.to_le_bytes())
    }

    fn i64(&mut self, value: i64) -> Result<(), TerminalRelativeError> {
        self.raw(&value.to_le_bytes())
    }

    fn len(&mut self, value: usize) -> Result<(), TerminalRelativeError> {
        self.u64(u64::try_from(value).map_err(|_| TerminalRelativeError::LengthOverflow)?)
    }

    fn bytes(&mut self, value: &[u8]) -> Result<(), TerminalRelativeError> {
        self.len(value.len())?;
        self.raw(value)
    }

    fn text(&mut self, value: &str) -> Result<(), TerminalRelativeError> {
        self.bytes(value.as_bytes())
    }

    fn dims(&mut self, dims: Dims) -> Result<(), TerminalRelativeError> {
        for exponent in dims.0 {
            self.u8(exponent.to_le_bytes()[0])?;
        }
        Ok(())
    }

    fn cell(&mut self, cell: CellRef) -> Result<(), TerminalRelativeError> {
        self.u8(cell.degree)?;
        self.u32(cell.ordinal)
    }

    fn subcomplex(&mut self, subcomplex: &CellularSubcomplex) -> Result<(), TerminalRelativeError> {
        self.text(subcomplex.id.as_str())?;
        self.len(subcomplex.cells.len())?;
        for cell in &subcomplex.cells {
            self.cell(*cell)?;
        }
        Ok(())
    }

    fn port_schema(&mut self, port: &PortSchema) -> Result<(), TerminalRelativeError> {
        self.u64(u64::from(port.version()))?;
        self.text(port.id().as_str())?;
        self.u8(port_kind_tag(port.kind()))?;
        self.dims(port.effort_dimensions())?;
        self.dims(port.flow_dimensions())?;
        match port.shape() {
            PortValueShape::Scalar => self.u8(0)?,
            PortValueShape::Vector(components) => {
                self.u8(1)?;
                self.len(components.get())?;
            }
            PortValueShape::Tensor { rows, columns } => {
                self.u8(2)?;
                self.len(rows.get())?;
                self.len(columns.get())?;
            }
            PortValueShape::Field {
                components,
                effort_space,
                flow_space,
            } => {
                self.u8(3)?;
                self.len(components.get())?;
                self.text(effort_space.name())?;
                self.text(flow_space.name())?;
            }
        }
        self.text(port.coordinates().basis().as_str())?;
        self.text(port.coordinates().frame().as_str())?;
        self.u8(port_orientation_tag(port.coordinates().orientation()))?;
        match port.power_pairing() {
            PowerPairing::ScalarProduct => self.u8(0)?,
            PowerPairing::EuclideanDot => self.u8(1)?,
            PowerPairing::FieldDuality {
                measure_dimensions,
                measure_side,
            } => {
                self.u8(2)?;
                self.dims(measure_dimensions)?;
                self.u8(match measure_side {
                    FieldMeasureSide::Effort => 0,
                    FieldMeasureSide::Flow => 1,
                })?;
            }
        }
        self.text(port.timestamp().clock().as_str())?;
        self.u64(port.timestamp().tick())?;
        self.len(port.conservation_roles().len())?;
        for role in port.conservation_roles() {
            self.u8(conservation_role_tag(*role))?;
        }
        Ok(())
    }
}

const fn port_kind_tag(kind: PortKind) -> u8 {
    match kind {
        PortKind::MechanicalForceVelocity => 0,
        PortKind::FluidPressureFlux => 1,
        PortKind::ThermalTemperatureEntropy => 2,
        PortKind::RotationalTorqueAngularVelocity => 3,
        PortKind::ElectricalVoltageCurrent => 4,
        PortKind::MagneticMmfFluxRate => 5,
        PortKind::ChemicalPotentialAmountFlow => 6,
    }
}

const fn port_orientation_tag(orientation: PortOrientation) -> u8 {
    match orientation {
        PortOrientation::OutwardFromOwner => 0,
        PortOrientation::AlongFrame => 1,
        PortOrientation::AgainstFrame => 2,
    }
}

const fn conservation_role_tag(role: ConservationRole) -> u8 {
    match role {
        ConservationRole::Energy => 0,
        ConservationRole::Mass => 1,
        ConservationRole::Amount => 2,
        ConservationRole::LinearMomentum => 3,
        ConservationRole::AngularMomentum => 4,
        ConservationRole::Entropy => 5,
        ConservationRole::ElectricCharge => 6,
    }
}

/// Structured fail-closed admission diagnostics.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum TerminalRelativeError {
    /// A nominal stable ID was malformed.
    InvalidIdentity {
        /// Identity role.
        role: &'static str,
        /// Rejected spelling.
        value: String,
    },
    /// Complex dimension exceeded the physical lane.
    DimensionTooLarge {
        /// Supplied dimension.
        actual: u8,
        /// Maximum dimension.
        max: u8,
    },
    /// Cell-count vector did not cover exactly `0..=dimension`.
    CellCountArity {
        /// Complex dimension.
        dimension: u8,
        /// Supplied vector length.
        actual: usize,
    },
    /// Cell budget was zero or exceeded.
    CellBudgetExceeded {
        /// Supplied total.
        actual: usize,
        /// Maximum total.
        max: usize,
    },
    /// Incidence budget was exceeded.
    IncidenceBudgetExceeded {
        /// Supplied entries.
        actual: usize,
        /// Maximum entries.
        max: usize,
    },
    /// A cell reference was outside its degree extent.
    CellOutOfRange {
        /// Rejected cell.
        cell: CellRef,
    },
    /// Incidence did not join adjacent dimensions.
    InvalidIncidenceDegree {
        /// Lower cell.
        lower: CellRef,
        /// Upper cell.
        upper: CellRef,
    },
    /// Same boundary-matrix coordinate appeared twice.
    DuplicateIncidence {
        /// Lower cell.
        lower: CellRef,
        /// Upper cell.
        upper: CellRef,
    },
    /// Exact `boundary * boundary` produced a nonzero coefficient.
    BoundarySquaredNonzero {
        /// Source cell.
        source: CellRef,
        /// Twice-lowered target cell.
        target: CellRef,
        /// Exact nonzero coefficient.
        value: i128,
    },
    /// Same subcomplex cell appeared twice in caller input.
    DuplicateSubcomplexCell {
        /// Subcomplex identity.
        id: String,
        /// Duplicate cell.
        cell: CellRef,
    },
    /// Declared support omitted a boundary face.
    NotASubcomplex {
        /// Subcomplex identity.
        id: String,
        /// Cell whose boundary is incomplete.
        cell: CellRef,
        /// Missing face.
        missing_boundary: CellRef,
    },
    /// A support required to be nonempty was empty.
    EmptySupport {
        /// Object role.
        object: &'static str,
        /// Object identity.
        id: String,
    },
    /// Component count was empty or exceeded.
    ComponentBudgetExceeded {
        /// Supplied count.
        actual: usize,
        /// Maximum count.
        max: usize,
    },
    /// Terminal count was below two or exceeded.
    TerminalBudgetExceeded {
        /// Supplied count.
        actual: usize,
        /// Maximum count.
        max: usize,
    },
    /// Two semantic objects reused one identity.
    DuplicateIdentity {
        /// Nominal role.
        role: &'static str,
        /// Duplicate ID.
        id: String,
    },
    /// Component support escaped the conductor.
    ComponentOutsideConductor {
        /// Component ID.
        component: String,
        /// Escaped cell.
        cell: CellRef,
    },
    /// A declared component contained no ambient top-dimensional cell.
    ComponentHasNoTopCell {
        /// Component ID.
        component: String,
        /// Required top degree.
        top_degree: u8,
    },
    /// Component support was not exactly the closure of its top cells.
    ComponentSupportNotTopClosure {
        /// Component ID.
        component: String,
        /// First differing cell.
        cell: CellRef,
    },
    /// The explicit quotient subcomplex was not contained in the conductor.
    RelativeOutsideConductor {
        /// Escaped relative cell.
        cell: CellRef,
    },
    /// Insulation was not contained in the explicit quotient subcomplex.
    InsulationOutsideRelativeSubcomplex {
        /// Escaped insulation cell.
        cell: CellRef,
    },
    /// Component supports overlapped.
    OverlappingComponents {
        /// First overlapping cell.
        cell: CellRef,
    },
    /// Components did not exactly partition the conductor.
    ComponentPartitionMismatch {
        /// First mismatched cell.
        cell: CellRef,
    },
    /// Terminal named a component not in the pair.
    UnknownTerminalComponent {
        /// Terminal ID.
        terminal: String,
        /// Missing component ID.
        component: String,
    },
    /// Terminal support escaped its component.
    TerminalOutsideComponent {
        /// Terminal ID.
        terminal: String,
        /// Component ID.
        component: String,
        /// Escaped cell.
        cell: CellRef,
    },
    /// Terminal support was not contained in the explicit quotient subcomplex.
    TerminalOutsideRelativeSubcomplex {
        /// Terminal ID.
        terminal: String,
        /// Escaped cell.
        cell: CellRef,
    },
    /// Terminal support was not a codimension-one patch with its closure.
    TerminalCodimension {
        /// Terminal ID.
        terminal: String,
        /// Ambient dimension.
        ambient_dimension: u8,
    },
    /// Terminal support contained cells outside its codimension-one closure.
    TerminalSupportNotPatchClosure {
        /// Terminal ID.
        terminal: String,
        /// First differing cell.
        cell: CellRef,
    },
    /// A codimension-one terminal cell was not on its component boundary.
    TerminalNotOnComponentBoundary {
        /// Terminal ID.
        terminal: String,
        /// Unsupported terminal cell.
        cell: CellRef,
    },
    /// Terminal and insulation supports overlapped.
    TerminalInsulationOverlap {
        /// Terminal ID.
        terminal: String,
        /// First overlapping cell.
        cell: CellRef,
    },
    /// Two physical terminals overlapped.
    OverlappingTerminals {
        /// First overlapping cell.
        cell: CellRef,
    },
    /// Winding terminal used a non-electrical port kind.
    TerminalRequiresElectricalPort {
        /// Terminal ID.
        terminal: String,
        /// Rejected kind.
        actual: PortKind,
    },
    /// A winding terminal selected effort instead of electrical flow/current.
    TerminalRequiresFlowCoordinate {
        /// Terminal ID.
        terminal: String,
        /// Rejected coordinate.
        actual: TerminalPortCoordinate,
    },
    /// Trivialization named a different port.
    TrivializationPortMismatch {
        /// Terminal ID.
        terminal: String,
        /// Required port ID.
        expected: String,
        /// Supplied port ID.
        actual: String,
    },
    /// Presented MachineGraph domain or version did not match the L6 schema.
    MachineGraphSchemaMismatch {
        /// Required identity domain.
        expected_domain: &'static str,
        /// Required schema version.
        expected_version: u32,
        /// Presented identity domain.
        actual_domain: String,
        /// Presented schema version.
        actual_version: u32,
    },
    /// Presented MachineGraph digest was the all-zero sentinel.
    ZeroMachineGraphDigest,
    /// Presented Machine-IR effort and flow terminal keys aliased.
    MachinePortTerminalAlias {
        /// Aliased terminal key.
        terminal: String,
    },
    /// Presented Machine-IR port key did not match the complete PortSchema ID.
    MachinePortSchemaMismatch {
        /// Physical terminal ID.
        terminal: String,
        /// Required PortSchema ID.
        expected: String,
        /// Presented Machine-IR port ID.
        actual: String,
    },
    /// This first physical lane only admits owner-outward port coordinates.
    UnsupportedPortOrientation {
        /// Physical terminal ID.
        terminal: String,
        /// Rejected orientation.
        actual: PortOrientation,
    },
    /// Port orientation, trivialization sign, and physical direction disagreed.
    TerminalOrientationMismatch {
        /// Physical terminal ID.
        terminal: String,
        /// Port coordinate orientation.
        port_orientation: PortOrientation,
        /// Explicit coordinate-map sign.
        trivialization: OrientationMapSign,
        /// Rejected physical direction.
        actual: TerminalOrientation,
    },
    /// A phase did not declare exactly one driven and one return terminal.
    PhaseTerminalCardinality {
        /// Phase ID.
        phase: String,
        /// Presented terminal count.
        actual: usize,
    },
    /// A phase omitted driven or return/reference semantics.
    MissingPhaseRole {
        /// Phase ID.
        phase: String,
        /// Missing role.
        role: TerminalRole,
    },
    /// A phase repeated a driven or return/reference role.
    DuplicatePhaseRole {
        /// Phase ID.
        phase: String,
    },
    /// Driven and return terminals of one phase named different components.
    PhaseComponentMismatch {
        /// Phase ID.
        phase: String,
        /// Component named by the driven terminal.
        driven_component: String,
        /// Component named by the return terminal.
        return_component: String,
    },
    /// Two electrical phases claimed the same conductor component.
    ComponentPhaseConflict {
        /// Contested component ID.
        component: String,
        /// First phase ID.
        first_phase: String,
        /// Second phase ID.
        second_phase: String,
    },
    /// A conductor component had no driven/return phase binding.
    UnboundConductorComponent {
        /// Component missing a phase binding.
        component: String,
    },
    /// An admitted pair's derived phase/component binding could not be resolved.
    PhaseComponentBindingLost {
        /// Phase ID.
        phase: String,
        /// Missing component ID.
        component: String,
    },
    /// Phase terminals did not include both current directions.
    PhaseOrientationDoesNotClose {
        /// Phase ID.
        phase: String,
    },
    /// Two terminals of one phase disagreed on a shared convention field.
    PhaseConventionMismatch {
        /// Phase ID.
        phase: String,
        /// First mismatching convention field.
        field: &'static str,
    },
    /// Canonical length did not fit the host representation.
    LengthOverflow,
    /// Canonical payload exceeded its explicit envelope.
    CanonicalBytesExceeded {
        /// Requested bytes.
        requested: usize,
        /// Maximum bytes.
        max: usize,
    },
    /// Strong identity encoder refused the payload.
    CanonicalIdentity(CanonicalError),
    /// Chain/cochain named a phase absent from the pair.
    UnknownPhase {
        /// Phase ID.
        phase: String,
    },
    /// Degree exceeded the complex.
    DegreeOutOfRange {
        /// Supplied degree.
        degree: u8,
        /// Complex dimension.
        dimension: u8,
    },
    /// Coefficient count did not match the relative basis.
    CoefficientArity {
        /// Expected coefficients.
        expected: usize,
        /// Supplied coefficients.
        actual: usize,
    },
    /// Chain/cochain belongs to another pair.
    PairIdentityMismatch,
    /// Degree zero has no boundary target.
    NoBoundaryPredecessor,
    /// Top degree has no coboundary target.
    NoCoboundarySuccessor,
    /// Exact integral accumulation overflowed i64 publication.
    CoefficientOverflow,
    /// Integral chain/cochain evaluation used different phases.
    PairingPhaseMismatch {
        /// Cochain phase.
        cochain: String,
        /// Chain phase.
        chain: String,
    },
    /// Integral chain/cochain evaluation used different degrees.
    PairingDegreeMismatch {
        /// Cochain degree.
        cochain: u8,
        /// Chain degree.
        chain: u8,
    },
    /// Exact integral pairing accumulation overflowed i128.
    PairingOverflow,
    /// A real coefficient was NaN or infinite.
    NonFiniteRealCoefficient {
        /// Coefficient index.
        index: usize,
    },
    /// Candidate winding chain had nonzero relative boundary.
    NotARelativeCycle {
        /// Boundary coefficient index.
        index: usize,
        /// Exact nonzero coefficient.
        coefficient: i64,
    },
    /// Current amplitude was NaN or infinite.
    NonFiniteCurrentAmplitude,
    /// Distributed-current cochain did not carry ampere dimensions.
    DistributedCurrentUnits {
        /// Actual dimensions.
        actual: Dims,
    },
    /// Geometric coil named no admitted component.
    UnknownCoilComponent {
        /// Missing component ID.
        component: String,
    },
    /// A geometric coil used a component owned by another phase.
    CoilPhaseComponentMismatch {
        /// Coil phase.
        phase: String,
        /// Component owned by the coil phase.
        expected_component: String,
        /// Component supplied by the coil declaration.
        actual_component: String,
    },
    /// Conversion endpoints used different phases.
    ConversionPhaseMismatch {
        /// Source phase.
        source: String,
        /// Target phase.
        target: String,
    },
    /// A conversion map did not cross nominal sectors.
    ConversionKindUnchanged {
        /// Repeated sector.
        kind: PhysicalObjectKind,
    },
    /// Declared map family did not match its endpoint sectors.
    ConversionKindMismatch {
        /// Declared map family.
        map: DeclaredPhysicalMapKind,
        /// Actual source sector.
        source: PhysicalObjectKind,
        /// Actual target sector.
        target: PhysicalObjectKind,
    },
}

impl From<CanonicalError> for TerminalRelativeError {
    fn from(value: CanonicalError) -> Self {
        Self::CanonicalIdentity(value)
    }
}

impl fmt::Display for TerminalRelativeError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "terminal-relative admission refused: {self:?}")
    }
}

impl core::error::Error for TerminalRelativeError {}

/// Fail-closed diagnostics for admission and use of an explicit signed
/// terminal-relative relabeling.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum TerminalRelativeSignedRelabelError {
    /// Source and target complexes have different top dimensions.
    ComplexDimensionMismatch {
        /// Source top dimension.
        source: u8,
        /// Target top dimension.
        target: u8,
    },
    /// Source and target have different cell counts at one degree.
    CellCountMismatch {
        /// Differing cell degree.
        degree: u8,
        /// Source cell count.
        source: u32,
        /// Target cell count.
        target: u32,
    },
    /// A declaration did not contain exactly one row per ambient source cell.
    EntryCountMismatch {
        /// Required row count.
        expected: usize,
        /// Supplied row count.
        actual: usize,
    },
    /// A declared source cell is outside the source complex.
    SourceCellOutOfRange {
        /// Rejected source cell.
        cell: CellRef,
    },
    /// A declared target cell is outside the target complex.
    TargetCellOutOfRange {
        /// Rejected target cell.
        cell: CellRef,
    },
    /// One map row changed cellular degree.
    CellDegreeMismatch {
        /// Source cell.
        source: CellRef,
        /// Target cell.
        target: CellRef,
    },
    /// Two declaration rows used the same source cell.
    DuplicateSourceCell {
        /// Repeated source cell.
        cell: CellRef,
    },
    /// Two declaration rows used the same target cell.
    DuplicateTargetCell {
        /// Repeated target cell.
        cell: CellRef,
    },
    /// An admitted operation could not find a required source-map row.
    MissingSourceCell {
        /// Missing source cell.
        cell: CellRef,
    },
    /// The signed source boundary did not equal the target boundary at one
    /// mapped matrix coordinate.
    MappedIncidenceMismatch {
        /// Target-space lower cell.
        lower: CellRef,
        /// Target-space upper cell.
        upper: CellRef,
        /// Sign required by the mapped source incidence, or `None` for an
        /// extra target incidence.
        expected: Option<IncidenceSign>,
        /// Actual target sign, or `None` for a missing target incidence.
        actual: Option<IncidenceSign>,
    },
    /// A preserved cellular support changed its semantic subcomplex identity.
    SubcomplexIdentityMismatch {
        /// Support role.
        role: &'static str,
        /// Owning component or terminal identity, when applicable.
        owner: Option<String>,
        /// Source subcomplex identity.
        source: String,
        /// Target subcomplex identity.
        target: String,
    },
    /// A source support did not map exactly onto its target support.
    MappedSupportMismatch {
        /// Support role.
        role: &'static str,
        /// Owning component or terminal identity, when applicable.
        owner: Option<String>,
        /// Target-space witness cell.
        cell: CellRef,
        /// Whether the mapped source support contains the witness.
        expected_mapped: bool,
        /// Whether the declared target support contains the witness.
        actual_target: bool,
    },
    /// Component or terminal semantic identity sets differ.
    SemanticIdentitySetMismatch {
        /// Identity role.
        role: &'static str,
        /// Differing identity.
        id: String,
        /// Whether the source contains it.
        source_present: bool,
        /// Whether the target contains it.
        target_present: bool,
    },
    /// A phase identity or its component binding changed.
    PhaseComponentBindingMismatch {
        /// Differing phase identity.
        phase: String,
        /// Source component, if the phase exists there.
        source_component: Option<String>,
        /// Target component, if the phase exists there.
        target_component: Option<String>,
    },
    /// Non-support terminal semantics changed under the relabeling.
    TerminalMetadataMismatch {
        /// Terminal identity.
        terminal: String,
        /// First differing metadata field.
        field: &'static str,
    },
    /// A pair or transported value was not bound to the expected endpoint.
    PairIdentityMismatch {
        /// Endpoint or value role.
        role: &'static str,
        /// Expected strong pair identity.
        expected: TerminalRelativePairId,
        /// Actual strong pair identity.
        actual: TerminalRelativePairId,
    },
    /// A phase-local source basis cell mapped outside the corresponding target
    /// basis.
    MappedBasisCellMissing {
        /// Preserved phase identity.
        phase: String,
        /// Basis degree.
        degree: u8,
        /// Source basis cell.
        source: CellRef,
        /// Mapped target cell.
        target: CellRef,
    },
    /// Negating an exact integral coefficient overflowed `i64`.
    CoefficientOverflow {
        /// Source cell whose coefficient overflowed.
        cell: CellRef,
    },
    /// Existing terminal-relative construction or validation refused.
    TerminalRelative(TerminalRelativeError),
    /// Strong canonical identity admission refused.
    CanonicalIdentity(CanonicalError),
}

impl From<TerminalRelativeError> for TerminalRelativeSignedRelabelError {
    fn from(value: TerminalRelativeError) -> Self {
        Self::TerminalRelative(value)
    }
}

impl From<CanonicalError> for TerminalRelativeSignedRelabelError {
    fn from(value: CanonicalError) -> Self {
        Self::CanonicalIdentity(value)
    }
}

impl fmt::Display for TerminalRelativeSignedRelabelError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "terminal-relative signed relabel refused: {self:?}")
    }
}

impl core::error::Error for TerminalRelativeSignedRelabelError {}

/// Fail-closed diagnostics for admission and use of an explicit physical
/// terminal-relative relabeling.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum TerminalRelativePhysicalRelabelError {
    /// The complete signed cell map or exact top-level support checks refused.
    CellRelabel(TerminalRelativeSignedRelabelError),
    /// Source and target semantic identity sets have different cardinality.
    SemanticCardinalityMismatch {
        /// Semantic identity role.
        role: &'static str,
        /// Source identity count.
        source: usize,
        /// Target identity count.
        target: usize,
    },
    /// A semantic map did not declare exactly one row per source identity.
    SemanticEntryCountMismatch {
        /// Semantic identity role.
        role: &'static str,
        /// Required row count.
        expected: usize,
        /// Supplied row count.
        actual: usize,
    },
    /// A semantic row names an identity absent from the source pair.
    UnknownSemanticSource {
        /// Semantic identity role.
        role: &'static str,
        /// Unknown source identity.
        id: String,
    },
    /// A semantic row names an identity absent from the target pair.
    UnknownSemanticTarget {
        /// Semantic identity role.
        role: &'static str,
        /// Unknown target identity.
        id: String,
    },
    /// Two semantic rows use the same source identity.
    DuplicateSemanticSource {
        /// Semantic identity role.
        role: &'static str,
        /// Repeated source identity.
        id: String,
    },
    /// Two semantic rows use the same target identity.
    DuplicateSemanticTarget {
        /// Semantic identity role.
        role: &'static str,
        /// Repeated target identity.
        id: String,
    },
    /// Explicitly mapped component or terminal support cells differ.
    MappedSemanticSupportMismatch {
        /// Mapped support role.
        role: &'static str,
        /// Source semantic owner.
        source_owner: String,
        /// Target semantic owner.
        target_owner: String,
        /// Target-space witness cell.
        cell: CellRef,
        /// Whether the mapped source support contains the witness.
        expected_mapped: bool,
        /// Whether the target support contains the witness.
        actual_target: bool,
    },
    /// Phase and component permutations do not form a commuting square.
    PhaseComponentSquareMismatch {
        /// Source phase identity.
        source_phase: String,
        /// Declared target phase identity.
        target_phase: String,
        /// Component required by the component permutation.
        expected_target_component: String,
        /// Component actually bound to the target phase.
        actual_target_component: String,
    },
    /// Terminal and phase permutations do not form a commuting square.
    TerminalPhaseSquareMismatch {
        /// Source terminal identity.
        source_terminal: String,
        /// Declared target terminal identity.
        target_terminal: String,
        /// Phase required by the phase permutation.
        expected_target_phase: String,
        /// Phase actually bound to the target terminal.
        actual_target_phase: String,
    },
    /// Terminal and component permutations do not form a commuting square.
    TerminalComponentSquareMismatch {
        /// Source terminal identity.
        source_terminal: String,
        /// Declared target terminal identity.
        target_terminal: String,
        /// Component required by the component permutation.
        expected_target_component: String,
        /// Component actually bound to the target terminal.
        actual_target_component: String,
    },
    /// Terminal parity or non-nominal port/Machine convention changed.
    TerminalConventionMismatch {
        /// Source terminal identity.
        source_terminal: String,
        /// Declared target terminal identity.
        target_terminal: String,
        /// First incompatible field or parity rule.
        field: &'static str,
    },
    /// A target distributed-current constraint receipt reused source evidence.
    ConstraintReceiptNotFresh {
        /// Target constraint role being supplied.
        role: &'static str,
        /// Reused source receipt identity.
        receipt: String,
    },
    /// A target geometric-coil artifact reused source realization evidence.
    GeometricCoilArtifactNotFresh {
        /// Target artifact role being supplied.
        role: &'static str,
        /// Reused source artifact identity.
        artifact: String,
    },
    /// A pair or transported value was not bound to the expected endpoint.
    PairIdentityMismatch {
        /// Endpoint or value role.
        role: &'static str,
        /// Expected pair identity.
        expected: TerminalRelativePairId,
        /// Actual pair identity.
        actual: TerminalRelativePairId,
    },
    /// An admitted or composed operation could not find a required image.
    MissingSemanticImage {
        /// Image role.
        role: &'static str,
        /// Missing source identity or cell coordinate.
        id: String,
    },
    /// A mapped phase-local basis cell was absent from the target basis.
    MappedBasisCellMissing {
        /// Source phase identity.
        source_phase: String,
        /// Target phase identity.
        target_phase: String,
        /// Basis degree.
        degree: u8,
        /// Source basis cell.
        source: CellRef,
        /// Mapped target cell.
        target: CellRef,
    },
    /// The single required exact coefficient negation overflowed `i64`.
    CoefficientOverflow {
        /// Source cell whose coefficient overflowed.
        cell: CellRef,
    },
    /// Existing terminal-relative construction or validation refused.
    TerminalRelative(TerminalRelativeError),
    /// Strong canonical identity admission refused.
    CanonicalIdentity(CanonicalError),
}

impl From<TerminalRelativeError> for TerminalRelativePhysicalRelabelError {
    fn from(value: TerminalRelativeError) -> Self {
        Self::TerminalRelative(value)
    }
}

impl From<CanonicalError> for TerminalRelativePhysicalRelabelError {
    fn from(value: CanonicalError) -> Self {
        Self::CanonicalIdentity(value)
    }
}

impl fmt::Display for TerminalRelativePhysicalRelabelError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "terminal-relative physical relabel refused: {self:?}")
    }
}

impl core::error::Error for TerminalRelativePhysicalRelabelError {}
