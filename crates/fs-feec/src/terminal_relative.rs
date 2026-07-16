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
    CanonicalEncoder, CanonicalError, CanonicalLimits, CanonicalSchema, Field, FieldSpec,
    IdentityReceipt, NeverCancel, SemanticId, StrongIdentity, WireType,
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
