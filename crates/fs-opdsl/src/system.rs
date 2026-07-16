//! Multi-field equation/type IR with a stable semantic identity
//! (bead i94v.1.1.1). This is the neutral compiler-facing layer above
//! the single-field [`crate::expr`] machinery: block fields declare
//! their differential-form degree, basis/frame/orientation, quantity
//! kind, state ownership, clock, spatial support, and conservation
//! roles EXPLICITLY, and every cross-field composition is admitted
//! before any lowering. Ill-typed contractions, mixed frames, affine
//! temperature misuse, clock mismatches, and non-power-conjugate port
//! pairings are structured refusals — they cannot reach code
//! generation.
//!
//! Identity discipline: a [`SystemDef`]'s [`SystemId`] hashes the
//! canonical STRUCTURE and semantic payloads only. Field and equation
//! display names are never hash inputs, and both tables are re-sorted
//! into a canonical order before encoding, so equivalent renaming or
//! serialization order preserves identity while any meaningful
//! convention change (degree, dims, quantity kind, frame, orientation,
//! clock, support, state slot, pairing, scalar convention) moves it.
//! Two fields whose complete semantic payloads are byte-identical are
//! refused as ambiguous rather than silently tie-broken.

use crate::expr::Space;
use fs_blake3::identity::{
    CanonicalEncoder, CanonicalLimits, CanonicalSchema, Field, FieldSpec, NeverCancel, SemanticId,
    WireType,
};
use fs_couple::{PortKind, PortOrientation};
use fs_qty::Dims;
use fs_qty::semantic::SemanticType;

/// IR language version for the multi-field system surface. Bound into
/// every [`SystemId`]; parsers/consumers refuse other versions unless
/// an explicit audited migration runs first.
pub const SYSTEM_IR_VERSION: u32 = 1;

/// Depth cap for [`SystemExpr`] trees: adversarial nesting refuses
/// structurally instead of exhausting the stack (all traversals here
/// are explicit-stack iterative; the cap bounds work, not recursion).
pub const MAX_SYSTEM_EXPR_DEPTH: usize = 512;

/// Maximum fields/equations per system: a bounded, auditable IR.
pub const MAX_SYSTEM_FIELDS: usize = 256;

/// Maximum opaque-extension bytes retained (and hashed) per system.
pub const MAX_SYSTEM_EXTENSION_BYTES: usize = 4096;

const IDENTITY_LIMITS: CanonicalLimits = CanonicalLimits::new(1 << 20, 1 << 18, 64, 16_384, 4096);
const POWER_DIMS: Dims = Dims([2, 1, -3, 0, 0, 0]);

/// A validated machine identifier for frames, bases, charts, and
/// clocks: 1..=64 ASCII graphic bytes. These are REFERENCES into
/// caller-owned registries — the IR binds which frame/clock a field
/// uses, never what the frame/clock is.
#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct ConventionRef(String);

impl ConventionRef {
    /// Validate and wrap a convention reference.
    ///
    /// # Errors
    /// [`SystemTypeError::InvalidConventionRef`] for empty, oversized,
    /// or non-graphic-ASCII input.
    pub fn new(raw: impl Into<String>) -> Result<Self, SystemTypeError> {
        let raw = raw.into();
        if raw.is_empty() || raw.len() > 64 || !raw.bytes().all(|b| b.is_ascii_graphic()) {
            return Err(SystemTypeError::InvalidConventionRef { raw });
        }
        Ok(Self(raw))
    }

    /// The validated reference text.
    #[must_use]
    pub fn as_str(&self) -> &str {
        &self.0
    }
}

/// What a field's stored values mean, beyond raw dimensions.
/// Dimensional equality is not semantic equality (fs-qty semantic
/// kinds): torque is not energy, absolute temperature is not a
/// temperature difference.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FieldQuantity {
    /// Plain dimensional values with no further semantic claim.
    Dimensional(Dims),
    /// A semantic quantity kind + value form; its expected dims are
    /// authoritative and must agree with the field's [`Space`].
    Semantic(SemanticType),
}

impl FieldQuantity {
    /// The dimensions this quantity carries.
    #[must_use]
    pub fn dims(&self) -> Dims {
        match self {
            FieldQuantity::Dimensional(dims) => *dims,
            FieldQuantity::Semantic(semantic) => semantic.expected_dims(),
        }
    }

    /// Affine absolute quantities (absolute temperature) admit neither
    /// scaling nor summation — only differences move them.
    #[must_use]
    pub fn is_affine_absolute(&self) -> bool {
        matches!(
            self,
            FieldQuantity::Semantic(semantic)
                if semantic.kind() == fs_qty::semantic::QuantityKind::AbsoluteTemperature
        )
    }
}

/// Where a field lives on the domain.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SpatialSupport {
    /// Values over the interior of the named chart's domain.
    Interior,
    /// A boundary trace (the k-form restricted to the boundary).
    BoundaryTrace,
}

/// Who owns a field's storage across time steps.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StateOwnership {
    /// Evolving state owned by this system in a distinct slot.
    Owned {
        /// The state-slot index; distinct per owned field.
        slot: u32,
    },
    /// Read-only values supplied per evaluation by the caller.
    External,
    /// A parameter (design/material/control) with a declared role tag.
    Parameter {
        /// Which parameter table this field indexes.
        role: ParameterRole,
    },
}

/// The role a parameter field plays (identity-bearing).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ParameterRole {
    /// Design/geometry parameters.
    Design,
    /// Material/constitutive parameters.
    Material,
    /// Control/actuation parameters.
    Control,
}

/// The coordinate convention a field's components are expressed in:
/// which basis and frame (as references) and the orientation sense.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct CoordinateConvention {
    /// Basis reference (e.g. a chart's coordinate basis id).
    pub basis: ConventionRef,
    /// Frame reference (observer/material frame id).
    pub frame: ConventionRef,
    /// Orientation sense, reusing the port-schema vocabulary so port
    /// pairings and field conventions cannot drift apart.
    pub orientation: PortOrientation,
}

/// Scalar convention for the whole system (identity-bearing).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScalarConvention {
    /// Real arithmetic only.
    RealOnly,
    /// Complex fields under the Hermitian inner-product convention
    /// (conjugate-linear in the FIRST argument).
    ComplexHermitian,
}

/// One declared block field. The `name` is DISPLAY ONLY — never a hash
/// input; every other member is identity-bearing.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FieldDecl {
    /// Human-facing display name (diagnostics only).
    pub name: String,
    /// The field's discrete space: cochain degree, dof count, dims.
    pub space: Space,
    /// What the values mean (dims or a semantic quantity kind).
    pub quantity: FieldQuantity,
    /// Basis/frame/orientation convention references.
    pub coordinates: CoordinateConvention,
    /// The time clock this field advances on (reference).
    pub clock: ConventionRef,
    /// Interior or boundary-trace support.
    pub support: SpatialSupport,
    /// State ownership.
    pub state: StateOwnership,
}

/// A checked handle into a [`SystemDef`]'s field table.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FieldId(pub usize);

/// The multi-field expression tree. Cross-field structure only —
/// single-field operator chains keep using [`crate::expr::Expr`]
/// inside each block via atoms.
#[derive(Debug, Clone, PartialEq)]
pub enum SystemExpr {
    /// A declared field's current value.
    FieldRef(FieldId),
    /// Apply a registered single-field atom (by index into the
    /// system's atom space table) to a sub-expression.
    Apply {
        /// Index into [`SystemDef`]'s atom-space table.
        atom: usize,
        /// The argument.
        arg: Box<SystemExpr>,
    },
    /// Scale by a real constant (refused on affine-absolute operands).
    Scale(f64, Box<SystemExpr>),
    /// Sum of two expressions in the same admitted space/convention.
    Add(Box<SystemExpr>, Box<SystemExpr>),
    /// A power-conjugate pairing of two sub-expressions (effort, flow)
    /// under a declared port kind: admitted only when the operand
    /// dimensions multiply (with the measure) to power.
    PortPair {
        /// The physical port family this pairing claims.
        kind: PortKind,
        /// Effort side.
        effort: Box<SystemExpr>,
        /// Flow side.
        flow: Box<SystemExpr>,
        /// Integration-measure dimensions applied by the pairing
        /// (e.g. volume for field duality; NONE for scalar ports).
        measure_dims: Dims,
    },
}

/// An atom signature registered on the system: the in/out spaces plus
/// the coordinate convention and clock it preserves. (The numerical
/// atom itself stays in the single-field layer; the system IR only
/// needs its type.)
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct AtomSignature {
    /// Display name (diagnostics only; not a hash input).
    pub name: String,
    /// Identity-bearing content reference into the caller's operator
    /// registry (e.g. `feec:d0:kuhn2-v1`). Two atoms with equal spaces
    /// are still DIFFERENT operators; this reference is what
    /// distinguishes them in the semantic identity, exactly like frame
    /// and clock references.
    pub content: ConventionRef,
    /// Input space.
    pub in_space: Space,
    /// Output space.
    pub out_space: Space,
}

/// One block equation: `residual(target) = rhs`, i.e. the rhs
/// expression contributes to the target field's residual block. The
/// rhs space must exactly match the target field's space, convention,
/// and clock.
#[derive(Debug, Clone, PartialEq)]
pub struct BlockEquation {
    /// Display name (diagnostics only; not a hash input).
    pub name: String,
    /// The field whose residual block this equation defines.
    pub target: FieldId,
    /// The right-hand side.
    pub rhs: SystemExpr,
}

/// Structured refusals: every admissibility failure names what failed
/// and both sides, with a teaching hint. These fire BEFORE lowering.
#[derive(Debug, Clone, PartialEq)]
pub enum SystemTypeError {
    /// A convention reference was empty/oversized/non-ASCII-graphic.
    InvalidConventionRef {
        /// The offending raw text.
        raw: String,
    },
    /// A field/atom/equation id out of range.
    UnknownId {
        /// Which table.
        what: &'static str,
        /// The offending index.
        id: usize,
    },
    /// A semantic quantity kind whose expected dims disagree with the
    /// field's space dims.
    QuantityDimsMismatch {
        /// Field display name.
        field: String,
        /// Dims declared on the space.
        space_dims: Dims,
        /// Dims the semantic kind requires.
        kind_dims: Dims,
    },
    /// Two owned fields share a state slot.
    DuplicateStateSlot {
        /// The shared slot.
        slot: u32,
        /// First field display name.
        first: String,
        /// Second field display name.
        second: String,
    },
    /// Cross-field composition across different clocks.
    ClockMismatch {
        /// Left clock reference.
        left: String,
        /// Right clock reference.
        right: String,
    },
    /// Cross-field composition across different frames/bases or
    /// non-composable orientations.
    ConventionMismatch {
        /// Left (basis, frame) references.
        left: (String, String),
        /// Right (basis, frame) references.
        right: (String, String),
    },
    /// Sum/apply across mismatched spaces (degree, dof count, dims).
    SpaceMismatch {
        /// What was being composed.
        context: &'static str,
        /// Expected space.
        expected: Space,
        /// Found space.
        found: Space,
    },
    /// Scaling or summing an affine absolute quantity (temperature):
    /// only differences move affine quantities.
    AffineQuantityMisuse {
        /// The operation attempted.
        operation: &'static str,
        /// The field display name carrying the affine kind.
        field: String,
    },
    /// A port pairing whose effort x flow (x measure) is not power.
    NonConjugatePairing {
        /// The claimed port family.
        kind: PortKind,
        /// Effort dims found.
        effort_dims: Dims,
        /// Flow dims found.
        flow_dims: Dims,
        /// Measure dims declared.
        measure_dims: Dims,
    },
    /// Expression nesting beyond [`MAX_SYSTEM_EXPR_DEPTH`].
    DepthExceeded {
        /// The cap that was crossed.
        cap: usize,
    },
    /// Table growth beyond [`MAX_SYSTEM_FIELDS`].
    TooManyFields {
        /// The cap.
        cap: usize,
    },
    /// Two fields whose complete semantic payloads are byte-identical:
    /// identity cannot order them canonically, so the system is
    /// ambiguous — distinguish them by slot, role, clock, or kind.
    IndistinguishableFields {
        /// First field display name.
        first: String,
        /// Second field display name.
        second: String,
    },
    /// Duplicate equations (identical target + canonical rhs).
    DuplicateEquation {
        /// The equation display name of the later duplicate.
        name: String,
    },
    /// Opaque extension payload beyond its byte bound.
    ExtensionTooLarge {
        /// Bytes supplied.
        len: usize,
        /// The cap.
        cap: usize,
    },
    /// The canonical identity encoder refused (resource bounds).
    IdentityEncoding {
        /// The encoder's message.
        detail: String,
    },
    /// A versioned payload from another IR version: run an explicit
    /// audited migration first.
    VersionMismatch {
        /// Version found.
        found: u32,
        /// Version this build reads/writes.
        supported: u32,
    },
    /// A non-finite scale constant.
    NonFiniteScale {
        /// The offending bits.
        bits: u64,
    },
    /// Two atoms whose complete semantic payloads (content reference
    /// plus spaces) are byte-identical: declare one atom once.
    IndistinguishableAtoms {
        /// First atom display name.
        first: String,
        /// Second atom display name.
        second: String,
    },
    /// A display name the transport cannot carry fail-closed.
    InvalidName {
        /// The offending name.
        name: String,
    },
    /// A transport line the strict parser refuses.
    TransportMalformed {
        /// 1-based line number.
        line: usize,
        /// What was wrong.
        detail: String,
    },
}

impl std::fmt::Display for SystemTypeError {
    fn fmt(&self, out: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::InvalidConventionRef { raw } => write!(
                out,
                "convention reference {raw:?} must be 1..=64 ASCII graphic bytes"
            ),
            Self::UnknownId { what, id } => write!(out, "unknown {what} id {id}"),
            Self::QuantityDimsMismatch {
                field,
                space_dims,
                kind_dims,
            } => write!(
                out,
                "field {field:?}: space dims {space_dims:?} disagree with the semantic kind's expected dims {kind_dims:?}"
            ),
            Self::DuplicateStateSlot {
                slot,
                first,
                second,
            } => write!(
                out,
                "fields {first:?} and {second:?} both own state slot {slot}"
            ),
            Self::ClockMismatch { left, right } => write!(
                out,
                "cross-field composition across clocks {left:?} and {right:?}: interpose an explicit clock-transfer operator"
            ),
            Self::ConventionMismatch { left, right } => write!(
                out,
                "cross-field composition across conventions (basis,frame) {left:?} vs {right:?}: interpose an explicit pullback"
            ),
            Self::SpaceMismatch {
                context,
                expected,
                found,
            } => write!(
                out,
                "{context}: expected space {expected:?}, found {found:?}"
            ),
            Self::AffineQuantityMisuse { operation, field } => write!(
                out,
                "{operation} on affine absolute quantity of field {field:?}: only differences move affine quantities"
            ),
            Self::NonConjugatePairing {
                kind,
                effort_dims,
                flow_dims,
                measure_dims,
            } => write!(
                out,
                "{kind:?} pairing is not power-conjugate: effort {effort_dims:?} x flow {flow_dims:?} x measure {measure_dims:?} != power"
            ),
            Self::DepthExceeded { cap } => write!(out, "expression depth exceeds the {cap} cap"),
            Self::TooManyFields { cap } => write!(out, "system exceeds the {cap}-field cap"),
            Self::IndistinguishableFields { first, second } => write!(
                out,
                "fields {first:?} and {second:?} have byte-identical semantic payloads: identity is ambiguous — distinguish by slot, role, clock, or kind"
            ),
            Self::DuplicateEquation { name } => {
                write!(out, "equation {name:?} duplicates an earlier one")
            }
            Self::ExtensionTooLarge { len, cap } => {
                write!(out, "opaque extension is {len} bytes; the cap is {cap}")
            }
            Self::IdentityEncoding { detail } => write!(out, "identity encoding refused: {detail}"),
            Self::VersionMismatch { found, supported } => write!(
                out,
                "system IR version {found} is not the supported {supported}: run an explicit audited migration"
            ),
            Self::NonFiniteScale { bits } => {
                write!(out, "scale constant is not finite (bits {bits:#018x})")
            }
            Self::IndistinguishableAtoms { first, second } => write!(
                out,
                "atoms {first:?} and {second:?} have byte-identical content/space payloads: declare one atom once"
            ),
            Self::InvalidName { name } => write!(
                out,
                "display name {name:?} must be non-empty printable text without tabs or newlines"
            ),
            Self::TransportMalformed { line, detail } => {
                write!(out, "transport line {line}: {detail}")
            }
        }
    }
}

impl std::error::Error for SystemTypeError {}

/// Identity schema for one admitted multi-field system.
pub struct SystemIdentitySchemaV1;

impl CanonicalSchema for SystemIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-opdsl.system.v1";
    const NAME: &'static str = "multi-field-system";
    const VERSION: u32 = 1;
    const CONTEXT: &'static str = "canonical multi-field equation/type IR: fields and equations in canonical order, display names excluded";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("system-ir-version", WireType::U64),
        FieldSpec::required("scalar-convention", WireType::Variant),
        FieldSpec::required("atom-signatures", WireType::OrderedBytes),
        FieldSpec::required("fields", WireType::OrderedBytes),
        FieldSpec::required("equations", WireType::OrderedBytes),
        FieldSpec::optional_bytes("opaque-extension"),
    ];
}

/// The stable semantic identity of an admitted system.
pub type SystemId = SemanticId<SystemIdentitySchemaV1>;

/// A validated multi-field system: fields, atom signatures, and block
/// equations, all admitted, plus the canonical identity.
#[derive(Debug)]
pub struct AdmittedSystem {
    fields: Vec<FieldDecl>,
    atoms: Vec<AtomSignature>,
    equations: Vec<BlockEquation>,
    scalar_convention: ScalarConvention,
    extension: Vec<u8>,
    identity: SystemId,
}

impl AdmittedSystem {
    /// The declared fields, in declaration order.
    #[must_use]
    pub fn fields(&self) -> &[FieldDecl] {
        &self.fields
    }

    /// The registered atom signatures.
    #[must_use]
    pub fn atoms(&self) -> &[AtomSignature] {
        &self.atoms
    }

    /// The block equations, in declaration order.
    #[must_use]
    pub fn equations(&self) -> &[BlockEquation] {
        &self.equations
    }

    /// The system's scalar convention.
    #[must_use]
    pub fn scalar_convention(&self) -> ScalarConvention {
        self.scalar_convention
    }

    /// The retained opaque extension bytes (identity-bearing).
    #[must_use]
    pub fn extension(&self) -> &[u8] {
        &self.extension
    }

    /// The stable semantic identity.
    #[must_use]
    pub fn identity(&self) -> SystemId {
        self.identity
    }
}

/// Builder for a multi-field system. Declaration order is preserved
/// for display; identity is order- and name-independent.
#[derive(Debug, Default)]
pub struct SystemDef {
    fields: Vec<FieldDecl>,
    atoms: Vec<AtomSignature>,
    equations: Vec<BlockEquation>,
    scalar_convention: Option<ScalarConvention>,
    extension: Vec<u8>,
}

impl SystemDef {
    /// An empty system under [`ScalarConvention::RealOnly`].
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Set the scalar convention (default [`ScalarConvention::RealOnly`]).
    #[must_use]
    pub fn scalar_convention(mut self, convention: ScalarConvention) -> Self {
        self.scalar_convention = Some(convention);
        self
    }

    /// Attach opaque extension bytes: retained verbatim and
    /// identity-bearing, so future dialects can extend without
    /// silently aliasing v1 identities.
    ///
    /// # Errors
    /// [`SystemTypeError::ExtensionTooLarge`] beyond the byte cap.
    pub fn with_extension(mut self, bytes: Vec<u8>) -> Result<Self, SystemTypeError> {
        if bytes.len() > MAX_SYSTEM_EXTENSION_BYTES {
            return Err(SystemTypeError::ExtensionTooLarge {
                len: bytes.len(),
                cap: MAX_SYSTEM_EXTENSION_BYTES,
            });
        }
        self.extension = bytes;
        Ok(self)
    }

    /// Declare a field.
    ///
    /// # Errors
    /// Refuses table overflow, semantic-kind/space dims disagreement,
    /// and duplicate owned state slots.
    pub fn declare_field(&mut self, field: FieldDecl) -> Result<FieldId, SystemTypeError> {
        if self.fields.len() >= MAX_SYSTEM_FIELDS {
            return Err(SystemTypeError::TooManyFields {
                cap: MAX_SYSTEM_FIELDS,
            });
        }
        if let FieldQuantity::Semantic(semantic) = &field.quantity {
            let kind_dims = semantic.expected_dims();
            if kind_dims != field.space.dims {
                return Err(SystemTypeError::QuantityDimsMismatch {
                    field: field.name.clone(),
                    space_dims: field.space.dims,
                    kind_dims,
                });
            }
        }
        if let StateOwnership::Owned { slot } = field.state
            && let Some(previous) = self.fields.iter().find(
                |existing| matches!(existing.state, StateOwnership::Owned { slot: s } if s == slot),
            )
        {
            return Err(SystemTypeError::DuplicateStateSlot {
                slot,
                first: previous.name.clone(),
                second: field.name.clone(),
            });
        }
        self.fields.push(field);
        Ok(FieldId(self.fields.len() - 1))
    }

    /// Register an atom signature.
    pub fn register_atom(&mut self, atom: AtomSignature) -> usize {
        self.atoms.push(atom);
        self.atoms.len() - 1
    }

    /// Add a block equation after full admissibility checking of its
    /// right-hand side against the target field.
    ///
    /// # Errors
    /// Any [`SystemTypeError`] the rhs or target binding produces.
    pub fn add_equation(&mut self, equation: BlockEquation) -> Result<(), SystemTypeError> {
        let target = self
            .fields
            .get(equation.target.0)
            .ok_or(SystemTypeError::UnknownId {
                what: "field",
                id: equation.target.0,
            })?;
        let admitted = self.admit_expr(&equation.rhs)?;
        if admitted.space != target.space {
            return Err(SystemTypeError::SpaceMismatch {
                context: "equation rhs vs target field",
                expected: target.space,
                found: admitted.space,
            });
        }
        if let Some(convention) = &admitted.convention {
            check_convention_match(convention, target)?;
        }
        self.equations.push(equation);
        Ok(())
    }

    /// Admit the whole system and mint its canonical identity.
    ///
    /// # Errors
    /// Any residual admissibility refusal (indistinguishable fields,
    /// duplicate equations, encoder bounds).
    pub fn admit(self) -> Result<AdmittedSystem, SystemTypeError> {
        let scalar_convention = self.scalar_convention.unwrap_or(ScalarConvention::RealOnly);

        // Canonical field payloads (names excluded).
        let mut field_payloads: Vec<(Vec<u8>, usize)> = self
            .fields
            .iter()
            .enumerate()
            .map(|(index, field)| (canonical_field_bytes(field), index))
            .collect();
        field_payloads.sort();
        for pair in field_payloads.windows(2) {
            if pair[0].0 == pair[1].0 {
                return Err(SystemTypeError::IndistinguishableFields {
                    first: self.fields[pair[0].1].name.clone(),
                    second: self.fields[pair[1].1].name.clone(),
                });
            }
        }
        // Declaration ordinal -> canonical ordinal.
        let mut remap = vec![0usize; self.fields.len()];
        for (canonical, (_, declaration)) in field_payloads.iter().enumerate() {
            remap[*declaration] = canonical;
        }

        let mut atom_payloads: Vec<(Vec<u8>, usize)> = self
            .atoms
            .iter()
            .enumerate()
            .map(|(index, atom)| (canonical_atom_bytes(atom), index))
            .collect();
        atom_payloads.sort();
        for pair in atom_payloads.windows(2) {
            if pair[0].0 == pair[1].0 {
                return Err(SystemTypeError::IndistinguishableAtoms {
                    first: self.atoms[pair[0].1].name.clone(),
                    second: self.atoms[pair[1].1].name.clone(),
                });
            }
        }
        let mut atom_remap = vec![0usize; self.atoms.len()];
        for (canonical, (_, declaration)) in atom_payloads.iter().enumerate() {
            atom_remap[*declaration] = canonical;
        }

        let mut equation_payloads: Vec<Vec<u8>> = self
            .equations
            .iter()
            .map(|equation| canonical_equation_bytes(equation, &remap, &atom_remap))
            .collect();
        equation_payloads.sort();
        for pair in equation_payloads.windows(2) {
            if pair[0] == pair[1] {
                let duplicate = self
                    .equations
                    .last()
                    .map(|equation| equation.name.clone())
                    .unwrap_or_default();
                return Err(SystemTypeError::DuplicateEquation { name: duplicate });
            }
        }

        let encode = || -> Result<_, fs_blake3::identity::CanonicalError> {
            let mut encoder = CanonicalEncoder::<SystemId, _>::new(IDENTITY_LIMITS, NeverCancel)?;
            encoder = encoder.u64(
                Field::new(0, "system-ir-version"),
                u64::from(SYSTEM_IR_VERSION),
            )?;
            encoder = encoder.variant(
                Field::new(1, "scalar-convention"),
                match scalar_convention {
                    ScalarConvention::RealOnly => 0,
                    ScalarConvention::ComplexHermitian => 1,
                },
                &[],
            )?;
            encoder = encoder.ordered_bytes(
                Field::new(2, "atom-signatures"),
                atom_payloads.len() as u64,
                atom_payloads.iter().map(|(bytes, _)| bytes.as_slice()),
            )?;
            encoder = encoder.ordered_bytes(
                Field::new(3, "fields"),
                field_payloads.len() as u64,
                field_payloads.iter().map(|(bytes, _)| bytes.as_slice()),
            )?;
            encoder = encoder.ordered_bytes(
                Field::new(4, "equations"),
                equation_payloads.len() as u64,
                equation_payloads.iter().map(Vec::as_slice),
            )?;
            encoder = encoder.optional_bytes(
                Field::new(5, "opaque-extension"),
                if self.extension.is_empty() {
                    None
                } else {
                    Some(self.extension.as_slice())
                },
            )?;
            encoder.finish()
        };
        let receipt = encode().map_err(|error| SystemTypeError::IdentityEncoding {
            detail: error.to_string(),
        })?;

        Ok(AdmittedSystem {
            fields: self.fields,
            atoms: self.atoms,
            equations: self.equations,
            scalar_convention,
            extension: self.extension,
            identity: receipt.id(),
        })
    }

    /// Admit an expression: explicit-stack two-pass (validate + infer)
    /// so adversarial nesting cannot overflow the call stack.
    fn admit_expr(&self, root: &SystemExpr) -> Result<AdmittedExpr, SystemTypeError> {
        // Pass 1: depth + reference validation, iteratively.
        let mut stack: Vec<(&SystemExpr, usize)> = vec![(root, 1)];
        while let Some((node, depth)) = stack.pop() {
            if depth > MAX_SYSTEM_EXPR_DEPTH {
                return Err(SystemTypeError::DepthExceeded {
                    cap: MAX_SYSTEM_EXPR_DEPTH,
                });
            }
            match node {
                SystemExpr::FieldRef(field) => {
                    if field.0 >= self.fields.len() {
                        return Err(SystemTypeError::UnknownId {
                            what: "field",
                            id: field.0,
                        });
                    }
                }
                SystemExpr::Apply { atom, arg } => {
                    if *atom >= self.atoms.len() {
                        return Err(SystemTypeError::UnknownId {
                            what: "atom",
                            id: *atom,
                        });
                    }
                    stack.push((arg, depth + 1));
                }
                SystemExpr::Scale(value, inner) => {
                    if !value.is_finite() {
                        return Err(SystemTypeError::NonFiniteScale {
                            bits: value.to_bits(),
                        });
                    }
                    stack.push((inner, depth + 1));
                }
                SystemExpr::Add(left, right) => {
                    stack.push((left, depth + 1));
                    stack.push((right, depth + 1));
                }
                SystemExpr::PortPair { effort, flow, .. } => {
                    stack.push((effort, depth + 1));
                    stack.push((flow, depth + 1));
                }
            }
        }
        // Pass 2: post-order type inference with an explicit value stack.
        let mut work = vec![Step::Enter(root)];
        let mut values: Vec<AdmittedExpr> = Vec::new();
        while let Some(step) = work.pop() {
            match step {
                Step::Enter(node) => {
                    work.push(Step::Exit(node));
                    match node {
                        SystemExpr::FieldRef(_) => {}
                        SystemExpr::Apply { arg, .. } => work.push(Step::Enter(arg)),
                        SystemExpr::Scale(_, inner) => work.push(Step::Enter(inner)),
                        SystemExpr::Add(left, right) => {
                            work.push(Step::Enter(right));
                            work.push(Step::Enter(left));
                        }
                        SystemExpr::PortPair { effort, flow, .. } => {
                            work.push(Step::Enter(flow));
                            work.push(Step::Enter(effort));
                        }
                    }
                }
                Step::Exit(node) => {
                    let admitted = self.admit_node(node, &mut values)?;
                    values.push(admitted);
                }
            }
        }
        debug_assert_eq!(values.len(), 1);
        values.pop().ok_or(SystemTypeError::IdentityEncoding {
            detail: "expression evaluation stack imbalance".to_string(),
        })
    }

    #[allow(clippy::too_many_lines)] // One exhaustive match: every SystemExpr node's admission rule in one place.
    fn admit_node(
        &self,
        node: &SystemExpr,
        values: &mut Vec<AdmittedExpr>,
    ) -> Result<AdmittedExpr, SystemTypeError> {
        match node {
            SystemExpr::FieldRef(field) => {
                let declared = &self.fields[field.0];
                Ok(AdmittedExpr {
                    space: declared.space,
                    convention: Some(FieldConvention {
                        basis: declared.coordinates.basis.clone(),
                        frame: declared.coordinates.frame.clone(),
                        clock: declared.clock.clone(),
                    }),
                    affine_field: declared
                        .quantity
                        .is_affine_absolute()
                        .then(|| declared.name.clone()),
                })
            }
            SystemExpr::Apply { atom, .. } => {
                let arg = values.pop().expect("apply argument admitted");
                let signature = &self.atoms[*atom];
                if arg.space != signature.in_space {
                    return Err(SystemTypeError::SpaceMismatch {
                        context: "atom application",
                        expected: signature.in_space,
                        found: arg.space,
                    });
                }
                // Applying an operator to an affine absolute quantity is a
                // difference-producing act only if the atom says so; v1 takes
                // the conservative route and refuses.
                if let Some(field) = arg.affine_field {
                    return Err(SystemTypeError::AffineQuantityMisuse {
                        operation: "atom application",
                        field,
                    });
                }
                Ok(AdmittedExpr {
                    space: signature.out_space,
                    convention: arg.convention,
                    affine_field: None,
                })
            }
            SystemExpr::Scale(_, _) => {
                let inner = values.pop().expect("scale operand admitted");
                if let Some(field) = inner.affine_field {
                    return Err(SystemTypeError::AffineQuantityMisuse {
                        operation: "scaling",
                        field,
                    });
                }
                Ok(inner)
            }
            SystemExpr::Add(_, _) => {
                let right = values.pop().expect("add rhs admitted");
                let left = values.pop().expect("add lhs admitted");
                if let Some(field) = left.affine_field.or(right.affine_field) {
                    return Err(SystemTypeError::AffineQuantityMisuse {
                        operation: "summation",
                        field,
                    });
                }
                if left.space != right.space {
                    return Err(SystemTypeError::SpaceMismatch {
                        context: "summation",
                        expected: left.space,
                        found: right.space,
                    });
                }
                let convention = match (&left.convention, &right.convention) {
                    (Some(a), Some(b)) => {
                        if a.clock != b.clock {
                            return Err(SystemTypeError::ClockMismatch {
                                left: a.clock.as_str().to_string(),
                                right: b.clock.as_str().to_string(),
                            });
                        }
                        if a.basis != b.basis || a.frame != b.frame {
                            return Err(SystemTypeError::ConventionMismatch {
                                left: (a.basis.as_str().to_string(), a.frame.as_str().to_string()),
                                right: (b.basis.as_str().to_string(), b.frame.as_str().to_string()),
                            });
                        }
                        Some(a.clone())
                    }
                    (Some(a), None) | (None, Some(a)) => Some(a.clone()),
                    (None, None) => None,
                };
                Ok(AdmittedExpr {
                    space: left.space,
                    convention,
                    affine_field: None,
                })
            }
            SystemExpr::PortPair {
                kind, measure_dims, ..
            } => {
                let flow = values.pop().expect("pair flow admitted");
                let effort = values.pop().expect("pair effort admitted");
                if let (Some(a), Some(b)) = (&effort.convention, &flow.convention)
                    && a.clock != b.clock
                {
                    return Err(SystemTypeError::ClockMismatch {
                        left: a.clock.as_str().to_string(),
                        right: b.clock.as_str().to_string(),
                    });
                }
                let product = effort
                    .space
                    .dims
                    .checked_plus(flow.space.dims)
                    .and_then(|sum| sum.checked_plus(*measure_dims));
                if product != Some(POWER_DIMS) {
                    return Err(SystemTypeError::NonConjugatePairing {
                        kind: *kind,
                        effort_dims: effort.space.dims,
                        flow_dims: flow.space.dims,
                        measure_dims: *measure_dims,
                    });
                }
                // A pairing yields a scalar power density block.
                Ok(AdmittedExpr {
                    space: Space {
                        degree: 0,
                        n: 1,
                        dims: POWER_DIMS,
                    },
                    convention: effort.convention,
                    affine_field: None,
                })
            }
        }
    }
}

enum Step<'e> {
    Enter(&'e SystemExpr),
    Exit(&'e SystemExpr),
}

struct FieldConvention {
    basis: ConventionRef,
    frame: ConventionRef,
    clock: ConventionRef,
}

impl Clone for FieldConvention {
    fn clone(&self) -> Self {
        Self {
            basis: self.basis.clone(),
            frame: self.frame.clone(),
            clock: self.clock.clone(),
        }
    }
}

struct AdmittedExpr {
    space: Space,
    convention: Option<FieldConvention>,
    affine_field: Option<String>,
}

fn check_convention_match(
    convention: &FieldConvention,
    target: &FieldDecl,
) -> Result<(), SystemTypeError> {
    if convention.clock != target.clock {
        return Err(SystemTypeError::ClockMismatch {
            left: convention.clock.as_str().to_string(),
            right: target.clock.as_str().to_string(),
        });
    }
    if convention.basis != target.coordinates.basis || convention.frame != target.coordinates.frame
    {
        return Err(SystemTypeError::ConventionMismatch {
            left: (
                convention.basis.as_str().to_string(),
                convention.frame.as_str().to_string(),
            ),
            right: (
                target.coordinates.basis.as_str().to_string(),
                target.coordinates.frame.as_str().to_string(),
            ),
        });
    }
    Ok(())
}

// ---- canonical byte payloads (display names excluded) ----

fn push_u32(bytes: &mut Vec<u8>, value: u32) {
    bytes.extend_from_slice(&value.to_be_bytes());
}

fn push_dims(bytes: &mut Vec<u8>, dims: Dims) {
    for exponent in dims.0 {
        bytes.push(exponent.cast_unsigned());
    }
}

fn push_ref(bytes: &mut Vec<u8>, reference: &ConventionRef) {
    push_u32(
        bytes,
        u32::try_from(reference.as_str().len()).expect("bounded ref"),
    );
    bytes.extend_from_slice(reference.as_str().as_bytes());
}

fn push_space(bytes: &mut Vec<u8>, space: &Space) {
    bytes.push(space.degree);
    bytes.extend_from_slice(&(space.n as u64).to_be_bytes());
    push_dims(bytes, space.dims);
}

fn canonical_field_bytes(field: &FieldDecl) -> Vec<u8> {
    let mut bytes = Vec::new();
    push_space(&mut bytes, &field.space);
    match &field.quantity {
        FieldQuantity::Dimensional(dims) => {
            bytes.push(0);
            push_dims(&mut bytes, *dims);
        }
        FieldQuantity::Semantic(semantic) => {
            bytes.push(1);
            // Debug rendering of the sealed kind/form pair is stable and
            // versioned by SYSTEM_IR_VERSION; fs-qty exposes no canonical
            // byte encoding for SemanticType yet (documented no-claim).
            let rendered = format!("{semantic:?}");
            push_u32(&mut bytes, u32::try_from(rendered.len()).expect("bounded"));
            bytes.extend_from_slice(rendered.as_bytes());
        }
    }
    push_ref(&mut bytes, &field.coordinates.basis);
    push_ref(&mut bytes, &field.coordinates.frame);
    bytes.push(match field.coordinates.orientation {
        PortOrientation::OutwardFromOwner => 0,
        PortOrientation::AlongFrame => 1,
        PortOrientation::AgainstFrame => 2,
    });
    push_ref(&mut bytes, &field.clock);
    bytes.push(match field.support {
        SpatialSupport::Interior => 0,
        SpatialSupport::BoundaryTrace => 1,
    });
    match field.state {
        StateOwnership::Owned { slot } => {
            bytes.push(0);
            push_u32(&mut bytes, slot);
        }
        StateOwnership::External => bytes.push(1),
        StateOwnership::Parameter { role } => {
            bytes.push(2);
            bytes.push(match role {
                ParameterRole::Design => 0,
                ParameterRole::Material => 1,
                ParameterRole::Control => 2,
            });
        }
    }
    bytes
}

fn canonical_atom_bytes(atom: &AtomSignature) -> Vec<u8> {
    let mut bytes = Vec::new();
    push_ref(&mut bytes, &atom.content);
    push_space(&mut bytes, &atom.in_space);
    push_space(&mut bytes, &atom.out_space);
    bytes
}

fn canonical_equation_bytes(
    equation: &BlockEquation,
    remap: &[usize],
    atom_remap: &[usize],
) -> Vec<u8> {
    let mut bytes = Vec::new();
    push_u32(
        &mut bytes,
        u32::try_from(remap[equation.target.0]).expect("bounded field table"),
    );
    // Iterative pre-order serialization with explicit child counts: the
    // tree shape is unambiguous without recursion.
    let mut stack = vec![&equation.rhs];
    while let Some(node) = stack.pop() {
        match node {
            SystemExpr::FieldRef(field) => {
                bytes.push(0);
                push_u32(
                    &mut bytes,
                    u32::try_from(remap[field.0]).expect("bounded field table"),
                );
            }
            SystemExpr::Apply { atom, arg } => {
                bytes.push(1);
                push_u32(
                    &mut bytes,
                    u32::try_from(atom_remap[*atom]).expect("bounded atom table"),
                );
                stack.push(arg);
            }
            SystemExpr::Scale(value, inner) => {
                bytes.push(2);
                bytes.extend_from_slice(&value.to_bits().to_be_bytes());
                stack.push(inner);
            }
            SystemExpr::Add(left, right) => {
                bytes.push(3);
                stack.push(right);
                stack.push(left);
            }
            SystemExpr::PortPair {
                kind,
                effort,
                flow,
                measure_dims,
            } => {
                bytes.push(4);
                let rendered = format!("{kind:?}");
                push_u32(&mut bytes, u32::try_from(rendered.len()).expect("bounded"));
                bytes.extend_from_slice(rendered.as_bytes());
                push_dims(&mut bytes, *measure_dims);
                stack.push(flow);
                stack.push(effort);
            }
        }
    }
    bytes
}

/// Versioned canonical text transport for admitted systems: one
/// LF-terminated line per record, tab-separated fields, strict
/// fail-closed parsing. The transport binds [`SYSTEM_IR_VERSION`];
/// any other version refuses with [`SystemTypeError::VersionMismatch`]
/// so migrations stay explicit and audited. Import NEVER trusts the
/// payload: it rebuilds a [`SystemDef`] and re-runs full admission, so
/// a tampered transport either refuses structurally or mints a
/// different [`SystemId`] — it cannot alias the original identity.
pub mod transport {
    use super::{
        AdmittedSystem, AtomSignature, BlockEquation, ConventionRef, CoordinateConvention,
        FieldDecl, FieldId, FieldQuantity, ParameterRole, SYSTEM_IR_VERSION, ScalarConvention,
        SpatialSupport, StateOwnership, SystemDef, SystemExpr, SystemTypeError,
    };
    use crate::expr::Space;
    use fs_couple::{PortKind, PortOrientation};
    use fs_qty::Dims;
    use fs_qty::semantic::{
        AngleDomain, CompositionBasis, QuantityKind, SemanticType, StrainBasis, StrainComponent,
        ValueForm,
    };
    use std::fmt::Write as _;

    const MAGIC: &str = "fs-opdsl-system-transport-v1";

    enum Frame {
        Apply {
            atom: usize,
        },
        Scale {
            bits: u64,
        },
        Add {
            left: Option<SystemExpr>,
        },
        Pair {
            kind: PortKind,
            measure: Dims,
            effort: Option<SystemExpr>,
        },
    }

    fn validate_name(name: &str) -> Result<(), SystemTypeError> {
        if name.is_empty() || name.len() > 256 || name.chars().any(|c| c.is_control() || c == '\t')
        {
            return Err(SystemTypeError::InvalidName {
                name: name.to_string(),
            });
        }
        Ok(())
    }

    fn render_space(space: &Space) -> String {
        let d = space.dims.0;
        format!(
            "{},{},{},{},{},{},{},{}",
            space.degree, space.n, d[0], d[1], d[2], d[3], d[4], d[5]
        )
    }

    fn render_dims(dims: Dims) -> String {
        let d = dims.0;
        format!("{},{},{},{},{},{}", d[0], d[1], d[2], d[3], d[4], d[5])
    }

    fn render_kind(kind: QuantityKind) -> String {
        match kind {
            QuantityKind::AbsoluteTemperature => "absolute-temperature".to_string(),
            QuantityKind::TemperatureDifference => "temperature-difference".to_string(),
            QuantityKind::Angle(domain) => format!("angle:{}", render_angle_domain(domain)),
            QuantityKind::AngularVelocity(domain) => {
                format!("angular-velocity:{}", render_angle_domain(domain))
            }
            QuantityKind::Torque => "torque".to_string(),
            QuantityKind::Energy => "energy".to_string(),
            QuantityKind::Pressure => "pressure".to_string(),
            QuantityKind::Stress => "stress".to_string(),
            QuantityKind::Strain { basis, component } => format!(
                "strain:{}:{}",
                match basis {
                    StrainBasis::Tensor => "tensor",
                    StrainBasis::Engineering => "engineering",
                },
                match component {
                    StrainComponent::Normal => "normal",
                    StrainComponent::Shear => "shear",
                }
            ),
            QuantityKind::Composition(basis) => format!(
                "composition:{}",
                match basis {
                    CompositionBasis::MassFraction => "mass",
                    CompositionBasis::MoleFraction => "mole",
                    CompositionBasis::VolumeFraction => "volume",
                }
            ),
            QuantityKind::Mass => "mass".to_string(),
            QuantityKind::Amount => "amount".to_string(),
            QuantityKind::MolarMass => "molar-mass".to_string(),
            QuantityKind::MassConcentration => "mass-concentration".to_string(),
            QuantityKind::AmountConcentration => "amount-concentration".to_string(),
            QuantityKind::Entropy => "entropy".to_string(),
            QuantityKind::HeatCapacity => "heat-capacity".to_string(),
            QuantityKind::AcousticPressure => "acoustic-pressure".to_string(),
            QuantityKind::AcousticPower => "acoustic-power".to_string(),
        }
    }

    fn render_angle_domain(domain: AngleDomain) -> &'static str {
        match domain {
            AngleDomain::Mechanical => "mechanical",
            AngleDomain::Electrical => "electrical",
        }
    }

    fn render_form(form: ValueForm) -> &'static str {
        match form {
            ValueForm::Static => "static",
            ValueForm::Instantaneous => "instantaneous",
            ValueForm::Peak => "peak",
            ValueForm::Rms => "rms",
        }
    }

    fn render_port_kind(kind: PortKind) -> &'static str {
        match kind {
            PortKind::MechanicalForceVelocity => "mechanical-force-velocity",
            PortKind::FluidPressureFlux => "fluid-pressure-flux",
            PortKind::ThermalTemperatureEntropy => "thermal-temperature-entropy",
            PortKind::RotationalTorqueAngularVelocity => "rotational-torque-angular-velocity",
            PortKind::ElectricalVoltageCurrent => "electrical-voltage-current",
            PortKind::MagneticMmfFluxRate => "magnetic-mmf-flux-rate",
            PortKind::ChemicalPotentialAmountFlow => "chemical-potential-amount-flow",
        }
    }

    fn render_expr(root: &SystemExpr, out: &mut String) {
        // Pre-order token stream; arity is implied by each token.
        let mut stack = vec![root];
        let mut first = true;
        while let Some(node) = stack.pop() {
            if !first {
                out.push(' ');
            }
            first = false;
            match node {
                SystemExpr::FieldRef(field) => {
                    let _ = write!(out, "f{}", field.0);
                }
                SystemExpr::Apply { atom, arg } => {
                    let _ = write!(out, "a{atom}");
                    stack.push(arg);
                }
                SystemExpr::Scale(value, inner) => {
                    let _ = write!(out, "x{:016x}", value.to_bits());
                    stack.push(inner);
                }
                SystemExpr::Add(left, right) => {
                    out.push('+');
                    stack.push(right);
                    stack.push(left);
                }
                SystemExpr::PortPair {
                    kind,
                    effort,
                    flow,
                    measure_dims,
                } => {
                    let _ = write!(
                        out,
                        "p{}:{}",
                        render_port_kind(*kind),
                        render_dims(*measure_dims)
                    );
                    stack.push(flow);
                    stack.push(effort);
                }
            }
        }
    }

    /// Serialize an admitted system as the canonical v1 transport text.
    ///
    /// # Errors
    /// [`SystemTypeError::InvalidName`] when a display name cannot be
    /// carried fail-closed (control characters or tabs).
    pub fn to_text(system: &AdmittedSystem) -> Result<String, SystemTypeError> {
        let mut out = String::new();
        out.push_str(MAGIC);
        out.push('\n');
        let _ = writeln!(out, "version\t{SYSTEM_IR_VERSION}");
        let _ = writeln!(
            out,
            "convention\t{}",
            match system.scalar_convention() {
                ScalarConvention::RealOnly => "real-only",
                ScalarConvention::ComplexHermitian => "complex-hermitian",
            }
        );
        if !system.extension().is_empty() {
            out.push_str("extension\t");
            for byte in system.extension() {
                let _ = write!(out, "{byte:02x}");
            }
            out.push('\n');
        }
        for atom in system.atoms() {
            validate_name(&atom.name)?;
            let _ = writeln!(
                out,
                "atom\t{}\t{}\t{}\t{}",
                atom.name,
                atom.content.as_str(),
                render_space(&atom.in_space),
                render_space(&atom.out_space)
            );
        }
        for field in system.fields() {
            validate_name(&field.name)?;
            let quantity = match &field.quantity {
                FieldQuantity::Dimensional(dims) => format!("dims:{}", render_dims(*dims)),
                FieldQuantity::Semantic(semantic) => format!(
                    "kind:{}:{}",
                    render_kind(semantic.kind()),
                    render_form(semantic.form())
                ),
            };
            let state = match field.state {
                StateOwnership::Owned { slot } => format!("owned:{slot}"),
                StateOwnership::External => "external".to_string(),
                StateOwnership::Parameter { role } => format!(
                    "parameter:{}",
                    match role {
                        ParameterRole::Design => "design",
                        ParameterRole::Material => "material",
                        ParameterRole::Control => "control",
                    }
                ),
            };
            let _ = writeln!(
                out,
                "field\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                field.name,
                render_space(&field.space),
                quantity,
                field.coordinates.basis.as_str(),
                field.coordinates.frame.as_str(),
                match field.coordinates.orientation {
                    PortOrientation::OutwardFromOwner => "outward",
                    PortOrientation::AlongFrame => "along",
                    PortOrientation::AgainstFrame => "against",
                },
                field.clock.as_str(),
                match field.support {
                    SpatialSupport::Interior => "interior",
                    SpatialSupport::BoundaryTrace => "boundary-trace",
                },
                state
            );
        }
        for equation in system.equations() {
            validate_name(&equation.name)?;
            let mut rendered = String::new();
            render_expr(&equation.rhs, &mut rendered);
            let _ = writeln!(
                out,
                "equation\t{}\t{}\t{}",
                equation.name, equation.target.0, rendered
            );
        }
        Ok(out)
    }

    fn malformed(line: usize, detail: impl Into<String>) -> SystemTypeError {
        SystemTypeError::TransportMalformed {
            line,
            detail: detail.into(),
        }
    }

    fn parse_i8(token: &str, line: usize) -> Result<i8, SystemTypeError> {
        token
            .parse::<i8>()
            .map_err(|_| malformed(line, format!("{token:?} is not an i8 exponent")))
    }

    fn parse_dims(token: &str, line: usize) -> Result<Dims, SystemTypeError> {
        let parts: Vec<&str> = token.split(',').collect();
        if parts.len() != 6 {
            return Err(malformed(
                line,
                format!("dims {token:?} must have 6 exponents"),
            ));
        }
        let mut dims = [0i8; 6];
        for (slot, part) in dims.iter_mut().zip(parts) {
            *slot = parse_i8(part, line)?;
        }
        Ok(Dims(dims))
    }

    fn parse_space(token: &str, line: usize) -> Result<Space, SystemTypeError> {
        let parts: Vec<&str> = token.split(',').collect();
        if parts.len() != 8 {
            return Err(malformed(
                line,
                format!("space {token:?} must have 8 members"),
            ));
        }
        let degree = parts[0]
            .parse::<u8>()
            .map_err(|_| malformed(line, "space degree must be u8"))?;
        let n = parts[1]
            .parse::<usize>()
            .map_err(|_| malformed(line, "space dof count must be usize"))?;
        let mut dims = [0i8; 6];
        for (slot, part) in dims.iter_mut().zip(&parts[2..]) {
            *slot = parse_i8(part, line)?;
        }
        Ok(Space {
            degree,
            n,
            dims: Dims(dims),
        })
    }

    fn parse_angle_domain(token: &str, line: usize) -> Result<AngleDomain, SystemTypeError> {
        match token {
            "mechanical" => Ok(AngleDomain::Mechanical),
            "electrical" => Ok(AngleDomain::Electrical),
            other => Err(malformed(line, format!("unknown angle domain {other:?}"))),
        }
    }

    fn parse_quantity(token: &str, line: usize) -> Result<FieldQuantity, SystemTypeError> {
        if let Some(rest) = token.strip_prefix("dims:") {
            return Ok(FieldQuantity::Dimensional(parse_dims(rest, line)?));
        }
        let Some(rest) = token.strip_prefix("kind:") else {
            return Err(malformed(
                line,
                format!("quantity {token:?} must be dims: or kind:"),
            ));
        };
        let Some((kind_token, form_token)) = rest.rsplit_once(':') else {
            return Err(malformed(
                line,
                format!("kind {rest:?} must end with :<form>"),
            ));
        };
        let kind = match kind_token {
            "absolute-temperature" => QuantityKind::AbsoluteTemperature,
            "temperature-difference" => QuantityKind::TemperatureDifference,
            "torque" => QuantityKind::Torque,
            "energy" => QuantityKind::Energy,
            "pressure" => QuantityKind::Pressure,
            "stress" => QuantityKind::Stress,
            "mass" => QuantityKind::Mass,
            "amount" => QuantityKind::Amount,
            "molar-mass" => QuantityKind::MolarMass,
            "mass-concentration" => QuantityKind::MassConcentration,
            "amount-concentration" => QuantityKind::AmountConcentration,
            "entropy" => QuantityKind::Entropy,
            "heat-capacity" => QuantityKind::HeatCapacity,
            "acoustic-pressure" => QuantityKind::AcousticPressure,
            "acoustic-power" => QuantityKind::AcousticPower,
            other => {
                if let Some(domain) = other.strip_prefix("angle:") {
                    QuantityKind::Angle(parse_angle_domain(domain, line)?)
                } else if let Some(domain) = other.strip_prefix("angular-velocity:") {
                    QuantityKind::AngularVelocity(parse_angle_domain(domain, line)?)
                } else if let Some(rest) = other.strip_prefix("strain:") {
                    let Some((basis, component)) = rest.split_once(':') else {
                        return Err(malformed(line, "strain kind needs basis and component"));
                    };
                    QuantityKind::Strain {
                        basis: match basis {
                            "tensor" => StrainBasis::Tensor,
                            "engineering" => StrainBasis::Engineering,
                            unknown => {
                                return Err(malformed(
                                    line,
                                    format!("unknown strain basis {unknown:?}"),
                                ));
                            }
                        },
                        component: match component {
                            "normal" => StrainComponent::Normal,
                            "shear" => StrainComponent::Shear,
                            unknown => {
                                return Err(malformed(
                                    line,
                                    format!("unknown strain component {unknown:?}"),
                                ));
                            }
                        },
                    }
                } else if let Some(basis) = other.strip_prefix("composition:") {
                    QuantityKind::Composition(match basis {
                        "mass" => CompositionBasis::MassFraction,
                        "mole" => CompositionBasis::MoleFraction,
                        "volume" => CompositionBasis::VolumeFraction,
                        unknown => {
                            return Err(malformed(
                                line,
                                format!("unknown composition basis {unknown:?}"),
                            ));
                        }
                    })
                } else {
                    return Err(malformed(line, format!("unknown quantity kind {other:?}")));
                }
            }
        };
        let form = match form_token {
            "static" => ValueForm::Static,
            "instantaneous" => ValueForm::Instantaneous,
            "peak" => ValueForm::Peak,
            "rms" => ValueForm::Rms,
            unknown => return Err(malformed(line, format!("unknown value form {unknown:?}"))),
        };
        Ok(FieldQuantity::Semantic(SemanticType::new(kind, form)))
    }

    fn parse_port_kind(token: &str, line: usize) -> Result<PortKind, SystemTypeError> {
        Ok(match token {
            "mechanical-force-velocity" => PortKind::MechanicalForceVelocity,
            "fluid-pressure-flux" => PortKind::FluidPressureFlux,
            "thermal-temperature-entropy" => PortKind::ThermalTemperatureEntropy,
            "rotational-torque-angular-velocity" => PortKind::RotationalTorqueAngularVelocity,
            "electrical-voltage-current" => PortKind::ElectricalVoltageCurrent,
            "magnetic-mmf-flux-rate" => PortKind::MagneticMmfFluxRate,
            "chemical-potential-amount-flow" => PortKind::ChemicalPotentialAmountFlow,
            other => return Err(malformed(line, format!("unknown port kind {other:?}"))),
        })
    }

    #[allow(clippy::too_many_lines)] // One arity-driven prefix parse: token dispatch and frame reduction are one loop.
    fn parse_expr(tokens: &[&str], line: usize) -> Result<SystemExpr, SystemTypeError> {
        // Arity-driven prefix parse with an explicit frame stack.
        let mut frames: Vec<Frame> = Vec::new();
        let mut tokens = tokens.iter();
        let mut pending: Option<SystemExpr> = None;
        loop {
            // Feed completed values upward until a frame still wants more.
            while let Some(value) = pending.take() {
                match frames.pop() {
                    None => {
                        if tokens.next().is_some() {
                            return Err(malformed(line, "trailing expression tokens"));
                        }
                        return Ok(value);
                    }
                    Some(Frame::Apply { atom }) => {
                        pending = Some(SystemExpr::Apply {
                            atom,
                            arg: Box::new(value),
                        });
                    }
                    Some(Frame::Scale { bits }) => {
                        pending = Some(SystemExpr::Scale(f64::from_bits(bits), Box::new(value)));
                    }
                    Some(Frame::Add { left: None }) => {
                        frames.push(Frame::Add { left: Some(value) });
                    }
                    Some(Frame::Add { left: Some(left) }) => {
                        pending = Some(SystemExpr::Add(Box::new(left), Box::new(value)));
                    }
                    Some(Frame::Pair {
                        kind,
                        measure,
                        effort: None,
                    }) => {
                        frames.push(Frame::Pair {
                            kind,
                            measure,
                            effort: Some(value),
                        });
                    }
                    Some(Frame::Pair {
                        kind,
                        measure,
                        effort: Some(effort),
                    }) => {
                        pending = Some(SystemExpr::PortPair {
                            kind,
                            effort: Box::new(effort),
                            flow: Box::new(value),
                            measure_dims: measure,
                        });
                    }
                }
            }
            let Some(token) = tokens.next() else {
                return Err(malformed(line, "expression ended before all operands"));
            };
            if let Some(index) = token.strip_prefix('f') {
                let index = index
                    .parse::<usize>()
                    .map_err(|_| malformed(line, format!("field token {token:?}")))?;
                pending = Some(SystemExpr::FieldRef(FieldId(index)));
            } else if let Some(index) = token.strip_prefix('a') {
                let atom = index
                    .parse::<usize>()
                    .map_err(|_| malformed(line, format!("atom token {token:?}")))?;
                frames.push(Frame::Apply { atom });
            } else if let Some(hex) = token.strip_prefix('x') {
                let bits = u64::from_str_radix(hex, 16)
                    .map_err(|_| malformed(line, format!("scale token {token:?}")))?;
                frames.push(Frame::Scale { bits });
            } else if *token == "+" {
                frames.push(Frame::Add { left: None });
            } else if let Some(rest) = token.strip_prefix('p') {
                let Some((kind_token, dims_token)) = rest.rsplit_once(':') else {
                    return Err(malformed(line, format!("pair token {token:?}")));
                };
                frames.push(Frame::Pair {
                    kind: parse_port_kind(kind_token, line)?,
                    measure: parse_dims(dims_token, line)?,
                    effort: None,
                });
            } else {
                return Err(malformed(
                    line,
                    format!("unknown expression token {token:?}"),
                ));
            }
            if frames.len() > super::MAX_SYSTEM_EXPR_DEPTH {
                return Err(SystemTypeError::DepthExceeded {
                    cap: super::MAX_SYSTEM_EXPR_DEPTH,
                });
            }
        }
    }

    /// Parse canonical v1 transport text back into a [`SystemDef`].
    /// Everything re-admits from scratch — a tampered payload either
    /// refuses here, refuses at admission, or mints a different
    /// identity.
    ///
    /// # Errors
    /// [`SystemTypeError::VersionMismatch`] for any other IR version;
    /// [`SystemTypeError::TransportMalformed`] for structural refusals;
    /// every admission error the rebuilt definition produces.
    #[allow(clippy::too_many_lines)] // One strict line-dispatch transaction: every record grammar in one place.
    pub fn from_text(text: &str) -> Result<SystemDef, SystemTypeError> {
        let mut lines = text.lines().enumerate();
        let Some((_, magic)) = lines.next() else {
            return Err(malformed(1, "empty transport"));
        };
        if magic != MAGIC {
            return Err(malformed(1, format!("magic {magic:?} is not {MAGIC:?}")));
        }
        let Some((_, version_line)) = lines.next() else {
            return Err(malformed(2, "missing version line"));
        };
        let version = version_line
            .strip_prefix("version\t")
            .and_then(|value| value.parse::<u32>().ok())
            .ok_or_else(|| malformed(2, "version line must be 'version\\t<u32>'"))?;
        if version != SYSTEM_IR_VERSION {
            return Err(SystemTypeError::VersionMismatch {
                found: version,
                supported: SYSTEM_IR_VERSION,
            });
        }
        let Some((_, convention_line)) = lines.next() else {
            return Err(malformed(3, "missing convention line"));
        };
        let convention = match convention_line.strip_prefix("convention\t") {
            Some("real-only") => ScalarConvention::RealOnly,
            Some("complex-hermitian") => ScalarConvention::ComplexHermitian,
            _ => {
                return Err(malformed(
                    3,
                    "convention must be real-only or complex-hermitian",
                ));
            }
        };
        let mut system = SystemDef::new().scalar_convention(convention);
        let mut equations: Vec<(usize, BlockEquation)> = Vec::new();
        for (index, raw) in lines {
            let line = index + 1;
            let mut parts = raw.split('\t');
            let record = parts.next().unwrap_or_default();
            let rest: Vec<&str> = parts.collect();
            match record {
                "extension" => {
                    let [hex] = rest.as_slice() else {
                        return Err(malformed(line, "extension takes exactly one hex payload"));
                    };
                    if hex.is_empty() || hex.len() % 2 != 0 {
                        return Err(malformed(line, "extension hex must be non-empty and even"));
                    }
                    let mut bytes = Vec::with_capacity(hex.len() / 2);
                    for pair in hex.as_bytes().chunks(2) {
                        let text = std::str::from_utf8(pair)
                            .map_err(|_| malformed(line, "extension hex must be ASCII"))?;
                        bytes.push(
                            u8::from_str_radix(text, 16)
                                .map_err(|_| malformed(line, "extension must be lowercase hex"))?,
                        );
                    }
                    system = system.with_extension(bytes)?;
                }
                "atom" => {
                    let [name, content, in_space, out_space] = rest.as_slice() else {
                        return Err(malformed(
                            line,
                            "atom takes name, content-ref, in-space, out-space",
                        ));
                    };
                    validate_name(name)?;
                    system.register_atom(AtomSignature {
                        name: (*name).to_string(),
                        content: ConventionRef::new(*content)?,
                        in_space: parse_space(in_space, line)?,
                        out_space: parse_space(out_space, line)?,
                    });
                }
                "field" => {
                    let [
                        name,
                        space,
                        quantity,
                        basis,
                        frame,
                        orientation,
                        clock,
                        support,
                        state,
                    ] = rest.as_slice()
                    else {
                        return Err(malformed(line, "field takes 9 members"));
                    };
                    validate_name(name)?;
                    let state = if let Some(slot) = state.strip_prefix("owned:") {
                        StateOwnership::Owned {
                            slot: slot
                                .parse::<u32>()
                                .map_err(|_| malformed(line, "owned slot must be u32"))?,
                        }
                    } else if *state == "external" {
                        StateOwnership::External
                    } else if let Some(role) = state.strip_prefix("parameter:") {
                        StateOwnership::Parameter {
                            role: match role {
                                "design" => ParameterRole::Design,
                                "material" => ParameterRole::Material,
                                "control" => ParameterRole::Control,
                                unknown => {
                                    return Err(malformed(
                                        line,
                                        format!("unknown parameter role {unknown:?}"),
                                    ));
                                }
                            },
                        }
                    } else {
                        return Err(malformed(line, format!("unknown state {state:?}")));
                    };
                    system.declare_field(FieldDecl {
                        name: (*name).to_string(),
                        space: parse_space(space, line)?,
                        quantity: parse_quantity(quantity, line)?,
                        coordinates: CoordinateConvention {
                            basis: ConventionRef::new(*basis)?,
                            frame: ConventionRef::new(*frame)?,
                            orientation: match *orientation {
                                "outward" => PortOrientation::OutwardFromOwner,
                                "along" => PortOrientation::AlongFrame,
                                "against" => PortOrientation::AgainstFrame,
                                unknown => {
                                    return Err(malformed(
                                        line,
                                        format!("unknown orientation {unknown:?}"),
                                    ));
                                }
                            },
                        },
                        clock: ConventionRef::new(*clock)?,
                        support: match *support {
                            "interior" => SpatialSupport::Interior,
                            "boundary-trace" => SpatialSupport::BoundaryTrace,
                            unknown => {
                                return Err(malformed(
                                    line,
                                    format!("unknown support {unknown:?}"),
                                ));
                            }
                        },
                        state,
                    })?;
                }
                "equation" => {
                    let [name, target, tokens @ ..] = rest.as_slice() else {
                        return Err(malformed(line, "equation takes name, target, expression"));
                    };
                    validate_name(name)?;
                    let target = target
                        .parse::<usize>()
                        .map_err(|_| malformed(line, "equation target must be an index"))?;
                    let [expression] = tokens else {
                        return Err(malformed(line, "equation takes exactly one expression"));
                    };
                    let tokens: Vec<&str> = expression.split(' ').collect();
                    equations.push((
                        line,
                        BlockEquation {
                            name: (*name).to_string(),
                            target: FieldId(target),
                            rhs: parse_expr(&tokens, line)?,
                        },
                    ));
                }
                unknown => {
                    return Err(malformed(line, format!("unknown record {unknown:?}")));
                }
            }
        }
        for (_, equation) in equations {
            system.add_equation(equation)?;
        }
        Ok(system)
    }
}
