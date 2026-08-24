//! The 1D elliptic testbed: P1 finite elements for `−u″ = f` on (0,1)
//! with homogeneous Dirichlet data, polynomial manufactured solutions
//! (so the mathematical Gauss rule is exact for the admitted integrands),
//! fallible bounded solves/oracles, and a Newton solver for the nonlinear
//! warm-start class. Stored antiderivative coefficients are rounded replay
//! metadata; the verifier intervalizes their exact construction separately.

use crate::interval::Iv;
use core::{fmt, mem::size_of};
use fs_obs::ident::{BoundedIdentityBuilder, IdentityBuildError, ReplayIdentity};

/// Largest mesh admitted by one synchronous v0 fem1d operation.
pub const MAX_FEM1D_MESH_NODES: usize = 1_000_000;
/// Polynomial exactness envelope: degree at most five.
pub const MAX_FEM1D_POLY_COEFFICIENTS: usize = 6;
/// Largest non-canonical coefficient vector inspected before zero trimming.
///
/// This bounds hostile-input work while allowing semantically redundant
/// trailing zeroes to normalize before the degree-five class cap is applied.
pub const MAX_FEM1D_RAW_POLY_COEFFICIENTS: usize = 4_096;
/// Largest Newton update budget accepted by the synchronous toy solver.
pub const MAX_FEM1D_NEWTON_ITERATIONS: u32 = 10_000;
/// Conservative scalar-work ceiling for one synchronous fem1d call.
pub const MAX_FEM1D_WORK_UNITS: usize = 50_000_000;
/// Largest manufactured-class name admitted at this boundary.
pub const MAX_FEM1D_CLASS_NAME_BYTES: usize = 4_096;
const IDENTITY_STREAM_FRAME_BYTES: usize = 4 + size_of::<u32>() + size_of::<u64>();
const IDENTITY_FIELD_FRAME_BYTES: usize = 1 + 2 * size_of::<u64>();

const fn identity_header_bytes(kind: &str) -> usize {
    IDENTITY_STREAM_FRAME_BYTES + kind.len()
}

const fn identity_field_bytes(key: &str, value_bytes: usize) -> usize {
    IDENTITY_FIELD_FRAME_BYTES + key.len() + value_bytes
}

/// Semantic schema for canonical manufactured-solution class identities.
///
/// Version 2 deliberately rotates the formerly unversioned artifact kind.
pub const MMS_CLASS_IDENTITY_VERSION: u32 = 2;
/// Versioned artifact domain for an admitted manufactured-solution class.
pub const MMS_CLASS_IDENTITY_DOMAIN: &str = "fs-verify/fem1d-mms-class.v2";
/// Semantic schema for canonical meshed manufactured-problem identities.
///
/// Version 2 deliberately rotates the formerly unversioned artifact kind and
/// embeds the v2 class identity.
pub const MMS_PROBLEM_IDENTITY_VERSION: u32 = 2;
/// Versioned artifact domain for an admitted meshed manufactured problem.
pub const MMS_PROBLEM_IDENTITY_DOMAIN: &str = "fs-verify/fem1d-mms-problem.v2";

/// Largest producer-v2 canonical replay identity for an admitted MMS class.
///
/// The polynomial payload maximum is 15 stored f64 coefficients: six for `u`,
/// four for `-u''`, and five for its antiderivative.
pub const MAX_FEM1D_CLASS_CANONICAL_IDENTITY_BYTES: usize =
    identity_header_bytes(MMS_CLASS_IDENTITY_DOMAIN)
        + identity_field_bytes("class_schema", size_of::<u64>())
        + identity_field_bytes("name", MAX_FEM1D_CLASS_NAME_BYTES)
        + identity_field_bytes(
            "exact_solution_f64_le",
            size_of::<u64>() * MAX_FEM1D_POLY_COEFFICIENTS,
        )
        + identity_field_bytes(
            "forcing_f64_le",
            size_of::<u64>() * (MAX_FEM1D_POLY_COEFFICIENTS - 2),
        )
        + identity_field_bytes(
            "rounded_forcing_antiderivative_f64_le",
            size_of::<u64>() * (MAX_FEM1D_POLY_COEFFICIENTS - 1),
        );
/// Largest producer-v2 canonical replay identity for an admitted meshed problem.
///
/// It binds the maximum class identity and eight exact bytes for each of the
/// 1,000,000 admitted mesh nodes.
pub const MAX_FEM1D_PROBLEM_CANONICAL_IDENTITY_BYTES: usize =
    identity_header_bytes(MMS_PROBLEM_IDENTITY_DOMAIN)
        + identity_field_bytes("problem_schema", size_of::<u64>())
        + identity_field_bytes("class", size_of::<u32>() + size_of::<u64>())
        + identity_field_bytes(
            "class_canonical_bytes",
            MAX_FEM1D_CLASS_CANONICAL_IDENTITY_BYTES,
        )
        + identity_field_bytes("mesh_f64_le", size_of::<u64>() * MAX_FEM1D_MESH_NODES);
/// Owner-local class-identity declaration consumed by `xtask check-identities`.
pub const MMS_CLASS_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-verify:fem1d-mms-class",
    "version_const=MMS_CLASS_IDENTITY_VERSION",
    "version=2",
    "domain=fs-verify/fem1d-mms-class.v2",
    "domain_const=MMS_CLASS_IDENTITY_DOMAIN",
    "encoder=class_identity",
    "encoder_helpers=class_identity_with_budget,f64_bytes,class_identity_from_parts_with_budget",
    "schema_constants=MMS_CLASS_IDENTITY_VERSION,MMS_CLASS_IDENTITY_DOMAIN,MAX_FEM1D_CLASS_CANONICAL_IDENTITY_BYTES,MAX_FEM1D_CLASS_NAME_BYTES,MAX_FEM1D_POLY_COEFFICIENTS,MAX_FEM1D_RAW_POLY_COEFFICIENTS,IDENTITY_STREAM_FRAME_BYTES,IDENTITY_FIELD_FRAME_BYTES",
    "schema_functions=identity_header_bytes,identity_field_bytes,identity_build_error,MmsClass::new,Poly::from_coefficients,Poly::is_exactly_zero_at_one,Poly::derive,Poly::neg,Poly::antiderive,canonicalize_zero,validate_identity,validate_retained_identity_parts,copy_identity_bytes,check_mms_identity_version,check_mms_class_identity_version,MmsClass::identity_receipt,MmsClass::admit_identity_receipt,crates/fs-obs/src/ident.rs#ReplayIdentity::root,crates/fs-obs/src/ident.rs#ReplayIdentity::canonical_bytes",
    "schema_types=crates/fs-verify/src/fem1d.rs#MmsIdentityAdmissionError",
    "schema_dependencies=fs-obs:replay-identity-frame",
    "digest=fnv1a64",
    "encoding=typed-binary",
    "sources=MmsClass,MmsClassIdentityReceipt",
    "source_fields=MmsClass.name:semantic,MmsClass.u:semantic,MmsClass.forcing:semantic,MmsClass.rounded_forcing_antiderivative:semantic,MmsClass.identity:derived:recomputed-from-class-fields,MmsClassIdentityReceipt.declared_identity_version:semantic,MmsClassIdentityReceipt.canonical_bytes:semantic,MmsClassIdentityReceipt.root:derived:validated-fnv-root-of-retained-canonical-bytes",
    "source_bindings=MmsClass.name>name-utf8,MmsClass.u>exact-solution-f64-le,MmsClass.forcing>forcing-f64-le,MmsClass.rounded_forcing_antiderivative>rounded-forcing-antiderivative-f64-le,MmsClassIdentityReceipt.declared_identity_version>retained-producer-version,MmsClassIdentityReceipt.canonical_bytes>retained-canonical-bytes",
    "external_semantic_fields=artifact-domain,class-schema-version",
    "semantic_fields=artifact-domain,class-schema-version,name-utf8,exact-solution-f64-le,forcing-f64-le,rounded-forcing-antiderivative-f64-le,retained-producer-version,retained-canonical-bytes",
    "excluded_fields=none",
    "consumers=MmsClass::new,MmsClassIdentityReceipt,MmsClass::admit_identity_receipt,MmsProblem::from_class,MmsProblem::with_mesh",
    "mutations=artifact-domain:crates/fs-verify/src/fem1d.rs#mms_class_identity_fields_move_independently,class-schema-version:crates/fs-verify/src/fem1d.rs#mms_class_identity_fields_move_independently,name-utf8:crates/fs-verify/src/fem1d.rs#mms_class_identity_fields_move_independently,exact-solution-f64-le:crates/fs-verify/src/fem1d.rs#mms_class_identity_fields_move_independently,forcing-f64-le:crates/fs-verify/src/fem1d.rs#mms_class_identity_fields_move_independently,rounded-forcing-antiderivative-f64-le:crates/fs-verify/src/fem1d.rs#mms_class_identity_fields_move_independently,retained-producer-version:crates/fs-verify/src/fem1d.rs#mms_class_identity_receipts_fail_closed,retained-canonical-bytes:crates/fs-verify/src/fem1d.rs#mms_class_identity_receipts_fail_closed",
    "nonsemantic_mutations=none",
    "field_guard=classify_mms_class_identity_fields",
    "transport_guard=MmsClassIdentityReceipt::from_retained_parts",
    "version_guard=crates/fs-verify/src/fem1d.rs#mms_class_identity_receipts_fail_closed",
    "coupling_surface=fs-verify:fem1d-mms-class",
];

/// Owner-local problem-identity declaration consumed by `xtask check-identities`.
pub const MMS_PROBLEM_IDENTITY_SCHEMA_DECLARATION: &[&str] = &[
    "frankensim-identity-schema-v1",
    "id=fs-verify:fem1d-mms-problem",
    "version_const=MMS_PROBLEM_IDENTITY_VERSION",
    "version=2",
    "domain=fs-verify/fem1d-mms-problem.v2",
    "domain_const=MMS_PROBLEM_IDENTITY_DOMAIN",
    "encoder=problem_identity",
    "encoder_helpers=problem_identity_with_budget,f64_bytes,problem_identity_from_parts_with_budget",
    "schema_constants=MMS_PROBLEM_IDENTITY_VERSION,MMS_PROBLEM_IDENTITY_DOMAIN,MAX_FEM1D_PROBLEM_CANONICAL_IDENTITY_BYTES,MAX_FEM1D_MESH_NODES,IDENTITY_STREAM_FRAME_BYTES,IDENTITY_FIELD_FRAME_BYTES",
    "schema_functions=identity_header_bytes,identity_field_bytes,identity_build_error,admit_mesh,validate_mesh,validate_mesh_size,canonicalize_zero,Poly::try_copy,validate_retained_identity_parts,copy_identity_bytes,check_mms_identity_version,check_mms_problem_identity_version,MmsClass::identity,MmsClass::canonical_bytes,MmsClass::name,MmsClass::exact_solution,MmsProblem::new,MmsProblem::from_class,MmsProblem::from_admitted_parts,MmsProblem::with_mesh,MmsProblem::identity_receipt,MmsProblem::admit_identity_receipt,crates/fs-obs/src/ident.rs#ReplayIdentity::root,crates/fs-obs/src/ident.rs#ReplayIdentity::canonical_bytes",
    "schema_dependencies=fs-obs:replay-identity-frame,fs-verify:fem1d-mms-class",
    "digest=fnv1a64",
    "encoding=typed-binary",
    "sources=MmsProblem,MmsProblemIdentityReceipt",
    "source_fields=MmsProblem.class:semantic,MmsProblem.mesh:semantic,MmsProblem.identity:derived:recomputed-from-class-and-mesh,MmsProblemIdentityReceipt.declared_identity_version:semantic,MmsProblemIdentityReceipt.canonical_bytes:semantic,MmsProblemIdentityReceipt.root:derived:validated-fnv-root-of-retained-canonical-bytes",
    "source_bindings=MmsProblem.class>class-child-frame-version-root+class-canonical-bytes,MmsProblem.mesh>mesh-f64-le,MmsProblemIdentityReceipt.declared_identity_version>retained-producer-version,MmsProblemIdentityReceipt.canonical_bytes>retained-canonical-bytes",
    "external_semantic_fields=artifact-domain,problem-schema-version",
    "semantic_fields=artifact-domain,problem-schema-version,class-child-frame-version-root,class-canonical-bytes,mesh-f64-le,retained-producer-version,retained-canonical-bytes",
    "excluded_fields=none",
    "consumers=MmsProblem::new,MmsProblem::from_class,MmsProblem::with_mesh,MmsProblemIdentityReceipt,MmsProblem::admit_identity_receipt",
    "mutations=artifact-domain:crates/fs-verify/src/fem1d.rs#mms_problem_identity_fields_move_independently,problem-schema-version:crates/fs-verify/src/fem1d.rs#mms_problem_identity_fields_move_independently,class-child-frame-version-root:crates/fs-verify/src/fem1d.rs#mms_problem_identity_fields_move_independently,class-canonical-bytes:crates/fs-verify/src/fem1d.rs#mms_problem_identity_fields_move_independently,mesh-f64-le:crates/fs-verify/src/fem1d.rs#mms_problem_identity_fields_move_independently,retained-producer-version:crates/fs-verify/src/fem1d.rs#mms_problem_identity_receipts_fail_closed,retained-canonical-bytes:crates/fs-verify/src/fem1d.rs#mms_problem_identity_receipts_fail_closed",
    "nonsemantic_mutations=none",
    "field_guard=classify_mms_problem_identity_fields",
    "transport_guard=MmsProblemIdentityReceipt::from_retained_parts",
    "version_guard=crates/fs-verify/src/fem1d.rs#mms_problem_identity_receipts_fail_closed",
    "coupling_surface=fs-verify:fem1d-mms-problem",
];

const NEWTON_RESIDUAL_TOLERANCE: f64 = 1.0e-10;

/// Structured failure at the bounded fem1d execution boundary.
#[derive(Debug, Clone, PartialEq)]
pub enum Fem1dError {
    /// A resource-driving count exceeds the synchronous v0 envelope.
    ResourceLimit {
        /// Stable resource name.
        resource: &'static str,
        /// Requested count.
        requested: usize,
        /// Maximum admitted count.
        limit: usize,
    },
    /// A bounded allocation failed before numerical work began.
    AllocationFailed {
        /// Stable allocation stage.
        stage: &'static str,
        /// Requested elements.
        requested: usize,
    },
    /// Conservative work accounting overflowed or exceeded its ceiling.
    WorkBudgetExceeded {
        /// Operation whose work was estimated.
        stage: &'static str,
        /// Estimated work, or `None` when checked multiplication overflowed.
        estimated: Option<usize>,
        /// Maximum admitted work.
        limit: usize,
    },
    /// A retained telemetry counter exhausted its fixed-width representation.
    CounterOverflow {
        /// Stable counter identity.
        counter: &'static str,
    },
    /// A scalar public input is not finite or outside its required range.
    InvalidScalar {
        /// Stable field name.
        field: &'static str,
        /// Stable refusal reason.
        reason: &'static str,
    },
    /// A polynomial is empty or outside the degree-five class.
    PolynomialCoefficientCount {
        /// `u`, `f`, or `big_f`.
        field: &'static str,
        /// Observed coefficient count.
        count: usize,
    },
    /// A polynomial coefficient is non-finite.
    NonFinitePolynomialCoefficient {
        /// `u`, `f`, or `big_f`.
        field: &'static str,
        /// Offending coefficient.
        index: usize,
    },
    /// The MMS exact-solution metadata is inconsistent with homogeneous BCs.
    ExactSolutionBoundary,
    /// Public derived polynomial data differs from the canonical value from `u`.
    DerivedPolynomialMismatch {
        /// `f` or `big_f`.
        field: &'static str,
    },
    /// Mesh endpoints are not bit-canonical `+0.0` and `1.0`.
    MeshDomain,
    /// A mesh node is non-finite.
    NonFiniteMeshNode {
        /// Offending node.
        index: usize,
    },
    /// A mesh cell is not strictly increasing.
    NonIncreasingMeshCell {
        /// Offending cell.
        cell: usize,
    },
    /// A positive mesh width has a non-finite reciprocal.
    NonFiniteReciprocal {
        /// Offending cell.
        cell: usize,
    },
    /// Nodal input length does not equal mesh length.
    CandidateLength {
        /// `candidate` or `start`.
        field: &'static str,
        /// Required length.
        expected: usize,
        /// Observed length.
        actual: usize,
    },
    /// A nodal public input is non-finite.
    NonFiniteCandidate {
        /// `candidate` or `start`.
        field: &'static str,
        /// Offending value.
        index: usize,
    },
    /// Nodal endpoints are not bit-canonical homogeneous `+0.0` values.
    CandidateBoundary {
        /// `candidate` or `start`.
        field: &'static str,
    },
    /// An internal tridiagonal system has inconsistent shapes.
    LinearSystemShape {
        /// Solver stage.
        stage: &'static str,
    },
    /// Thomas elimination encountered a zero or non-finite pivot.
    SingularPivot {
        /// Solver stage.
        stage: &'static str,
        /// Pivot row.
        row: usize,
    },
    /// A derived numerical value was non-finite.
    NonFiniteIntermediate {
        /// Stable computation stage.
        stage: &'static str,
        /// Element or row, when applicable.
        index: Option<usize>,
    },
    /// A caller required convergence but the bounded Newton run did not converge.
    NonConverged {
        /// `cold`, `warm`, or another stable caller role.
        stage: &'static str,
        /// Updates performed.
        iterations: u32,
        /// Final finite residual norm.
        residual_norm: f64,
    },
}

impl fmt::Display for Fem1dError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::ResourceLimit {
                resource,
                requested,
                limit,
            } => write!(
                f,
                "fem1d resource `{resource}` requested {requested}, limit is {limit}"
            ),
            Self::AllocationFailed { stage, requested } => {
                write!(
                    f,
                    "fem1d allocation failed during {stage} for {requested} elements"
                )
            }
            Self::WorkBudgetExceeded {
                stage,
                estimated,
                limit,
            } => write!(
                f,
                "fem1d work for {stage} was {estimated:?}, limit is {limit}"
            ),
            Self::CounterOverflow { counter } => {
                write!(f, "fem1d telemetry counter `{counter}` overflowed")
            }
            Self::InvalidScalar { field, reason } => {
                write!(f, "invalid fem1d scalar `{field}`: {reason}")
            }
            Self::PolynomialCoefficientCount { field, count } => write!(
                f,
                "fem1d polynomial `{field}` has {count} coefficients; expected 1..={MAX_FEM1D_POLY_COEFFICIENTS}"
            ),
            Self::NonFinitePolynomialCoefficient { field, index } => {
                write!(
                    f,
                    "fem1d polynomial `{field}` coefficient {index} is non-finite"
                )
            }
            Self::ExactSolutionBoundary => {
                f.write_str("fem1d exact-solution metadata is not homogeneous at both endpoints")
            }
            Self::DerivedPolynomialMismatch { field } => {
                write!(f, "fem1d derived polynomial `{field}` is not canonical")
            }
            Self::MeshDomain => f.write_str("fem1d mesh endpoints must be canonical +0.0 and 1.0"),
            Self::NonFiniteMeshNode { index } => {
                write!(f, "fem1d mesh node {index} is non-finite")
            }
            Self::NonIncreasingMeshCell { cell } => {
                write!(f, "fem1d mesh cell {cell} is not strictly increasing")
            }
            Self::NonFiniteReciprocal { cell } => {
                write!(
                    f,
                    "fem1d mesh cell {cell} has a non-finite reciprocal width"
                )
            }
            Self::CandidateLength {
                field,
                expected,
                actual,
            } => write!(
                f,
                "fem1d `{field}` length {actual} differs from mesh length {expected}"
            ),
            Self::NonFiniteCandidate { field, index } => {
                write!(f, "fem1d `{field}` value {index} is non-finite")
            }
            Self::CandidateBoundary { field } => {
                write!(
                    f,
                    "fem1d `{field}` endpoints must be canonical homogeneous +0.0"
                )
            }
            Self::LinearSystemShape { stage } => {
                write!(
                    f,
                    "fem1d {stage} tridiagonal system has inconsistent shapes"
                )
            }
            Self::SingularPivot { stage, row } => {
                write!(f, "fem1d {stage} pivot {row} is zero or non-finite")
            }
            Self::NonFiniteIntermediate { stage, index } => match index {
                Some(index) => write!(f, "fem1d {stage} is non-finite at index {index}"),
                None => write!(f, "fem1d {stage} is non-finite"),
            },
            Self::NonConverged {
                stage,
                iterations,
                residual_norm,
            } => write!(
                f,
                "fem1d {stage} did not converge after {iterations} updates (residual {residual_norm:.6e})"
            ),
        }
    }
}

impl std::error::Error for Fem1dError {}

/// Why retained MMS identity evidence was refused.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MmsIdentityAdmissionError {
    /// The producer-specific identity schema is stale or unknown.
    UnsupportedVersion {
        /// Identity surface being admitted.
        identity: &'static str,
        /// Version declared by retained state.
        declared: u32,
        /// Exact version supported by this build.
        supported: u32,
    },
    /// Retained canonical bytes exceed the owner-specific admission cap.
    CanonicalBytesExceeded {
        /// Identity surface being admitted.
        identity: &'static str,
        /// Bytes supplied by retained state.
        requested: usize,
        /// Maximum canonical bytes admitted by this build.
        limit: usize,
    },
    /// The retained root is not the FNV-1a root of the retained bytes.
    RootMismatch {
        /// Identity surface being admitted.
        identity: &'static str,
    },
    /// Self-consistent retained evidence does not equal the current object.
    IdentityMismatch {
        /// Identity surface being admitted.
        identity: &'static str,
    },
}

impl fmt::Display for MmsIdentityAdmissionError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnsupportedVersion {
                identity,
                declared,
                supported,
            } => write!(
                f,
                "retained {identity} identity v{declared} is unsupported; this build accepts exactly v{supported}"
            ),
            Self::CanonicalBytesExceeded {
                identity,
                requested,
                limit,
            } => write!(
                f,
                "retained {identity} identity has {requested} canonical bytes; limit is {limit}"
            ),
            Self::RootMismatch { identity } => {
                write!(
                    f,
                    "retained {identity} root does not match its canonical bytes"
                )
            }
            Self::IdentityMismatch { identity } => {
                write!(f, "retained {identity} does not match the admitted object")
            }
        }
    }
}

impl std::error::Error for MmsIdentityAdmissionError {}

fn check_mms_identity_version(
    identity: &'static str,
    declared: u32,
    supported: u32,
) -> Result<(), MmsIdentityAdmissionError> {
    if declared == supported {
        Ok(())
    } else {
        Err(MmsIdentityAdmissionError::UnsupportedVersion {
            identity,
            declared,
            supported,
        })
    }
}

/// Refuse retained MMS-class evidence from stale or future producer schemas.
///
/// # Errors
/// [`MmsIdentityAdmissionError::UnsupportedVersion`] unless `declared` is
/// exactly [`MMS_CLASS_IDENTITY_VERSION`].
pub fn check_mms_class_identity_version(declared: u32) -> Result<(), MmsIdentityAdmissionError> {
    check_mms_identity_version("MMS class", declared, MMS_CLASS_IDENTITY_VERSION)
}

/// Refuse retained MMS-problem evidence from stale or future producer schemas.
///
/// # Errors
/// [`MmsIdentityAdmissionError::UnsupportedVersion`] unless `declared` is
/// exactly [`MMS_PROBLEM_IDENTITY_VERSION`].
pub fn check_mms_problem_identity_version(declared: u32) -> Result<(), MmsIdentityAdmissionError> {
    check_mms_identity_version("MMS problem", declared, MMS_PROBLEM_IDENTITY_VERSION)
}

/// An immutable canonical polynomial in monomial coefficients
/// (`c[0] + c[1]x + …`).
///
/// Construction normalizes every signed zero to `+0.0` and removes trailing
/// zero coefficients while retaining one coefficient for the zero polynomial.
#[derive(Debug, Clone, PartialEq)]
pub struct Poly(Vec<f64>);

impl Poly {
    /// Construct one finite polynomial inside the degree-five MMS envelope.
    ///
    /// # Errors
    /// Returns [`Fem1dError`] for an empty, oversized, or non-finite input.
    pub fn new(coefficients: Vec<f64>) -> Result<Self, Fem1dError> {
        Self::from_coefficients(coefficients, "poly")
    }

    fn from_coefficients(
        mut coefficients: Vec<f64>,
        field: &'static str,
    ) -> Result<Self, Fem1dError> {
        if coefficients.is_empty() {
            return Err(Fem1dError::PolynomialCoefficientCount {
                field,
                count: coefficients.len(),
            });
        }
        if coefficients.len() > MAX_FEM1D_RAW_POLY_COEFFICIENTS {
            return Err(Fem1dError::ResourceLimit {
                resource: "raw polynomial coefficients",
                requested: coefficients.len(),
                limit: MAX_FEM1D_RAW_POLY_COEFFICIENTS,
            });
        }
        if let Some(index) = coefficients.iter().position(|value| !value.is_finite()) {
            return Err(Fem1dError::NonFinitePolynomialCoefficient { field, index });
        }
        for coefficient in &mut coefficients {
            *coefficient = canonicalize_zero(*coefficient);
        }
        while coefficients.len() > 1
            && coefficients
                .last()
                .is_some_and(|coefficient| *coefficient == 0.0)
        {
            coefficients.pop();
        }
        if coefficients.len() > MAX_FEM1D_POLY_COEFFICIENTS {
            return Err(Fem1dError::PolynomialCoefficientCount {
                field,
                count: coefficients.len(),
            });
        }
        Ok(Self(coefficients))
    }

    /// Canonical monomial coefficients.
    #[must_use]
    pub fn coefficients(&self) -> &[f64] {
        &self.0
    }

    fn try_copy(&self, stage: &'static str) -> Result<Self, Fem1dError> {
        let mut coefficients = Vec::new();
        coefficients
            .try_reserve_exact(self.0.len())
            .map_err(|_| Fem1dError::AllocationFailed {
                stage,
                requested: self.0.len(),
            })?;
        coefficients.extend_from_slice(&self.0);
        Ok(Self(coefficients))
    }

    /// Whether the finite binary64 coefficients sum to exactly zero.
    ///
    /// At `x = 1`, every monomial equals one, so this is the exact
    /// homogeneous-boundary predicate for the stored polynomial. A fixed-size
    /// binary superaccumulator compares positive and negative magnitudes in
    /// units of `2^-1074`; ordinary or interval Horner evaluation is not enough
    /// because rounding can respectively hide a residue or merely enclose zero.
    #[must_use]
    pub fn is_exactly_zero_at_one(&self) -> bool {
        const LIMBS: usize = 34;
        const FRACTION_MASK: u64 = (1_u64 << 52) - 1;
        const EXPONENT_MASK: u64 = 0x7ff;

        fn add_word(accumulator: &mut [u64; LIMBS], mut index: usize, word: u64) {
            let (sum, mut carry) = accumulator[index].overflowing_add(word);
            accumulator[index] = sum;
            index += 1;
            while carry && index < accumulator.len() {
                let (sum, next_carry) = accumulator[index].overflowing_add(1);
                accumulator[index] = sum;
                carry = next_carry;
                index += 1;
            }
        }

        fn add_significand(accumulator: &mut [u64; LIMBS], significand: u64, shift: usize) {
            let limb = shift / 64;
            let offset = shift % 64;
            add_word(accumulator, limb, significand << offset);
            if offset != 0 {
                add_word(accumulator, limb + 1, significand >> (64 - offset));
            }
        }

        if self.0.is_empty() || self.0.len() > MAX_FEM1D_POLY_COEFFICIENTS {
            return false;
        }
        let mut positive = [0_u64; LIMBS];
        let mut negative = [0_u64; LIMBS];
        for coefficient in &self.0 {
            let bits = coefficient.to_bits();
            let exponent = ((bits >> 52) & EXPONENT_MASK) as usize;
            if exponent == EXPONENT_MASK as usize {
                return false;
            }
            let fraction = bits & FRACTION_MASK;
            let significand = if exponent == 0 {
                fraction
            } else {
                fraction | (1_u64 << 52)
            };
            if significand == 0 {
                continue;
            }
            // In units of 2^-1074, subnormal significands start at bit zero
            // and a normal value with biased exponent e starts at bit e - 1.
            let shift = exponent.saturating_sub(1);
            let accumulator = if bits >> 63 == 0 {
                &mut positive
            } else {
                &mut negative
            };
            add_significand(accumulator, significand, shift);
        }
        positive == negative
    }

    /// Derivative.
    pub fn derive(&self) -> Result<Poly, Fem1dError> {
        if self.0.len() <= 1 {
            let mut derived = Vec::new();
            derived
                .try_reserve_exact(1)
                .map_err(|_| Fem1dError::AllocationFailed {
                    stage: "polynomial derivative",
                    requested: 1,
                })?;
            derived.push(0.0);
            return Poly::from_coefficients(derived, "derived polynomial");
        }
        let mut derived = Vec::new();
        derived
            .try_reserve_exact(self.0.len() - 1)
            .map_err(|_| Fem1dError::AllocationFailed {
                stage: "polynomial derivative",
                requested: self.0.len() - 1,
            })?;
        for (k, &coefficient) in self.0[1..].iter().enumerate() {
            let value = coefficient * (k + 1) as f64;
            if !value.is_finite() {
                return Err(Fem1dError::NonFiniteIntermediate {
                    stage: "polynomial derivative",
                    index: Some(k),
                });
            }
            derived.push(value);
        }
        Poly::from_coefficients(derived, "derived polynomial")
    }

    /// Antiderivative with zero constant term.
    pub fn antiderive(&self) -> Result<Poly, Fem1dError> {
        let requested =
            self.0
                .len()
                .checked_add(1)
                .ok_or(Fem1dError::PolynomialCoefficientCount {
                    field: "antiderivative",
                    count: usize::MAX,
                })?;
        if requested > MAX_FEM1D_POLY_COEFFICIENTS {
            return Err(Fem1dError::PolynomialCoefficientCount {
                field: "antiderivative",
                count: requested,
            });
        }
        let mut out = Vec::new();
        out.try_reserve_exact(requested)
            .map_err(|_| Fem1dError::AllocationFailed {
                stage: "polynomial antiderivative",
                requested,
            })?;
        out.push(0.0);
        for (k, &coefficient) in self.0.iter().enumerate() {
            let value = coefficient / (k + 1) as f64;
            if !value.is_finite() {
                return Err(Fem1dError::NonFiniteIntermediate {
                    stage: "polynomial antiderivative",
                    index: Some(k),
                });
            }
            out.push(value);
        }
        Poly::from_coefficients(out, "antiderivative")
    }

    /// Negate this canonical polynomial without another allocation.
    #[must_use]
    #[allow(clippy::should_implement_trait)] // consuming canonicalization helper, not algebraic `Neg`
    pub fn neg(mut self) -> Poly {
        for coefficient in &mut self.0 {
            *coefficient = canonicalize_zero(-*coefficient);
        }
        // Negating in place preserves finiteness, size, and canonical shape.
        self
    }

    /// Horner evaluation (f64).
    #[must_use]
    pub fn eval(&self, x: f64) -> f64 {
        self.0.iter().rev().fold(0.0, |acc, &c| acc * x + c)
    }

    /// Horner evaluation with outward-rounded intervals.
    #[must_use]
    pub fn eval_iv(&self, x: Iv) -> Iv {
        self.0
            .iter()
            .rev()
            .fold(Iv::zero(), |acc, &c| acc.mul(x).add(Iv::point(c)))
    }

    /// Polynomial degree (0 for constants).
    #[must_use]
    pub fn degree(&self) -> usize {
        self.0.len().saturating_sub(1)
    }
}

/// Retained, self-consistent identity evidence for one admitted MMS class.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MmsClassIdentityReceipt {
    declared_identity_version: u32,
    canonical_bytes: Vec<u8>,
    root: u64,
}

impl MmsClassIdentityReceipt {
    /// Admit retained receipt parts after enforcing the class byte cap and
    /// verifying that `root` authenticates the exact supplied bytes.
    ///
    /// # Errors
    /// [`MmsIdentityAdmissionError`] when the producer version is not current,
    /// the bytes exceed the class cap, or the supplied root is not
    /// self-consistent. Version refusal precedes byte hashing.
    pub fn from_retained_parts(
        declared_identity_version: u32,
        canonical_bytes: Vec<u8>,
        root: u64,
    ) -> Result<Self, MmsIdentityAdmissionError> {
        check_mms_class_identity_version(declared_identity_version)?;
        validate_retained_identity_parts(
            "MMS class",
            &canonical_bytes,
            root,
            MAX_FEM1D_CLASS_CANONICAL_IDENTITY_BYTES,
        )?;
        Ok(Self {
            declared_identity_version,
            canonical_bytes,
            root,
        })
    }

    /// Producer-specific schema declared by retained evidence.
    #[must_use]
    pub const fn declared_identity_version(&self) -> u32 {
        self.declared_identity_version
    }

    /// Exact retained canonical bytes.
    #[must_use]
    pub fn canonical_bytes(&self) -> &[u8] {
        &self.canonical_bytes
    }

    /// Retained FNV-1a root over [`Self::canonical_bytes`].
    #[must_use]
    pub const fn root(&self) -> u64 {
        self.root
    }
}

/// Retained, self-consistent identity evidence for one meshed MMS problem.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MmsProblemIdentityReceipt {
    declared_identity_version: u32,
    canonical_bytes: Vec<u8>,
    root: u64,
}

impl MmsProblemIdentityReceipt {
    /// Admit retained receipt parts after enforcing the problem byte cap and
    /// verifying that `root` authenticates the exact supplied bytes.
    ///
    /// # Errors
    /// [`MmsIdentityAdmissionError`] when the producer version is not current,
    /// the bytes exceed the problem cap, or the supplied root is not
    /// self-consistent. Version refusal precedes byte hashing.
    pub fn from_retained_parts(
        declared_identity_version: u32,
        canonical_bytes: Vec<u8>,
        root: u64,
    ) -> Result<Self, MmsIdentityAdmissionError> {
        check_mms_problem_identity_version(declared_identity_version)?;
        validate_retained_identity_parts(
            "MMS problem",
            &canonical_bytes,
            root,
            MAX_FEM1D_PROBLEM_CANONICAL_IDENTITY_BYTES,
        )?;
        Ok(Self {
            declared_identity_version,
            canonical_bytes,
            root,
        })
    }

    /// Producer-specific schema declared by retained evidence.
    #[must_use]
    pub const fn declared_identity_version(&self) -> u32 {
        self.declared_identity_version
    }

    /// Exact retained canonical bytes.
    #[must_use]
    pub fn canonical_bytes(&self) -> &[u8] {
        &self.canonical_bytes
    }

    /// Retained FNV-1a root over [`Self::canonical_bytes`].
    #[must_use]
    pub const fn root(&self) -> u64 {
        self.root
    }
}

fn validate_retained_identity_parts(
    identity: &'static str,
    canonical_bytes: &[u8],
    root: u64,
    limit: usize,
) -> Result<(), MmsIdentityAdmissionError> {
    if canonical_bytes.len() > limit {
        return Err(MmsIdentityAdmissionError::CanonicalBytesExceeded {
            identity,
            requested: canonical_bytes.len(),
            limit,
        });
    }
    if fs_obs::fnv1a64(canonical_bytes) != root {
        return Err(MmsIdentityAdmissionError::RootMismatch { identity });
    }
    Ok(())
}

fn copy_identity_bytes(bytes: &[u8], stage: &'static str) -> Result<Vec<u8>, Fem1dError> {
    let mut retained = Vec::new();
    retained
        .try_reserve_exact(bytes.len())
        .map_err(|_| Fem1dError::AllocationFailed {
            stage,
            requested: bytes.len(),
        })?;
    retained.extend_from_slice(bytes);
    Ok(retained)
}

/// One admitted manufactured-solution class, independent of discretization.
///
/// The exact solution is canonical and immutable. The forcing and its rounded
/// zero-constant antiderivative are derived once during construction and cannot
/// drift independently. The replay identity binds all stored semantic fields.
#[derive(Debug, Clone, PartialEq)]
pub struct MmsClass {
    name: String,
    u: Poly,
    forcing: Poly,
    rounded_forcing_antiderivative: Poly,
    identity: ReplayIdentity,
}

impl MmsClass {
    /// Construct and admit one canonical degree-five manufactured class.
    ///
    /// # Errors
    /// Returns [`Fem1dError`] for invalid identity text, non-homogeneous
    /// endpoint data, or non-finite derived arithmetic.
    pub fn new(name: &str, u: Poly) -> Result<Self, Fem1dError> {
        validate_identity(name, "problem.name")?;
        if u.0[0].to_bits() != 0.0_f64.to_bits() || !u.is_exactly_zero_at_one() {
            return Err(Fem1dError::ExactSolutionBoundary);
        }
        let forcing = u.derive()?.derive()?.neg();
        let rounded_forcing_antiderivative = forcing.antiderive()?;

        let mut owned_name = String::new();
        owned_name
            .try_reserve_exact(name.len())
            .map_err(|_| Fem1dError::AllocationFailed {
                stage: "manufactured class name",
                requested: name.len(),
            })?;
        owned_name.push_str(name);
        let identity = class_identity(&owned_name, &u, &forcing, &rounded_forcing_antiderivative)?;
        Ok(Self {
            name: owned_name,
            u,
            forcing,
            rounded_forcing_antiderivative,
            identity,
        })
    }

    /// Stable kernel/problem-class name.
    #[must_use]
    pub fn name(&self) -> &str {
        &self.name
    }

    /// Canonical exact-solution polynomial.
    #[must_use]
    pub fn exact_solution(&self) -> &Poly {
        &self.u
    }

    /// Canonical forcing `f = -u''`.
    #[must_use]
    pub fn forcing(&self) -> &Poly {
        &self.forcing
    }

    /// Rounded coefficients of the zero-constant forcing antiderivative.
    ///
    /// These bytes are replay/tightness metadata. The certified estimator still
    /// intervalizes exact coefficient division instead of treating them as an
    /// exact mathematical antiderivative.
    #[must_use]
    pub fn rounded_forcing_antiderivative(&self) -> &Poly {
        &self.rounded_forcing_antiderivative
    }

    /// Versioned canonical replay identity for the class.
    #[must_use]
    pub fn identity(&self) -> &ReplayIdentity {
        &self.identity
    }

    /// Exact canonical bytes bound by [`Self::identity`].
    #[must_use]
    pub fn canonical_bytes(&self) -> &[u8] {
        self.identity.canonical_bytes()
    }

    /// Fallibly capture retained identity evidence for this exact class.
    ///
    /// # Errors
    /// [`Fem1dError::AllocationFailed`] if the canonical-byte copy cannot be
    /// reserved before any receipt escapes.
    pub fn identity_receipt(&self) -> Result<MmsClassIdentityReceipt, Fem1dError> {
        let canonical_bytes = copy_identity_bytes(
            self.identity.canonical_bytes(),
            "MMS class identity receipt",
        )?;
        Ok(MmsClassIdentityReceipt {
            declared_identity_version: MMS_CLASS_IDENTITY_VERSION,
            canonical_bytes,
            root: self.identity.root(),
        })
    }

    /// Admit a retained class receipt only when its producer version, exact
    /// canonical bytes, and root equal this admitted class.
    ///
    /// # Errors
    /// [`MmsIdentityAdmissionError`] for stale/future versions, hostile sizes,
    /// inconsistent roots, or a different class identity.
    pub fn admit_identity_receipt(
        &self,
        receipt: &MmsClassIdentityReceipt,
    ) -> Result<(), MmsIdentityAdmissionError> {
        check_mms_class_identity_version(receipt.declared_identity_version)?;
        validate_retained_identity_parts(
            "MMS class",
            &receipt.canonical_bytes,
            receipt.root,
            MAX_FEM1D_CLASS_CANONICAL_IDENTITY_BYTES,
        )?;
        if receipt.root != self.identity.root()
            || receipt.canonical_bytes != self.identity.canonical_bytes()
        {
            return Err(MmsIdentityAdmissionError::IdentityMismatch {
                identity: "MMS class",
            });
        }
        Ok(())
    }
}

/// A canonical manufactured problem: one admitted class on one admitted mesh.
#[derive(Debug, Clone, PartialEq)]
pub struct MmsProblem {
    class: MmsClass,
    mesh: Vec<f64>,
    identity: ReplayIdentity,
}

impl MmsProblem {
    /// Build from an exact polynomial solution and a mesh.
    ///
    /// # Errors
    /// Returns [`Fem1dError`] when either the class or mesh is inadmissible.
    pub fn new(name: &str, u: Poly, mesh: Vec<f64>) -> Result<Self, Fem1dError> {
        let mesh = admit_mesh(mesh)?;
        Self::from_admitted_parts(MmsClass::new(name, u)?, mesh)
    }

    /// Place an admitted class on an admitted canonical mesh.
    ///
    /// # Errors
    /// Returns [`Fem1dError`] for invalid mesh size, coordinates, ordering, or
    /// reciprocal widths.
    pub fn from_class(class: MmsClass, mesh: Vec<f64>) -> Result<Self, Fem1dError> {
        let mesh = admit_mesh(mesh)?;
        Self::from_admitted_parts(class, mesh)
    }

    fn from_admitted_parts(class: MmsClass, mesh: Vec<f64>) -> Result<Self, Fem1dError> {
        let identity = problem_identity(&class, &mesh)?;
        Ok(Self {
            class,
            mesh,
            identity,
        })
    }

    /// Reuse the same admitted class on a new mesh.
    ///
    /// # Errors
    /// Returns [`Fem1dError`] when the new mesh is inadmissible or rebuilding
    /// the immutable class/problem identity cannot reserve bounded storage.
    pub fn with_mesh(&self, mesh: Vec<f64>) -> Result<Self, Fem1dError> {
        let mesh = admit_mesh(mesh)?;
        let exact_solution = self
            .class
            .exact_solution()
            .try_copy("MMS exact solution reuse")?;
        let class = MmsClass::new(self.class.name(), exact_solution)?;
        Self::from_admitted_parts(class, mesh)
    }

    /// Admitted manufactured-solution class.
    #[must_use]
    pub fn class(&self) -> &MmsClass {
        &self.class
    }

    /// Stable problem name.
    #[must_use]
    pub fn name(&self) -> &str {
        self.class.name()
    }

    /// Canonical exact-solution polynomial.
    #[must_use]
    pub fn exact_solution(&self) -> &Poly {
        self.class.exact_solution()
    }

    /// Canonical forcing polynomial.
    #[must_use]
    pub fn forcing(&self) -> &Poly {
        self.class.forcing()
    }

    /// Rounded forcing-antiderivative replay metadata.
    #[must_use]
    pub fn rounded_forcing_antiderivative(&self) -> &Poly {
        self.class.rounded_forcing_antiderivative()
    }

    /// Canonical mesh nodes.
    #[must_use]
    pub fn mesh(&self) -> &[f64] {
        &self.mesh
    }

    /// Versioned canonical replay identity for class plus mesh.
    #[must_use]
    pub fn identity(&self) -> &ReplayIdentity {
        &self.identity
    }

    /// Exact canonical bytes bound by [`Self::identity`].
    #[must_use]
    pub fn canonical_bytes(&self) -> &[u8] {
        self.identity.canonical_bytes()
    }

    /// Fallibly capture retained identity evidence for this exact problem.
    ///
    /// # Errors
    /// [`Fem1dError::AllocationFailed`] if the canonical-byte copy cannot be
    /// reserved before any receipt escapes.
    pub fn identity_receipt(&self) -> Result<MmsProblemIdentityReceipt, Fem1dError> {
        let canonical_bytes = copy_identity_bytes(
            self.identity.canonical_bytes(),
            "MMS problem identity receipt",
        )?;
        Ok(MmsProblemIdentityReceipt {
            declared_identity_version: MMS_PROBLEM_IDENTITY_VERSION,
            canonical_bytes,
            root: self.identity.root(),
        })
    }

    /// Admit a retained problem receipt only when its producer version, exact
    /// canonical bytes, and root equal this admitted problem.
    ///
    /// # Errors
    /// [`MmsIdentityAdmissionError`] for stale/future versions, hostile sizes,
    /// inconsistent roots, or a different problem identity.
    pub fn admit_identity_receipt(
        &self,
        receipt: &MmsProblemIdentityReceipt,
    ) -> Result<(), MmsIdentityAdmissionError> {
        check_mms_problem_identity_version(receipt.declared_identity_version)?;
        validate_retained_identity_parts(
            "MMS problem",
            &receipt.canonical_bytes,
            receipt.root,
            MAX_FEM1D_PROBLEM_CANONICAL_IDENTITY_BYTES,
        )?;
        if receipt.root != self.identity.root()
            || receipt.canonical_bytes != self.identity.canonical_bytes()
        {
            return Err(MmsIdentityAdmissionError::IdentityMismatch {
                identity: "MMS problem",
            });
        }
        Ok(())
    }
}

#[allow(dead_code)]
fn classify_mms_class_identity_fields(class: &MmsClass, receipt: &MmsClassIdentityReceipt) {
    let MmsClass {
        name: _,
        u: _,
        forcing: _,
        rounded_forcing_antiderivative: _,
        identity: _,
    } = class;
    let MmsClassIdentityReceipt {
        declared_identity_version: _,
        canonical_bytes: _,
        root: _,
    } = receipt;
}

#[allow(dead_code)]
fn classify_mms_problem_identity_fields(problem: &MmsProblem, receipt: &MmsProblemIdentityReceipt) {
    let MmsProblem {
        class: _,
        mesh: _,
        identity: _,
    } = problem;
    let MmsProblemIdentityReceipt {
        declared_identity_version: _,
        canonical_bytes: _,
        root: _,
    } = receipt;
}

fn canonicalize_zero(value: f64) -> f64 {
    if value == 0.0 { 0.0 } else { value }
}

fn admit_mesh(mut mesh: Vec<f64>) -> Result<Vec<f64>, Fem1dError> {
    validate_mesh_size(&mesh)?;
    for node in &mut mesh {
        *node = canonicalize_zero(*node);
    }
    validate_mesh(&mesh)?;
    Ok(mesh)
}

fn f64_bytes(values: &[f64], stage: &'static str) -> Result<Vec<u8>, Fem1dError> {
    let requested =
        values
            .len()
            .checked_mul(size_of::<u64>())
            .ok_or(Fem1dError::AllocationFailed {
                stage,
                requested: usize::MAX,
            })?;
    let mut bytes = Vec::new();
    bytes
        .try_reserve_exact(requested)
        .map_err(|_| Fem1dError::AllocationFailed { stage, requested })?;
    for value in values {
        bytes.extend_from_slice(&value.to_bits().to_le_bytes());
    }
    Ok(bytes)
}

fn identity_build_error(
    error: IdentityBuildError,
    stage: &'static str,
    resource: &'static str,
    limit: usize,
) -> Fem1dError {
    match error {
        IdentityBuildError::CanonicalBytesExceeded { requested, .. }
        | IdentityBuildError::FramedLengthNotRepresentable { length: requested } => {
            Fem1dError::ResourceLimit {
                resource,
                requested,
                limit,
            }
        }
        IdentityBuildError::CanonicalLengthOverflow => Fem1dError::ResourceLimit {
            resource,
            requested: usize::MAX,
            limit,
        },
        IdentityBuildError::AllocationFailed { requested } => {
            Fem1dError::AllocationFailed { stage, requested }
        }
    }
}

fn class_identity_with_budget(
    name: &str,
    exact_solution: &Poly,
    forcing: &Poly,
    rounded_forcing_antiderivative: &Poly,
    max_canonical_bytes: usize,
) -> Result<ReplayIdentity, Fem1dError> {
    let exact_solution_bytes = f64_bytes(&exact_solution.0, "MMS exact-solution identity")?;
    let forcing_bytes = f64_bytes(&forcing.0, "MMS forcing identity")?;
    let antiderivative_bytes = f64_bytes(
        &rounded_forcing_antiderivative.0,
        "MMS forcing-antiderivative identity",
    )?;
    class_identity_from_parts_with_budget(
        MMS_CLASS_IDENTITY_DOMAIN,
        MMS_CLASS_IDENTITY_VERSION,
        name,
        &exact_solution_bytes,
        &forcing_bytes,
        &antiderivative_bytes,
        max_canonical_bytes,
    )
}

#[allow(clippy::too_many_arguments)]
fn class_identity_from_parts_with_budget(
    domain: &str,
    producer_version: u32,
    name: &str,
    exact_solution_bytes: &[u8],
    forcing_bytes: &[u8],
    antiderivative_bytes: &[u8],
    max_canonical_bytes: usize,
) -> Result<ReplayIdentity, Fem1dError> {
    let identity = (|| -> Result<ReplayIdentity, IdentityBuildError> {
        let builder = BoundedIdentityBuilder::new(domain, max_canonical_bytes)?;
        let builder = builder.u64("class_schema", u64::from(producer_version))?;
        let builder = builder.str("name", name)?;
        let builder = builder.bytes("exact_solution_f64_le", exact_solution_bytes)?;
        let builder = builder.bytes("forcing_f64_le", forcing_bytes)?;
        let builder = builder.bytes(
            "rounded_forcing_antiderivative_f64_le",
            antiderivative_bytes,
        )?;
        Ok(builder.finish())
    })()
    .map_err(|error| {
        identity_build_error(
            error,
            "MMS class replay identity",
            "MMS class canonical identity bytes",
            max_canonical_bytes,
        )
    })?;
    Ok(identity)
}

fn class_identity(
    name: &str,
    exact_solution: &Poly,
    forcing: &Poly,
    rounded_forcing_antiderivative: &Poly,
) -> Result<ReplayIdentity, Fem1dError> {
    class_identity_with_budget(
        name,
        exact_solution,
        forcing,
        rounded_forcing_antiderivative,
        MAX_FEM1D_CLASS_CANONICAL_IDENTITY_BYTES,
    )
}

fn problem_identity_with_budget(
    class: &MmsClass,
    mesh: &[f64],
    max_canonical_bytes: usize,
) -> Result<ReplayIdentity, Fem1dError> {
    let mesh_bytes = f64_bytes(mesh, "MMS problem mesh identity")?;
    problem_identity_from_parts_with_budget(
        MMS_PROBLEM_IDENTITY_DOMAIN,
        MMS_PROBLEM_IDENTITY_VERSION,
        class.identity(),
        class.canonical_bytes(),
        &mesh_bytes,
        max_canonical_bytes,
    )
}

fn problem_identity_from_parts_with_budget(
    domain: &str,
    producer_version: u32,
    class_identity: &ReplayIdentity,
    class_canonical_bytes: &[u8],
    mesh_bytes: &[u8],
    max_canonical_bytes: usize,
) -> Result<ReplayIdentity, Fem1dError> {
    let identity = (|| -> Result<ReplayIdentity, IdentityBuildError> {
        let builder = BoundedIdentityBuilder::new(domain, max_canonical_bytes)?;
        let builder = builder.u64("problem_schema", u64::from(producer_version))?;
        let builder = builder.child("class", class_identity)?;
        let builder = builder.bytes("class_canonical_bytes", class_canonical_bytes)?;
        let builder = builder.bytes("mesh_f64_le", mesh_bytes)?;
        Ok(builder.finish())
    })()
    .map_err(|error| {
        identity_build_error(
            error,
            "MMS problem replay identity",
            "MMS problem canonical identity bytes",
            max_canonical_bytes,
        )
    })?;
    Ok(identity)
}

fn problem_identity(class: &MmsClass, mesh: &[f64]) -> Result<ReplayIdentity, Fem1dError> {
    problem_identity_with_budget(class, mesh, MAX_FEM1D_PROBLEM_CANONICAL_IDENTITY_BYTES)
}

fn poly_bits_equal(left: &Poly, right: &Poly) -> bool {
    left.0.len() == right.0.len()
        && left
            .0
            .iter()
            .zip(&right.0)
            .all(|(left, right)| left.to_bits() == right.to_bits())
}

/// Validate the complete v0 MMS problem before allocation or numerical work.
pub(crate) fn validate_problem(problem: &MmsProblem) -> Result<(), Fem1dError> {
    validate_identity(problem.name(), "problem.name")?;
    validate_mesh(problem.mesh())?;
    for (field, polynomial) in [
        ("u", problem.exact_solution()),
        ("f", problem.forcing()),
        ("big_f", problem.rounded_forcing_antiderivative()),
    ] {
        if !(1..=MAX_FEM1D_POLY_COEFFICIENTS).contains(&polynomial.0.len()) {
            return Err(Fem1dError::PolynomialCoefficientCount {
                field,
                count: polynomial.0.len(),
            });
        }
        if let Some(index) = polynomial.0.iter().position(|value| !value.is_finite()) {
            return Err(Fem1dError::NonFinitePolynomialCoefficient { field, index });
        }
    }
    if problem.exact_solution().0[0].to_bits() != 0.0_f64.to_bits()
        || !problem.exact_solution().is_exactly_zero_at_one()
    {
        return Err(Fem1dError::ExactSolutionBoundary);
    }
    let expected_f = problem.exact_solution().derive()?.derive()?.neg();
    if expected_f.0.iter().any(|value| !value.is_finite()) {
        return Err(Fem1dError::NonFiniteIntermediate {
            stage: "canonical forcing",
            index: None,
        });
    }
    if !poly_bits_equal(problem.forcing(), &expected_f) {
        return Err(Fem1dError::DerivedPolynomialMismatch { field: "f" });
    }
    let expected_big_f = expected_f.antiderive()?;
    if expected_big_f.0.iter().any(|value| !value.is_finite()) {
        return Err(Fem1dError::NonFiniteIntermediate {
            stage: "canonical forcing antiderivative",
            index: None,
        });
    }
    if !poly_bits_equal(problem.rounded_forcing_antiderivative(), &expected_big_f) {
        return Err(Fem1dError::DerivedPolynomialMismatch { field: "big_f" });
    }
    Ok(())
}

fn validate_mesh(mesh: &[f64]) -> Result<(), Fem1dError> {
    validate_mesh_size(mesh)?;
    if mesh[0].to_bits() != 0.0_f64.to_bits() || mesh[mesh.len() - 1].to_bits() != 1.0_f64.to_bits()
    {
        return Err(Fem1dError::MeshDomain);
    }
    if let Some(index) = mesh.iter().position(|value| !value.is_finite()) {
        return Err(Fem1dError::NonFiniteMeshNode { index });
    }
    for (cell, nodes) in mesh.windows(2).enumerate() {
        if nodes[0] >= nodes[1] {
            return Err(Fem1dError::NonIncreasingMeshCell { cell });
        }
        let reciprocal = 1.0 / (nodes[1] - nodes[0]);
        if !reciprocal.is_finite() {
            return Err(Fem1dError::NonFiniteReciprocal { cell });
        }
    }
    Ok(())
}

fn validate_mesh_size(mesh: &[f64]) -> Result<(), Fem1dError> {
    if !(2..=MAX_FEM1D_MESH_NODES).contains(&mesh.len()) {
        return Err(Fem1dError::ResourceLimit {
            resource: "mesh nodes",
            requested: mesh.len(),
            limit: MAX_FEM1D_MESH_NODES,
        });
    }
    Ok(())
}

pub(crate) fn validate_identity(value: &str, field: &'static str) -> Result<(), Fem1dError> {
    if value.is_empty() {
        return Err(Fem1dError::InvalidScalar {
            field,
            reason: "must not be empty",
        });
    }
    if value.len() > MAX_FEM1D_CLASS_NAME_BYTES {
        return Err(Fem1dError::ResourceLimit {
            resource: field,
            requested: value.len(),
            limit: MAX_FEM1D_CLASS_NAME_BYTES,
        });
    }
    if value.chars().any(char::is_control) {
        return Err(Fem1dError::InvalidScalar {
            field,
            reason: "must not contain control characters",
        });
    }
    Ok(())
}

/// Validate one conforming nodal vector against an already validated problem.
pub(crate) fn validate_candidate(
    problem: &MmsProblem,
    candidate: &[f64],
    field: &'static str,
) -> Result<(), Fem1dError> {
    if candidate.len() != problem.mesh().len() {
        return Err(Fem1dError::CandidateLength {
            field,
            expected: problem.mesh().len(),
            actual: candidate.len(),
        });
    }
    if let Some(index) = candidate.iter().position(|value| !value.is_finite()) {
        return Err(Fem1dError::NonFiniteCandidate { field, index });
    }
    if candidate[0].to_bits() != 0.0_f64.to_bits()
        || candidate[candidate.len() - 1].to_bits() != 0.0_f64.to_bits()
    {
        return Err(Fem1dError::CandidateBoundary { field });
    }
    Ok(())
}

pub(crate) fn validate_tolerance(tolerance: f64) -> Result<(), Fem1dError> {
    if !tolerance.is_finite() {
        Err(Fem1dError::InvalidScalar {
            field: "tolerance",
            reason: "must be finite",
        })
    } else if tolerance <= 0.0 {
        Err(Fem1dError::InvalidScalar {
            field: "tolerance",
            reason: "must be strictly positive",
        })
    } else {
        Ok(())
    }
}

pub(crate) fn validate_finite_scalar(value: f64, field: &'static str) -> Result<(), Fem1dError> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(Fem1dError::InvalidScalar {
            field,
            reason: "must be finite",
        })
    }
}

fn validate_work(stage: &'static str, factors: &[usize]) -> Result<(), Fem1dError> {
    let estimated = factors
        .iter()
        .try_fold(1usize, |work, factor| work.checked_mul(*factor));
    if estimated.is_none_or(|work| work > MAX_FEM1D_WORK_UNITS) {
        Err(Fem1dError::WorkBudgetExceeded {
            stage,
            estimated,
            limit: MAX_FEM1D_WORK_UNITS,
        })
    } else {
        Ok(())
    }
}

pub(crate) fn try_zeroed(stage: &'static str, length: usize) -> Result<Vec<f64>, Fem1dError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(length)
        .map_err(|_| Fem1dError::AllocationFailed {
            stage,
            requested: length,
        })?;
    values.resize(length, 0.0);
    Ok(values)
}

fn try_copy_slice(stage: &'static str, values: &[f64]) -> Result<Vec<f64>, Fem1dError> {
    let mut copy = Vec::new();
    copy.try_reserve_exact(values.len())
        .map_err(|_| Fem1dError::AllocationFailed {
            stage,
            requested: values.len(),
        })?;
    copy.extend_from_slice(values);
    Ok(copy)
}

fn non_finite(stage: &'static str, index: Option<usize>) -> Fem1dError {
    Fem1dError::NonFiniteIntermediate { stage, index }
}

fn add_finite(
    slot: &mut f64,
    delta: f64,
    stage: &'static str,
    index: usize,
) -> Result<(), Fem1dError> {
    let updated = *slot + delta;
    if !delta.is_finite() || !updated.is_finite() {
        return Err(non_finite(stage, Some(index)));
    }
    *slot = updated;
    Ok(())
}

fn solve_tridiagonal_into(
    off: &[f64],
    diag: &[f64],
    rhs: &[f64],
    c: &mut [f64],
    d: &mut [f64],
    solution: &mut [f64],
    stage: &'static str,
) -> Result<(), Fem1dError> {
    let n = diag.len();
    if rhs.len() != n
        || c.len() != n
        || d.len() != n
        || solution.len() != n
        || off.len() != n.saturating_sub(1)
    {
        return Err(Fem1dError::LinearSystemShape { stage });
    }
    if n == 0 {
        return Ok(());
    }
    let first_pivot = diag[0];
    if !first_pivot.is_finite() || first_pivot == 0.0 {
        return Err(Fem1dError::SingularPivot { stage, row: 0 });
    }
    c[0] = if n > 1 { off[0] / first_pivot } else { 0.0 };
    d[0] = rhs[0] / first_pivot;
    if !c[0].is_finite() || !d[0].is_finite() {
        return Err(non_finite("Thomas first row", Some(0)));
    }
    for row in 1..n {
        let pivot = diag[row] - off[row - 1] * c[row - 1];
        if !pivot.is_finite() || pivot == 0.0 {
            return Err(Fem1dError::SingularPivot { stage, row });
        }
        if row < n - 1 {
            c[row] = off[row] / pivot;
            if !c[row].is_finite() {
                return Err(non_finite("Thomas upper factor", Some(row)));
            }
        }
        d[row] = (rhs[row] - off[row - 1] * d[row - 1]) / pivot;
        if !d[row].is_finite() {
            return Err(non_finite("Thomas forward substitution", Some(row)));
        }
    }
    solution[n - 1] = d[n - 1];
    for row in (0..n - 1).rev() {
        solution[row] = d[row] - c[row] * solution[row + 1];
        if !solution[row].is_finite() {
            return Err(non_finite("Thomas back substitution", Some(row)));
        }
    }
    Ok(())
}

/// Solve `−u″ = f` with P1 elements on the problem's mesh (Thomas
/// algorithm; deterministic). Returns interior+boundary nodal values.
///
/// # Errors
/// Returns [`Fem1dError`] before issuing a solution when the public problem,
/// bounded work envelope, allocation, assembly, or elimination is unusable.
pub fn solve_p1(problem: &MmsProblem) -> Result<Vec<f64>, Fem1dError> {
    validate_problem(problem)?;
    let m = problem.mesh();
    let n = m.len();
    let interior = n - 2;
    validate_work(
        "P1 solve",
        &[
            n - 1,
            problem
                .forcing()
                .coefficients()
                .len()
                .saturating_mul(5)
                .saturating_add(16),
        ],
    )?;
    if interior == 0 {
        return try_zeroed("P1 boundary-only solution", n);
    }
    // Tridiagonal stiffness + exact load via 5-pt Gauss per element.
    let mut diag = try_zeroed("P1 diagonal", interior)?;
    let mut off = try_zeroed("P1 off-diagonal", interior.saturating_sub(1))?;
    let mut rhs = try_zeroed("P1 right-hand side", interior)?;
    for e in 0..n - 1 {
        let (x0, x1) = (m[e], m[e + 1]);
        let h = x1 - x0;
        let k = 1.0 / h;
        // Assemble stiffness.
        if e >= 1 {
            add_finite(&mut diag[e - 1], k, "P1 diagonal assembly", e - 1)?;
        }
        if e < interior {
            add_finite(&mut diag[e], k, "P1 diagonal assembly", e)?;
        }
        if e >= 1 && e < interior {
            add_finite(&mut off[e - 1], -k, "P1 off-diagonal assembly", e - 1)?;
        }
        // Load: ∫ f φ_a over the element, exact Gauss.
        for (gx, gw) in gauss5(x0, x1) {
            let fv = problem.forcing().eval(gx);
            let phi_left = (x1 - gx) / h;
            let phi_right = (gx - x0) / h;
            if !gx.is_finite()
                || !gw.is_finite()
                || !fv.is_finite()
                || !phi_left.is_finite()
                || !phi_right.is_finite()
            {
                return Err(non_finite("P1 quadrature", Some(e)));
            }
            if e >= 1 {
                add_finite(
                    &mut rhs[e - 1],
                    gw * fv * phi_left,
                    "P1 load assembly",
                    e - 1,
                )?;
            }
            if e < interior {
                add_finite(&mut rhs[e], gw * fv * phi_right, "P1 load assembly", e)?;
            }
        }
    }
    // Thomas solve.
    let mut c = try_zeroed("P1 Thomas upper factors", interior)?;
    let mut d = try_zeroed("P1 Thomas right-hand side", interior)?;
    let mut x = try_zeroed("P1 interior solution", interior)?;
    solve_tridiagonal_into(&off, &diag, &rhs, &mut c, &mut d, &mut x, "P1 solve")?;
    let mut full = try_zeroed("P1 full solution", n)?;
    full[1..=interior].copy_from_slice(&x);
    Ok(full)
}

/// Correctly-rounded 5-point Gauss–Legendre nodes/weights mapped to `[x0, x1]`.
/// The mathematical rule is exact for polynomial degree at most nine; this
/// ordinary-f64 helper is for solves, oracles, and Estimated diagnostics. The
/// verifier separately intervalizes the same constants and mapping.
#[must_use]
pub fn gauss5(x0: f64, x1: f64) -> [(f64, f64); 5] {
    const N: [f64; 5] = [
        -0.906_179_845_938_664,
        -0.538_469_310_105_683_1,
        0.0,
        0.538_469_310_105_683_1,
        0.906_179_845_938_664,
    ];
    const W: [f64; 5] = [
        0.236_926_885_056_189_08,
        0.478_628_670_499_366_47,
        0.568_888_888_888_888_9,
        0.478_628_670_499_366_47,
        0.236_926_885_056_189_08,
    ];
    let mid = f64::midpoint(x0, x1);
    let half = 0.5 * (x1 - x0);
    core::array::from_fn(|i| (mid + half * N[i], half * W[i]))
}

/// True energy-norm error `‖u′ − u_h′‖` (the ORACLE; high-resolution
/// f64 quadrature — the oracle needs accuracy, not rigor).
///
/// # Errors
/// Returns [`Fem1dError`] for malformed/nonconforming inputs, excessive work,
/// or any non-finite derived value. Oracle failure is never encoded as a
/// plausible scalar.
pub fn true_energy_error(problem: &MmsProblem, candidate: &[f64]) -> Result<f64, Fem1dError> {
    validate_problem(problem)?;
    validate_candidate(problem, candidate, "candidate")?;
    let m = problem.mesh();
    let du = problem.exact_solution().derive()?;
    validate_work(
        "energy-error oracle",
        &[m.len() - 1, 32, 5, du.0.len().saturating_add(8)],
    )?;
    let mut acc = 0.0;
    for e in 0..m.len() - 1 {
        let (x0, x1) = (m[e], m[e + 1]);
        let h = x1 - x0;
        let slope = (candidate[e + 1] - candidate[e]) / h;
        if !slope.is_finite() {
            return Err(non_finite("energy-error candidate slope", Some(e)));
        }
        // Subdivide each element for oracle accuracy.
        let sub = 32;
        for s in 0..sub {
            let a = x0 + h * f64::from(s) / f64::from(sub);
            let b = x0 + h * f64::from(s + 1) / f64::from(sub);
            for (gx, gw) in gauss5(a, b) {
                let d = du.eval(gx) - slope;
                let contribution = gw * d * d;
                let updated = acc + contribution;
                if !a.is_finite()
                    || !b.is_finite()
                    || !gx.is_finite()
                    || !gw.is_finite()
                    || !d.is_finite()
                    || !contribution.is_finite()
                    || !updated.is_finite()
                {
                    return Err(non_finite("energy-error quadrature", Some(e)));
                }
                acc = updated;
            }
        }
    }
    let error = acc.sqrt();
    if error.is_finite() {
        Ok(error)
    } else {
        Err(non_finite("energy-error result", None))
    }
}

/// One bounded Newton run, including an explicit convergence identity.
#[derive(Debug, Clone, PartialEq)]
pub struct NonlinearSolveReport {
    /// Final conforming nodal vector (also present on finite nonconvergence).
    pub solution: Vec<f64>,
    /// Newton updates performed.
    pub iterations: u32,
    /// Residual norm recomputed after the final update.
    pub residual_norm: f64,
    /// Whether `residual_norm < 1e-10`.
    pub converged: bool,
}

/// Convert finite nonconvergence into a structured caller refusal.
pub(crate) fn require_converged(
    report: &NonlinearSolveReport,
    stage: &'static str,
) -> Result<(), Fem1dError> {
    if report.converged {
        Ok(())
    } else {
        Err(Fem1dError::NonConverged {
            stage,
            iterations: report.iterations,
            residual_norm: report.residual_norm,
        })
    }
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn assemble_nonlinear(
    problem: &MmsProblem,
    u: &[f64],
    resid: &mut [f64],
    jac_diag: &mut [f64],
    jac_off: &mut [f64],
) -> Result<f64, Fem1dError> {
    let m = problem.mesh();
    let interior = m.len() - 2;
    if resid.len() != interior
        || jac_diag.len() != interior
        || jac_off.len() != interior.saturating_sub(1)
    {
        return Err(Fem1dError::LinearSystemShape {
            stage: "nonlinear assembly",
        });
    }
    resid.fill(0.0);
    jac_diag.fill(0.0);
    jac_off.fill(0.0);
    for e in 0..m.len() - 1 {
        let (x0, x1) = (m[e], m[e + 1]);
        let h = x1 - x0;
        let k = 1.0 / h;
        let slope_term = (u[e + 1] - u[e]) * k;
        if !slope_term.is_finite() {
            return Err(non_finite("nonlinear stiffness slope", Some(e)));
        }
        if e >= 1 {
            add_finite(
                &mut resid[e - 1],
                -slope_term,
                "nonlinear stiffness residual",
                e - 1,
            )?;
            add_finite(
                &mut jac_diag[e - 1],
                k,
                "nonlinear Jacobian diagonal",
                e - 1,
            )?;
        }
        if e < interior {
            add_finite(&mut resid[e], slope_term, "nonlinear stiffness residual", e)?;
            add_finite(&mut jac_diag[e], k, "nonlinear Jacobian diagonal", e)?;
        }
        if e >= 1 && e < interior {
            add_finite(
                &mut jac_off[e - 1],
                -k,
                "nonlinear Jacobian off-diagonal",
                e - 1,
            )?;
        }
        // Lumped nonlinear + load terms.
        let half_width = 0.5 * h;
        if e >= 1 {
            add_finite(
                &mut resid[e - 1],
                half_width * u[e].powi(3),
                "nonlinear cubic residual",
                e - 1,
            )?;
            add_finite(
                &mut jac_diag[e - 1],
                half_width * 3.0 * u[e] * u[e],
                "nonlinear cubic Jacobian",
                e - 1,
            )?;
        }
        if e < interior {
            add_finite(
                &mut resid[e],
                half_width * u[e + 1].powi(3),
                "nonlinear cubic residual",
                e,
            )?;
            add_finite(
                &mut jac_diag[e],
                half_width * 3.0 * u[e + 1] * u[e + 1],
                "nonlinear cubic Jacobian",
                e,
            )?;
        }
        for (gx, gw) in gauss5(x0, x1) {
            let forcing = problem.forcing().eval(gx) + problem.exact_solution().eval(gx).powi(3);
            let phi_left = (x1 - gx) / h;
            let phi_right = (gx - x0) / h;
            if !gx.is_finite()
                || !gw.is_finite()
                || !forcing.is_finite()
                || !phi_left.is_finite()
                || !phi_right.is_finite()
            {
                return Err(non_finite("nonlinear load quadrature", Some(e)));
            }
            if e >= 1 {
                add_finite(
                    &mut resid[e - 1],
                    -gw * forcing * phi_left,
                    "nonlinear load residual",
                    e - 1,
                )?;
            }
            if e < interior {
                add_finite(
                    &mut resid[e],
                    -gw * forcing * phi_right,
                    "nonlinear load residual",
                    e,
                )?;
            }
        }
    }
    let mut squared_norm = 0.0;
    for (index, value) in resid.iter().enumerate() {
        let updated = value.mul_add(*value, squared_norm);
        if !updated.is_finite() {
            return Err(non_finite("nonlinear residual norm", Some(index)));
        }
        squared_norm = updated;
    }
    let norm = squared_norm.sqrt();
    if norm.is_finite() {
        Ok(norm)
    } else {
        Err(non_finite("nonlinear residual norm", None))
    }
}

/// Newton solve of the toy nonlinear class `−u″ + u³ = f` from a conforming
/// `start`. Finite nonconvergence is explicit in the returned report.
///
/// # Errors
/// Returns [`Fem1dError`] for malformed inputs, excessive work, allocation
/// failure, non-finite arithmetic, or an unusable Thomas pivot.
#[allow(clippy::too_many_lines)]
pub fn solve_nonlinear(
    problem: &MmsProblem,
    start: &[f64],
    max_iter: u32,
) -> Result<NonlinearSolveReport, Fem1dError> {
    validate_problem(problem)?;
    validate_candidate(problem, start, "start")?;
    if max_iter > MAX_FEM1D_NEWTON_ITERATIONS {
        return Err(Fem1dError::ResourceLimit {
            resource: "Newton iterations",
            requested: max_iter as usize,
            limit: MAX_FEM1D_NEWTON_ITERATIONS as usize,
        });
    }
    let passes = (max_iter as usize)
        .checked_add(1)
        .ok_or(Fem1dError::WorkBudgetExceeded {
            stage: "nonlinear solve",
            estimated: None,
            limit: MAX_FEM1D_WORK_UNITS,
        })?;
    let per_cell = problem
        .exact_solution()
        .coefficients()
        .len()
        .checked_add(problem.forcing().coefficients().len())
        .and_then(|polynomial_work| polynomial_work.checked_mul(5))
        .and_then(|polynomial_work| polynomial_work.checked_add(32))
        .ok_or(Fem1dError::WorkBudgetExceeded {
            stage: "nonlinear solve",
            estimated: None,
            limit: MAX_FEM1D_WORK_UNITS,
        })?;
    validate_work(
        "nonlinear solve",
        &[problem.mesh().len() - 1, passes, per_cell],
    )?;

    let interior = problem.mesh().len() - 2;
    let mut u = try_copy_slice("nonlinear start", start)?;
    let mut resid = try_zeroed("nonlinear residual", interior)?;
    let mut jac_diag = try_zeroed("nonlinear Jacobian diagonal", interior)?;
    let mut jac_off = try_zeroed(
        "nonlinear Jacobian off-diagonal",
        interior.saturating_sub(1),
    )?;
    let mut rhs = try_zeroed("nonlinear right-hand side", interior)?;
    let mut c = try_zeroed("nonlinear Thomas upper factors", interior)?;
    let mut d = try_zeroed("nonlinear Thomas right-hand side", interior)?;
    let mut delta = try_zeroed("nonlinear Newton update", interior)?;

    let mut residual_norm =
        assemble_nonlinear(problem, &u, &mut resid, &mut jac_diag, &mut jac_off)?;
    if residual_norm < NEWTON_RESIDUAL_TOLERANCE {
        return Ok(NonlinearSolveReport {
            solution: u,
            iterations: 0,
            residual_norm,
            converged: true,
        });
    }
    for iteration in 1..=max_iter {
        for (right, residual) in rhs.iter_mut().zip(&resid) {
            *right = -*residual;
        }
        solve_tridiagonal_into(
            &jac_off,
            &jac_diag,
            &rhs,
            &mut c,
            &mut d,
            &mut delta,
            "nonlinear Newton solve",
        )?;
        for (index, update) in delta.iter().copied().enumerate() {
            let updated = u[index + 1] + update;
            if !updated.is_finite() {
                return Err(non_finite("nonlinear solution update", Some(index + 1)));
            }
            u[index + 1] = updated;
        }
        residual_norm = assemble_nonlinear(problem, &u, &mut resid, &mut jac_diag, &mut jac_off)?;
        if residual_norm < NEWTON_RESIDUAL_TOLERANCE {
            return Ok(NonlinearSolveReport {
                solution: u,
                iterations: iteration,
                residual_norm,
                converged: true,
            });
        }
    }
    Ok(NonlinearSolveReport {
        solution: u,
        iterations: max_iter,
        residual_norm,
        converged: false,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn poly(coefficients: Vec<f64>) -> Poly {
        Poly::new(coefficients).expect("valid test polynomial")
    }

    fn test_problem(mesh: Vec<f64>) -> MmsProblem {
        MmsProblem::new("fem1d-test", poly(vec![0.0, 1.0, -1.0]), mesh).expect("valid test problem")
    }

    #[test]
    fn malformed_public_inputs_refuse_before_indexing() {
        for mesh in [Vec::new(), vec![0.0]] {
            let error = MmsProblem::new("fem1d-test", poly(vec![0.0, 1.0, -1.0]), mesh)
                .expect_err("short mesh must refuse during construction");
            assert!(matches!(
                error,
                Fem1dError::ResourceLimit {
                    resource: "mesh nodes",
                    ..
                }
            ));
        }

        let problem = test_problem(vec![0.0, 0.5, 1.0]);
        for values in [&[0.0, 0.0][..], &[0.0, 0.0, 0.0, 0.0][..]] {
            assert!(matches!(
                solve_nonlinear(&problem, values, 1),
                Err(Fem1dError::CandidateLength { field: "start", .. })
            ));
            assert!(matches!(
                true_energy_error(&problem, values),
                Err(Fem1dError::CandidateLength {
                    field: "candidate",
                    ..
                })
            ));
        }
        assert!(matches!(
            solve_nonlinear(&problem, &[1.0, 0.0, 0.0], 1),
            Err(Fem1dError::CandidateBoundary { field: "start" })
        ));
        assert!(matches!(
            solve_nonlinear(&problem, &[0.0, f64::NAN, 0.0], 1),
            Err(Fem1dError::NonFiniteCandidate { field: "start", .. })
        ));

        let duplicate = MmsProblem::new(
            "duplicate",
            poly(vec![0.0, 1.0, -1.0]),
            vec![0.0, 0.5, 0.5, 1.0],
        );
        assert!(matches!(
            duplicate,
            Err(Fem1dError::NonIncreasingMeshCell { .. })
        ));
    }

    #[test]
    fn nonlinear_convergence_identity_is_explicit() {
        let problem = test_problem(vec![0.0, 0.5, 1.0]);
        let zero = [0.0, 0.0, 0.0];
        let zero_budget = solve_nonlinear(&problem, &zero, 0).expect("finite initial residual");
        assert!(!zero_budget.converged);
        assert_eq!(zero_budget.iterations, 0);
        assert!(zero_budget.residual_norm >= NEWTON_RESIDUAL_TOLERANCE);
        assert!(matches!(
            require_converged(&zero_budget, "test solve"),
            Err(Fem1dError::NonConverged {
                stage: "test solve",
                iterations: 0,
                ..
            })
        ));

        let solved = solve_nonlinear(&problem, &zero, 50).expect("bounded Newton solve");
        assert!(solved.converged);
        assert!(solved.residual_norm < NEWTON_RESIDUAL_TOLERANCE);
        assert_eq!(solved.solution.len(), problem.mesh().len());
        assert_eq!(solved.solution[0].to_bits(), 0.0_f64.to_bits());
        assert_eq!(
            solved.solution[solved.solution.len() - 1].to_bits(),
            0.0_f64.to_bits()
        );

        assert!(matches!(
            solve_nonlinear(&problem, &zero, MAX_FEM1D_NEWTON_ITERATIONS + 1),
            Err(Fem1dError::ResourceLimit {
                resource: "Newton iterations",
                ..
            })
        ));
    }

    #[test]
    fn gauss_constants_cover_independent_truth_brackets() {
        let rule = gauss5(-1.0, 1.0);
        let positive_constants = [
            (
                rule[4].0,
                0x3fec_ff6c_e053_3a69,
                0x3fec_ff6c_e053_3a69,
                0x3fec_ff6c_e053_3a6a,
            ),
            (
                rule[3].0,
                0x3fe1_3b23_fd99_b705,
                0x3fe1_3b23_fd99_b704,
                0x3fe1_3b23_fd99_b705,
            ),
            (
                rule[0].1,
                0x3fce_539e_c36e_038c,
                0x3fce_539e_c36e_038c,
                0x3fce_539e_c36e_038d,
            ),
            (
                rule[1].1,
                0x3fde_a1da_25ae_415b,
                0x3fde_a1da_25ae_415a,
                0x3fde_a1da_25ae_415b,
            ),
            (
                rule[2].1,
                0x3fe2_3456_789a_bcdf,
                0x3fe2_3456_789a_bcdf,
                0x3fe2_3456_789a_bce0,
            ),
        ];
        for (rounded, expected_bits, lower_bits, upper_bits) in positive_constants {
            let lower = f64::from_bits(lower_bits);
            let upper = f64::from_bits(upper_bits);
            assert_eq!(rounded.to_bits(), expected_bits);
            assert!(lower <= rounded && rounded <= upper);
        }
    }

    #[test]
    fn exact_boundary_sum_distinguishes_hidden_residue_from_true_cancellation() {
        let hidden_residue = poly(vec![0.0, 1.0e16, -1.0e16, 1.0]);
        assert_eq!(hidden_residue.eval(1.0).to_bits(), 0.0_f64.to_bits());
        assert!(
            !hidden_residue.is_exactly_zero_at_one(),
            "an absorbed exact residue is not a homogeneous trace"
        );

        let rounded_point = poly(vec![0.0, 1.0, 1.0e16, -1.0e16, -1.0]);
        assert_ne!(rounded_point.eval(1.0).to_bits(), 0.0_f64.to_bits());
        assert!(
            rounded_point.is_exactly_zero_at_one(),
            "the exact binary-rational coefficient sum is authoritative"
        );

        let extremes = poly(vec![
            0.0,
            f64::MAX,
            -f64::MAX,
            f64::MIN_POSITIVE,
            -f64::MIN_POSITIVE,
        ]);
        assert!(extremes.is_exactly_zero_at_one());
        assert!(matches!(
            Poly::new(vec![0.0, f64::INFINITY, f64::NEG_INFINITY]),
            Err(Fem1dError::NonFinitePolynomialCoefficient { .. })
        ));
    }

    #[test]
    fn identity_builder_failures_map_to_owner_errors() {
        assert_eq!(
            identity_build_error(
                IdentityBuildError::AllocationFailed { requested: 17 },
                "identity stage",
                "identity resource",
                64,
            ),
            Fem1dError::AllocationFailed {
                stage: "identity stage",
                requested: 17,
            }
        );
        assert_eq!(
            identity_build_error(
                IdentityBuildError::CanonicalBytesExceeded {
                    requested: 65,
                    limit: 64,
                },
                "identity stage",
                "identity resource",
                64,
            ),
            Fem1dError::ResourceLimit {
                resource: "identity resource",
                requested: 65,
                limit: 64,
            }
        );
        assert_eq!(
            identity_build_error(
                IdentityBuildError::FramedLengthNotRepresentable { length: 65 },
                "identity stage",
                "identity resource",
                64,
            ),
            Fem1dError::ResourceLimit {
                resource: "identity resource",
                requested: 65,
                limit: 64,
            }
        );
        assert_eq!(
            identity_build_error(
                IdentityBuildError::CanonicalLengthOverflow,
                "identity stage",
                "identity resource",
                64,
            ),
            Fem1dError::ResourceLimit {
                resource: "identity resource",
                requested: usize::MAX,
                limit: 64,
            }
        );
    }

    #[test]
    fn maximum_mms_identity_reaches_the_derived_canonical_caps() {
        let name = "m".repeat(MAX_FEM1D_CLASS_NAME_BYTES);
        let class = MmsClass::new(&name, poly(vec![0.0, 1.0, 1.0, 1.0, 1.0, -4.0]))
            .expect("maximum-size degree-five class is admitted");
        assert_eq!(
            class.canonical_bytes().len(),
            MAX_FEM1D_CLASS_CANONICAL_IDENTITY_BYTES
        );

        let denominator = (MAX_FEM1D_MESH_NODES - 1) as f64;
        let mesh: Vec<f64> = (0..MAX_FEM1D_MESH_NODES)
            .map(|index| index as f64 / denominator)
            .collect();
        let problem = MmsProblem::from_class(class, mesh)
            .expect("maximum-size mesh identity reaches the exact admitted cap");
        assert_eq!(
            problem.canonical_bytes().len(),
            MAX_FEM1D_PROBLEM_CANONICAL_IDENTITY_BYTES
        );
    }

    #[test]
    #[allow(clippy::too_many_lines)] // one identity-field mutation matrix
    fn mms_class_identity_fields_move_independently() {
        let exact = poly(vec![0.0, 1.0, -1.0]);
        let forcing = exact.derive().unwrap().derive().unwrap().neg();
        let antiderivative = forcing.antiderive().unwrap();
        let exact_bytes = f64_bytes(&exact.0, "class identity test").unwrap();
        let forcing_bytes = f64_bytes(&forcing.0, "class identity test").unwrap();
        let antiderivative_bytes = f64_bytes(&antiderivative.0, "class identity test").unwrap();
        let cap = MAX_FEM1D_CLASS_CANONICAL_IDENTITY_BYTES + 64;
        let root = |domain: &str,
                    version: u32,
                    name: &str,
                    exact: &[u8],
                    forcing: &[u8],
                    antiderivative: &[u8]| {
            class_identity_from_parts_with_budget(
                domain,
                version,
                name,
                exact,
                forcing,
                antiderivative,
                cap,
            )
            .unwrap()
            .root()
        };
        let baseline = root(
            MMS_CLASS_IDENTITY_DOMAIN,
            MMS_CLASS_IDENTITY_VERSION,
            "elliptic",
            &exact_bytes,
            &forcing_bytes,
            &antiderivative_bytes,
        );
        assert_ne!(
            baseline,
            root(
                "fs-verify/fem1d-mms-class-shadow.v2",
                MMS_CLASS_IDENTITY_VERSION,
                "elliptic",
                &exact_bytes,
                &forcing_bytes,
                &antiderivative_bytes,
            )
        );
        assert_ne!(
            baseline,
            root(
                MMS_CLASS_IDENTITY_DOMAIN,
                MMS_CLASS_IDENTITY_VERSION + 1,
                "elliptic",
                &exact_bytes,
                &forcing_bytes,
                &antiderivative_bytes,
            )
        );
        assert_ne!(
            baseline,
            root(
                MMS_CLASS_IDENTITY_DOMAIN,
                MMS_CLASS_IDENTITY_VERSION,
                "elliptic-renamed",
                &exact_bytes,
                &forcing_bytes,
                &antiderivative_bytes,
            )
        );

        let value = exact.0[1];
        let adjacent = f64::from_bits(value.to_bits() + 1);
        assert_eq!(format!("{value:.6e}"), format!("{adjacent:.6e}"));
        assert_ne!(value.to_bits(), adjacent.to_bits());
        let mut changed_exact = exact_bytes.clone();
        changed_exact[8..16].copy_from_slice(&adjacent.to_bits().to_le_bytes());
        assert_ne!(
            baseline,
            root(
                MMS_CLASS_IDENTITY_DOMAIN,
                MMS_CLASS_IDENTITY_VERSION,
                "elliptic",
                &changed_exact,
                &forcing_bytes,
                &antiderivative_bytes,
            )
        );
        let mut changed_forcing = forcing_bytes.clone();
        changed_forcing[0] ^= 1;
        assert_ne!(
            baseline,
            root(
                MMS_CLASS_IDENTITY_DOMAIN,
                MMS_CLASS_IDENTITY_VERSION,
                "elliptic",
                &exact_bytes,
                &changed_forcing,
                &antiderivative_bytes,
            )
        );
        let mut changed_antiderivative = antiderivative_bytes.clone();
        changed_antiderivative[0] ^= 1;
        assert_ne!(
            baseline,
            root(
                MMS_CLASS_IDENTITY_DOMAIN,
                MMS_CLASS_IDENTITY_VERSION,
                "elliptic",
                &exact_bytes,
                &forcing_bytes,
                &changed_antiderivative,
            )
        );
    }

    #[test]
    fn mms_problem_identity_fields_move_independently() {
        let class = MmsClass::new("elliptic", poly(vec![0.0, 1.0, -1.0])).unwrap();
        let alternate = MmsClass::new("elliptic-alt", poly(vec![0.0, 1.0, -1.0])).unwrap();
        let mesh_bytes = f64_bytes(&[0.0, 0.5, 1.0], "problem identity test").unwrap();
        let cap = MAX_FEM1D_PROBLEM_CANONICAL_IDENTITY_BYTES;
        let root = |domain: &str,
                    version: u32,
                    child: &ReplayIdentity,
                    class_bytes: &[u8],
                    mesh: &[u8]| {
            problem_identity_from_parts_with_budget(domain, version, child, class_bytes, mesh, cap)
                .unwrap()
                .root()
        };
        let baseline = root(
            MMS_PROBLEM_IDENTITY_DOMAIN,
            MMS_PROBLEM_IDENTITY_VERSION,
            class.identity(),
            class.canonical_bytes(),
            &mesh_bytes,
        );
        assert_ne!(
            baseline,
            root(
                "fs-verify/fem1d-mms-problem-shadow.v2",
                MMS_PROBLEM_IDENTITY_VERSION,
                class.identity(),
                class.canonical_bytes(),
                &mesh_bytes,
            )
        );
        assert_ne!(
            baseline,
            root(
                MMS_PROBLEM_IDENTITY_DOMAIN,
                MMS_PROBLEM_IDENTITY_VERSION + 1,
                class.identity(),
                class.canonical_bytes(),
                &mesh_bytes,
            )
        );
        assert_ne!(
            baseline,
            root(
                MMS_PROBLEM_IDENTITY_DOMAIN,
                MMS_PROBLEM_IDENTITY_VERSION,
                alternate.identity(),
                class.canonical_bytes(),
                &mesh_bytes,
            )
        );
        let mut changed_class_bytes = class.canonical_bytes().to_vec();
        changed_class_bytes[0] ^= 1;
        assert_ne!(
            baseline,
            root(
                MMS_PROBLEM_IDENTITY_DOMAIN,
                MMS_PROBLEM_IDENTITY_VERSION,
                class.identity(),
                &changed_class_bytes,
                &mesh_bytes,
            )
        );
        let mut changed_mesh = mesh_bytes.clone();
        changed_mesh[8] ^= 1;
        assert_ne!(
            baseline,
            root(
                MMS_PROBLEM_IDENTITY_DOMAIN,
                MMS_PROBLEM_IDENTITY_VERSION,
                class.identity(),
                class.canonical_bytes(),
                &changed_mesh,
            )
        );
    }

    #[test]
    fn mms_class_identity_receipts_fail_closed() {
        let class = MmsClass::new("elliptic", poly(vec![0.0, 1.0, -1.0])).unwrap();
        let current = class.identity_receipt().unwrap();
        assert_eq!(class.admit_identity_receipt(&current), Ok(()));
        for declared in [
            MMS_CLASS_IDENTITY_VERSION - 1,
            MMS_CLASS_IDENTITY_VERSION + 1,
        ] {
            assert_eq!(
                MmsClassIdentityReceipt::from_retained_parts(
                    declared,
                    current.canonical_bytes().to_vec(),
                    current.root(),
                ),
                Err(MmsIdentityAdmissionError::UnsupportedVersion {
                    identity: "MMS class",
                    declared,
                    supported: MMS_CLASS_IDENTITY_VERSION,
                })
            );
        }
        assert!(matches!(
            MmsClassIdentityReceipt::from_retained_parts(
                MMS_CLASS_IDENTITY_VERSION,
                current.canonical_bytes().to_vec(),
                current.root() ^ 1,
            ),
            Err(MmsIdentityAdmissionError::RootMismatch {
                identity: "MMS class"
            })
        ));
        let mut changed_bytes = current.canonical_bytes().to_vec();
        changed_bytes[0] ^= 1;
        assert!(matches!(
            MmsClassIdentityReceipt::from_retained_parts(
                MMS_CLASS_IDENTITY_VERSION,
                changed_bytes,
                current.root(),
            ),
            Err(MmsIdentityAdmissionError::RootMismatch {
                identity: "MMS class"
            })
        ));
        assert!(matches!(
            MmsClassIdentityReceipt::from_retained_parts(
                MMS_CLASS_IDENTITY_VERSION,
                vec![0; MAX_FEM1D_CLASS_CANONICAL_IDENTITY_BYTES + 1],
                0,
            ),
            Err(MmsIdentityAdmissionError::CanonicalBytesExceeded {
                identity: "MMS class",
                requested,
                limit: MAX_FEM1D_CLASS_CANONICAL_IDENTITY_BYTES,
            }) if requested == MAX_FEM1D_CLASS_CANONICAL_IDENTITY_BYTES + 1
        ));
        let foreign = MmsClass::new("foreign", poly(vec![0.0, 1.0, -1.0])).unwrap();
        assert!(matches!(
            class.admit_identity_receipt(&foreign.identity_receipt().unwrap()),
            Err(MmsIdentityAdmissionError::IdentityMismatch {
                identity: "MMS class"
            })
        ));
    }

    #[test]
    fn mms_problem_identity_receipts_fail_closed() {
        let class = MmsClass::new("elliptic", poly(vec![0.0, 1.0, -1.0])).unwrap();
        let problem = MmsProblem::from_class(class, vec![0.0, 0.5, 1.0]).unwrap();
        let current = problem.identity_receipt().unwrap();
        assert_eq!(problem.admit_identity_receipt(&current), Ok(()));
        for declared in [
            MMS_PROBLEM_IDENTITY_VERSION - 1,
            MMS_PROBLEM_IDENTITY_VERSION + 1,
        ] {
            assert_eq!(
                MmsProblemIdentityReceipt::from_retained_parts(
                    declared,
                    current.canonical_bytes().to_vec(),
                    current.root(),
                ),
                Err(MmsIdentityAdmissionError::UnsupportedVersion {
                    identity: "MMS problem",
                    declared,
                    supported: MMS_PROBLEM_IDENTITY_VERSION,
                })
            );
        }
        assert!(matches!(
            MmsProblemIdentityReceipt::from_retained_parts(
                MMS_PROBLEM_IDENTITY_VERSION,
                current.canonical_bytes().to_vec(),
                current.root() ^ 1,
            ),
            Err(MmsIdentityAdmissionError::RootMismatch {
                identity: "MMS problem"
            })
        ));
        let mut changed_bytes = current.canonical_bytes().to_vec();
        changed_bytes[0] ^= 1;
        assert!(matches!(
            MmsProblemIdentityReceipt::from_retained_parts(
                MMS_PROBLEM_IDENTITY_VERSION,
                changed_bytes,
                current.root(),
            ),
            Err(MmsIdentityAdmissionError::RootMismatch {
                identity: "MMS problem"
            })
        ));
        assert!(matches!(
            MmsProblemIdentityReceipt::from_retained_parts(
                MMS_PROBLEM_IDENTITY_VERSION,
                vec![0; MAX_FEM1D_PROBLEM_CANONICAL_IDENTITY_BYTES + 1],
                0,
            ),
            Err(MmsIdentityAdmissionError::CanonicalBytesExceeded {
                identity: "MMS problem",
                requested,
                limit: MAX_FEM1D_PROBLEM_CANONICAL_IDENTITY_BYTES,
            }) if requested == MAX_FEM1D_PROBLEM_CANONICAL_IDENTITY_BYTES + 1
        ));
        let foreign =
            MmsProblem::new("elliptic", poly(vec![0.0, 1.0, -1.0]), vec![0.0, 0.25, 1.0]).unwrap();
        assert!(matches!(
            problem.admit_identity_receipt(&foreign.identity_receipt().unwrap()),
            Err(MmsIdentityAdmissionError::IdentityMismatch {
                identity: "MMS problem"
            })
        ));
    }

    #[test]
    #[allow(clippy::too_many_lines)] // one canonical identity and semantic-field matrix
    fn canonical_class_and_problem_identity_bind_every_semantic_field() {
        assert_eq!(MAX_FEM1D_CLASS_CANONICAL_IDENTITY_BYTES, 4_441);
        assert_eq!(MAX_FEM1D_PROBLEM_CANONICAL_IDENTITY_BYTES, 8_004_626);
        let normalized = MmsClass::new(
            "elliptic",
            poly(vec![-0.0, 1.0, -1.0, -0.0, 0.0, -0.0, 0.0]),
        )
        .expect("canonical class");
        let ordinary =
            MmsClass::new("elliptic", poly(vec![0.0, 1.0, -1.0])).expect("same canonical class");
        assert_eq!(normalized, ordinary);
        assert_eq!(normalized.exact_solution().coefficients().len(), 3);
        assert_eq!(
            normalized.exact_solution().coefficients()[0].to_bits(),
            0.0_f64.to_bits()
        );
        assert_eq!(normalized.forcing().coefficients(), &[2.0]);
        assert_eq!(
            normalized.rounded_forcing_antiderivative().coefficients(),
            &[0.0, 2.0]
        );
        assert_eq!(normalized.identity(), ordinary.identity());
        assert_eq!(normalized.canonical_bytes(), ordinary.canonical_bytes());
        assert_eq!(normalized.identity().root(), 0x959a_7771_9f30_8c27);
        assert_eq!(normalized.canonical_bytes().len(), 281);
        let exact_class = class_identity_with_budget(
            normalized.name(),
            normalized.exact_solution(),
            normalized.forcing(),
            normalized.rounded_forcing_antiderivative(),
            normalized.canonical_bytes().len(),
        )
        .expect("the exact class identity cap is admitted");
        assert_eq!(&exact_class, normalized.identity());
        assert!(matches!(
            class_identity_with_budget(
                normalized.name(),
                normalized.exact_solution(),
                normalized.forcing(),
                normalized.rounded_forcing_antiderivative(),
                normalized.canonical_bytes().len() - 1,
            ),
            Err(Fem1dError::ResourceLimit {
                resource: "MMS class canonical identity bytes",
                requested: 281,
                limit: 280,
            })
        ));

        let renamed =
            MmsClass::new("elliptic-renamed", poly(vec![0.0, 1.0, -1.0])).expect("renamed class");
        let rescaled =
            MmsClass::new("elliptic", poly(vec![0.0, 2.0, -2.0])).expect("rescaled class");
        assert_ne!(normalized.identity(), renamed.identity());
        assert_ne!(normalized.identity(), rescaled.identity());
        assert_ne!(normalized.forcing(), rescaled.forcing());

        let coarse = MmsProblem::from_class(normalized.clone(), vec![-0.0, 0.5, 1.0])
            .expect("canonicalized mesh");
        let same =
            MmsProblem::from_class(ordinary, vec![0.0, 0.5, 1.0]).expect("same canonical mesh");
        let shifted = MmsProblem::from_class(normalized.clone(), vec![0.0, 0.75, 1.0])
            .expect("same-size different mesh");
        let refined =
            MmsProblem::from_class(normalized, vec![0.0, 0.25, 0.5, 1.0]).expect("different mesh");
        let reused = coarse
            .with_mesh(vec![0.0, 0.25, 0.5, 1.0])
            .expect("fallible class reuse on a new mesh");
        assert_eq!(coarse.identity(), same.identity());
        assert_eq!(coarse.canonical_bytes(), same.canonical_bytes());
        assert_eq!(coarse.identity().root(), 0x7148_ea04_d660_5664);
        assert_eq!(coarse.canonical_bytes().len(), 490);
        let exact_problem = problem_identity_with_budget(
            coarse.class(),
            coarse.mesh(),
            coarse.canonical_bytes().len(),
        )
        .expect("the exact problem identity cap is admitted");
        assert_eq!(&exact_problem, coarse.identity());
        assert!(matches!(
            problem_identity_with_budget(
                coarse.class(),
                coarse.mesh(),
                coarse.canonical_bytes().len() - 1,
            ),
            Err(Fem1dError::ResourceLimit {
                resource: "MMS problem canonical identity bytes",
                requested: 490,
                limit: 489,
            })
        ));
        assert_ne!(coarse.identity(), shifted.identity());
        assert_ne!(coarse.identity(), refined.identity());
        assert_eq!(reused.class().identity(), coarse.class().identity());
        assert_eq!(
            reused.class().canonical_bytes(),
            coarse.class().canonical_bytes()
        );
        assert_eq!(reused.identity(), refined.identity());
        assert_eq!(reused.canonical_bytes(), refined.canonical_bytes());

        let rounded = MmsClass::new(
            "rounded-antiderivative",
            poly(vec![0.0, 0.1, 0.0, 0.0, -0.1]),
        )
        .expect("nontrivial rounded antiderivative class");
        let rounded_coefficient = rounded.rounded_forcing_antiderivative().coefficients()[3];
        assert_ne!(rounded_coefficient.to_bits(), 0.4_f64.to_bits());
        assert_eq!(
            rounded_coefficient.to_bits(),
            (rounded.forcing().coefficients()[2] / 3.0).to_bits()
        );

        assert!(matches!(
            Poly::new(Vec::new()),
            Err(Fem1dError::PolynomialCoefficientCount { count: 0, .. })
        ));
        assert!(matches!(
            Poly::new(vec![0.0; MAX_FEM1D_RAW_POLY_COEFFICIENTS + 1]),
            Err(Fem1dError::ResourceLimit {
                resource: "raw polynomial coefficients",
                ..
            })
        ));
        assert!(matches!(
            MmsClass::new("not-homogeneous", poly(vec![0.0, 1.0])),
            Err(Fem1dError::ExactSolutionBoundary)
        ));
        assert!(matches!(
            MmsProblem::new("", poly(vec![0.0]), Vec::new()),
            Err(Fem1dError::ResourceLimit {
                resource: "mesh nodes",
                requested: 0,
                ..
            })
        ));
    }
}
