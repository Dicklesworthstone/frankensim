//! Admitted finite derived and stratified machine-geometry objects (RD.1a).
//!
//! This moonshot-gated module defines the *objects* consumed by later theorem
//! lanes.  Admission is deliberately structural: it records the algebraic
//! category, coefficient semantics, local chart, equality/inequality/contact
//! roles, finite complexes, stratification, units, frames, computability, and
//! proof-state metadata without claiming that a retained theorem is true.
//!
//! Equality, inequality, contact, and constitutive identities are nominally
//! distinct.  Digest equality cannot erase that distinction:
//!
//! ```compile_fail
//! use fs_geom::derived::{EqualityConstraintIdV1, InequalityConstraintIdV1};
//!
//! fn needs_equality(_: EqualityConstraintIdV1) {}
//! let inequality = InequalityConstraintIdV1::from_bytes([7; 32]);
//! needs_equality(inequality);
//! ```

#![allow(clippy::too_many_lines)]

use core::fmt;

use fs_blake3::identity::{
    CanonicalEncoder, CanonicalError, CanonicalLimits, CanonicalSchema, Field, FieldSpec,
    IdentityReceipt, ProblemSemanticId, WireType,
};
use fs_exec::Cx;

/// Current wire and semantic version for RD.1a object admission.
pub const DERIVED_GEOMETRY_SCHEMA_VERSION_V1: u32 = 1;
/// Absolute number of objects admitted in any one top-level collection.
pub const DERIVED_GEOMETRY_HARD_MAX_OBJECTS_V1: usize = 4096;
/// Absolute sum of all finite graded-space dimensions.
pub const DERIVED_GEOMETRY_HARD_MAX_TOTAL_RANK_V1: u64 = 1 << 24;
/// Absolute finite coordinate/ambient dimension.
pub const DERIVED_GEOMETRY_HARD_MAX_DIMENSION_V1: u32 = 1 << 20;
/// Absolute canonical-frame byte ceiling.
pub const DERIVED_GEOMETRY_HARD_MAX_CANONICAL_BYTES_V1: u64 = 1 << 24;
/// Absolute single-field byte ceiling.
pub const DERIVED_GEOMETRY_HARD_MAX_FIELD_BYTES_V1: u64 = 1 << 23;
/// Maximum number of nested rank entries processed between cancellation polls.
const DERIVED_GEOMETRY_CANCELLATION_STRIDE_V1: usize = 64;

trait DigestBytes {
    fn digest_bytes(&self) -> &[u8; 32];
    fn is_zero(&self) -> bool {
        self.digest_bytes().iter().all(|byte| *byte == 0)
    }
}

macro_rules! opaque_id {
    ($(#[$meta:meta])* $name:ident) => {
        $(#[$meta])*
        #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name([u8; 32]);

        impl $name {
            /// Construct an explicitly typed identity from exact digest bytes.
            /// The bytes identify content; they do not prove scientific truth.
            #[must_use]
            pub const fn from_bytes(bytes: [u8; 32]) -> Self {
                Self(bytes)
            }

            /// Exact digest bytes.
            #[must_use]
            pub const fn as_bytes(&self) -> &[u8; 32] {
                &self.0
            }
        }

        impl DigestBytes for $name {
            fn digest_bytes(&self) -> &[u8; 32] {
                self.as_bytes()
            }
        }
    };
}

opaque_id!(
    /// Physical subject whose machine geometry is represented.
    DerivedSubjectIdV1
);
opaque_id!(
    /// Immutable source-model version.
    DerivedModelVersionIdV1
);
opaque_id!(
    /// One finite configuration chart.
    ConfigurationChartIdV1
);
opaque_id!(
    /// Coordinate-frame convention.
    DerivedFrameIdV1
);
opaque_id!(
    /// Unit-system convention.
    DerivedUnitSystemIdV1
);
opaque_id!(
    /// Semantic physical quantity kind.
    DerivedQuantityKindIdV1
);
opaque_id!(
    /// Exact polynomial artifact.
    PolynomialIdV1
);
opaque_id!(
    /// Admitted restricted-analytic program.
    AnalyticProgramIdV1
);
opaque_id!(
    /// Retained mathematical or admission witness.
    DerivedWitnessIdV1
);
opaque_id!(
    /// Explicit no-claim or non-applicability artifact.
    DerivedNoClaimIdV1
);
opaque_id!(
    /// Equality-constraint germ identity.
    EqualityConstraintIdV1
);
opaque_id!(
    /// Inequality-constraint germ identity.
    InequalityConstraintIdV1
);
opaque_id!(
    /// Relative-boundary identity.
    RelativeBoundaryIdV1
);
opaque_id!(
    /// Unilateral-contact identity.
    ContactConstraintIdV1
);
opaque_id!(
    /// Constitutive relation identity.
    ConstitutiveDatumIdV1
);
opaque_id!(
    /// Finite tangent, cotangent, or deformation complex.
    DerivedComplexIdV1
);
opaque_id!(
    /// Linear differential/map artifact.
    DerivedLinearMapIdV1
);
opaque_id!(
    /// Finite resolution/truncation artifact.
    DerivedResolutionIdV1
);
opaque_id!(
    /// Derived local-model identity.
    DerivedLocalModelIdV1
);
opaque_id!(
    /// Finite stratification identity.
    StratificationIdV1
);
opaque_id!(
    /// One stratum identity.
    StratumIdV1
);
opaque_id!(
    /// One compact local-link identity.
    LocalLinkIdV1
);
opaque_id!(
    /// External checker implementation identity.
    DerivedCheckerIdV1
);
opaque_id!(
    /// Retained theorem-card identity.
    DerivedTheoremIdV1
);

/// Domain-separated semantic identity schema for one complete admitted object.
pub enum DerivedGeometryIdentitySchemaV1 {}

impl CanonicalSchema for DerivedGeometryIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-geom.derived-machine-geometry.v1";
    const NAME: &'static str = "derived-stratified-machine-geometry";
    const VERSION: u32 = DERIVED_GEOMETRY_SCHEMA_VERSION_V1;
    const CONTEXT: &'static str = "finite admitted category, charts, typed constraints, complexes, local models, stratification, and proof metadata";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("subject", WireType::Bytes),
        FieldSpec::required("model-version", WireType::Bytes),
        FieldSpec::required("category", WireType::Bytes),
        FieldSpec::required("global-context", WireType::Bytes),
        FieldSpec::required("charts", WireType::CanonicalSet),
        FieldSpec::required("equalities", WireType::CanonicalSet),
        FieldSpec::required("inequalities", WireType::CanonicalSet),
        FieldSpec::required("boundaries", WireType::CanonicalSet),
        FieldSpec::required("contacts", WireType::CanonicalSet),
        FieldSpec::required("constitutive-data", WireType::CanonicalSet),
        FieldSpec::required("complexes", WireType::CanonicalSet),
        FieldSpec::required("local-models", WireType::CanonicalSet),
        FieldSpec::required("stratification", WireType::Bytes),
        FieldSpec::required("proof-state", WireType::Bytes),
    ];
}

/// Typed semantic identity of a structurally admitted RD.1a object.
pub type DerivedGeometryIdV1 = ProblemSemanticId<DerivedGeometryIdentitySchemaV1>;

/// Declared coefficient field/ring and real-versus-complex semantics.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CoefficientSystemV1 {
    /// Exact rational coefficients with ordered real semantics.
    RationalReal,
    /// Exact real-algebraic coefficients.
    AlgebraicReal,
    /// Outward dyadic interval coefficients at a fixed nonzero precision.
    DyadicIntervalReal {
        /// Stored mantissa precision.
        precision_bits: u16,
    },
    /// Exact rational coefficients over the complex numbers.
    RationalComplex,
    /// Exact complex-algebraic coefficients.
    AlgebraicComplex,
}

impl CoefficientSystemV1 {
    const fn tag(self) -> u8 {
        match self {
            Self::RationalReal => 0,
            Self::AlgebraicReal => 1,
            Self::DyadicIntervalReal { .. } => 2,
            Self::RationalComplex => 3,
            Self::AlgebraicComplex => 4,
        }
    }

    /// Whether sign/order, inequality, and normal-cone semantics are defined.
    #[must_use]
    pub const fn is_ordered_real(self) -> bool {
        matches!(
            self,
            Self::RationalReal | Self::AlgebraicReal | Self::DyadicIntervalReal { .. }
        )
    }
}

/// Mathematical category in which local objects and encodings live.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GeometricCategoryV1 {
    /// Finite semialgebraic real geometry.
    Semialgebraic,
    /// Finite algebraic geometry over the declared real or complex field.
    Algebraic,
    /// Restricted real-analytic local geometry over an admitted primitive set.
    RestrictedAnalytic,
    /// Finite subanalytic presentation with a retained construction witness.
    Subanalytic {
        /// Witness for the admitted projection/closure construction.
        construction: DerivedWitnessIdV1,
    },
}

impl GeometricCategoryV1 {
    const fn tag(self) -> u8 {
        match self {
            Self::Semialgebraic => 0,
            Self::Algebraic => 1,
            Self::RestrictedAnalytic => 2,
            Self::Subanalytic { .. } => 3,
        }
    }
}

/// Locality of an object. Unbounded and infinite-dimensional variants are
/// decodable refusal states, never admitted v1 values.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LocalityScopeV1 {
    /// Germ at an explicitly identified point in one chart.
    GermAt {
        /// Chart containing the base point.
        chart: ConfigurationChartIdV1,
        /// Exact/certified base-point artifact.
        point: DerivedWitnessIdV1,
    },
    /// Relatively compact local neighborhood with retained closure witness.
    CompactNeighborhood {
        /// Chart containing the neighborhood.
        chart: ConfigurationChartIdV1,
        /// Closure/containment witness.
        witness: DerivedWitnessIdV1,
    },
    /// Complete compact finite object.
    GlobalCompact {
        /// Compactness witness.
        witness: DerivedWitnessIdV1,
    },
    /// Explicitly unbounded draft input; v1 refuses it.
    GlobalUnbounded,
    /// Infinite-dimensional formal model; v1 refuses it.
    InfiniteDimensional,
}

/// Compactness evidence attached to charts, strata, and local objects.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CompactnessV1 {
    /// Global compactness was externally established and retained.
    Proved {
        /// Exact witness artifact.
        witness: DerivedWitnessIdV1,
    },
    /// A local model has relatively compact closure in its containing chart.
    RelativelyCompact {
        /// Exact closure/containment witness.
        witness: DerivedWitnessIdV1,
    },
    /// Assumption only; retained for honest decoding but refused by admission.
    Assumed {
        /// Explicit no-claim artifact.
        no_claim: DerivedNoClaimIdV1,
    },
    /// Known unbounded object; refused.
    Unbounded,
    /// No compactness information; refused.
    Unknown,
}

/// Regularity class for an encoded germ or stratum.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RegularityClassV1 {
    /// Polynomial regularity.
    Polynomial,
    /// Real/complex analytic regularity in the declared category.
    Analytic,
    /// Finite C^k regularity.
    Differentiable {
        /// Nonzero differentiability order.
        order: u16,
    },
    /// Missing regularity; refused by admission.
    Unknown,
}

/// Finite computability envelope for a chart/function/complex.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FiniteComputabilityV1 {
    /// Exact finite evaluation through a pinned safe kernel.
    ExactFinite {
        /// Kernel/version witness.
        kernel: DerivedWitnessIdV1,
    },
    /// Finite outward interval evaluation.
    IntervalFinite {
        /// Primitive set and enclosure proof.
        enclosure: DerivedWitnessIdV1,
    },
    /// Finite truncation with a retained remainder enclosure.
    TruncatedFinite {
        /// Resolution identity.
        resolution: DerivedResolutionIdV1,
        /// Remainder witness.
        remainder: DerivedWitnessIdV1,
    },
    /// External opaque evaluation; refused.
    ExternalOpaque {
        /// Explicit no-claim artifact.
        no_claim: DerivedNoClaimIdV1,
    },
    /// Unbounded/infinite computation; refused.
    Infinite,
}

/// Machine-readable local function encoding.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LocalFunctionEncodingV1 {
    /// Exact finite multivariate polynomial.
    Polynomial {
        /// Polynomial content identity.
        polynomial: PolynomialIdV1,
        /// Number of coordinate variables.
        variables: u32,
        /// Total degree.
        degree: u32,
    },
    /// Restricted analytic program over an admitted primitive set.
    RestrictedAnalytic {
        /// Program content identity.
        program: AnalyticProgramIdV1,
        /// Primitive-set admission witness.
        primitives: DerivedWitnessIdV1,
        /// Derivatives retained through this order.
        derivative_order: u16,
    },
    /// Opaque analytic or black-box callback; refused.
    OpaqueExternal {
        /// Exact no-claim artifact.
        no_claim: DerivedNoClaimIdV1,
    },
}

/// Semantic binding of one scalar/vector output to units.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct UnitBindingV1 {
    /// Unit-system identity.
    pub system: DerivedUnitSystemIdV1,
    /// Semantic quantity kind, distinct even for dimensionally equal values.
    pub quantity: DerivedQuantityKindIdV1,
    /// Finite positive scale to the system's canonical unit.
    pub scale_to_canonical: f64,
}

/// Admitted chart family. This is not the runtime [`crate::Chart`] trait: it
/// is a finite mathematical presentation carried by the theorem IR.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConfigurationChartClassV1 {
    /// Smooth finite-dimensional manifold chart.
    SmoothManifold,
    /// Semialgebraic chart.
    Semialgebraic,
    /// Algebraic chart.
    Algebraic,
    /// Restricted analytic chart.
    RestrictedAnalytic,
    /// Chart already equipped with the declared finite stratification.
    Stratified,
}

/// One finite configuration-manifold/chart presentation.
#[derive(Debug, Clone, PartialEq)]
pub struct ConfigurationChartV1 {
    /// Stable chart identity.
    pub id: ConfigurationChartIdV1,
    /// Presentation family.
    pub class: ConfigurationChartClassV1,
    /// Finite coordinate dimension.
    pub coordinate_dimension: u32,
    /// Finite ambient dimension.
    pub ambient_dimension: u32,
    /// Frame convention.
    pub frame: DerivedFrameIdV1,
    /// Coordinate-unit binding.
    pub coordinates: UnitBindingV1,
    /// Scope of the presentation.
    pub locality: LocalityScopeV1,
    /// Compactness evidence.
    pub compactness: CompactnessV1,
    /// Regularity of the chart maps.
    pub regularity: RegularityClassV1,
    /// Finite computation contract.
    pub computability: FiniteComputabilityV1,
}

/// Equality-constraint germ. Its ID cannot be used as any other constraint.
#[derive(Debug, Clone, PartialEq)]
pub struct EqualityConstraintGermV1 {
    /// Stable equality identity.
    pub id: EqualityConstraintIdV1,
    /// Chart on which the germ is defined.
    pub chart: ConfigurationChartIdV1,
    /// Positive codomain dimension.
    pub codomain_dimension: u32,
    /// Exact/interval-admitted equation encoding.
    pub equation: LocalFunctionEncodingV1,
    /// Regularity of the equation map.
    pub regularity: RegularityClassV1,
    /// Equation-output quantity.
    pub units: UnitBindingV1,
    /// Finite computation contract.
    pub computability: FiniteComputabilityV1,
}

/// Sign convention for a real inequality.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InequalitySenseV1 {
    /// Admissible where g(q) >= 0.
    NonNegative,
    /// Admissible where g(q) <= 0.
    NonPositive,
}

/// Explicit active-set state. Candidate activity is not silently promoted.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ActiveSetStateV1 {
    /// Strictly inactive with a retained separation witness.
    Inactive {
        /// Sign-separation witness.
        witness: DerivedWitnessIdV1,
    },
    /// Active at the local model.
    Active {
        /// Zero/sign and active-set witness.
        witness: DerivedWitnessIdV1,
    },
    /// Numerically or semantically unresolved candidate.
    Candidate {
        /// Explicit no-claim artifact.
        no_claim: DerivedNoClaimIdV1,
    },
}

/// Normal-cone presentation for an active inequality or contact.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NormalConeClassV1 {
    /// One generated ray.
    Ray,
    /// Finite polyhedral cone.
    Polyhedral {
        /// Number of nonzero finite generators.
        generators: u32,
    },
    /// Smooth dual cone described by an admitted map.
    SmoothDual,
    /// Omitted/unknown cone; refused for active constraints.
    Unknown,
}

/// One inequality germ with ordering and active-set semantics.
#[derive(Debug, Clone, PartialEq)]
pub struct InequalityConstraintGermV1 {
    /// Stable inequality identity.
    pub id: InequalityConstraintIdV1,
    /// Chart on which the germ is defined.
    pub chart: ConfigurationChartIdV1,
    /// Ordered sign convention.
    pub sense: InequalitySenseV1,
    /// Scalar gap/inequality function.
    pub function: LocalFunctionEncodingV1,
    /// Active-set state.
    pub state: ActiveSetStateV1,
    /// Normal-cone representation.
    pub normal_cone: NormalConeClassV1,
    /// Function-output quantity.
    pub units: UnitBindingV1,
    /// Finite computation contract.
    pub computability: FiniteComputabilityV1,
}

/// Orientation of a relative boundary.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BoundaryOrientationV1 {
    /// Outward orientation.
    Outward,
    /// Inward orientation.
    Inward,
    /// No orientation theorem is claimed.
    Unoriented {
        /// Explicit no-claim artifact.
        no_claim: DerivedNoClaimIdV1,
    },
}

/// Relative boundary connecting an incident parent and boundary stratum.
#[derive(Debug, Clone, PartialEq)]
pub struct RelativeBoundaryV1 {
    /// Stable boundary identity.
    pub id: RelativeBoundaryIdV1,
    /// Chart containing both strata.
    pub chart: ConfigurationChartIdV1,
    /// Higher-dimensional stratum.
    pub parent: StratumIdV1,
    /// Boundary stratum.
    pub boundary: StratumIdV1,
    /// Orientation convention.
    pub orientation: BoundaryOrientationV1,
    /// Incidence/relative-boundary witness.
    pub witness: DerivedWitnessIdV1,
    /// Boundary-defining quantity.
    pub units: UnitBindingV1,
}

/// Contact-law class, kept separate from generic inequalities.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum ContactLawV1 {
    /// Frictionless unilateral contact.
    Frictionless,
    /// Coulomb cone with a finite nonnegative coefficient. Negative zero is a
    /// noncanonical input and is refused.
    Coulomb {
        /// Friction coefficient.
        friction_coefficient: f64,
    },
    /// Set-valued relay/contact relation with retained graph witness.
    SetValued {
        /// Closed-graph/constitutive witness.
        graph: DerivedWitnessIdV1,
    },
}

/// Typed unilateral contact between two relative boundaries.
#[derive(Debug, Clone, PartialEq)]
pub struct ContactConstraintV1 {
    /// Stable contact identity.
    pub id: ContactConstraintIdV1,
    /// Chart containing the local contact model.
    pub chart: ConfigurationChartIdV1,
    /// First relative boundary.
    pub side_a: RelativeBoundaryIdV1,
    /// Second relative boundary.
    pub side_b: RelativeBoundaryIdV1,
    /// Ordered signed gap; contact is not an equality relation.
    pub gap: LocalFunctionEncodingV1,
    /// Active-set state.
    pub state: ActiveSetStateV1,
    /// Normal-cone class.
    pub normal_cone: NormalConeClassV1,
    /// Contact/relay law.
    pub law: ContactLawV1,
    /// Gap units.
    pub units: UnitBindingV1,
    /// Finite computation contract.
    pub computability: FiniteComputabilityV1,
}

/// Constitutive semantic role. It is not a geometric constraint.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConstitutiveRoleV1 {
    /// Stored-energy relation.
    Energy,
    /// Dissipation potential or monotone relation.
    Dissipation,
    /// General finite state relation with no energy claim.
    GeneralRelation,
}

/// Typed constitutive datum retained beside, but never confused with,
/// geometry constraints.
#[derive(Debug, Clone, PartialEq)]
pub struct ConstitutiveDatumV1 {
    /// Stable constitutive identity.
    pub id: ConstitutiveDatumIdV1,
    /// Chart/state neighborhood where the law is admitted.
    pub chart: ConfigurationChartIdV1,
    /// Semantic role.
    pub role: ConstitutiveRoleV1,
    /// Positive finite state dimension.
    pub state_dimension: u32,
    /// Law encoding.
    pub law: LocalFunctionEncodingV1,
    /// Output quantity.
    pub units: UnitBindingV1,
    /// Finite computation contract.
    pub computability: FiniteComputabilityV1,
}

/// Role of a finite graded complex in one local model.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DerivedComplexRoleV1 {
    /// Infinitesimal motions and linearized constraints.
    Tangent,
    /// Dual/reaction and cotangent semantics.
    Cotangent,
    /// Deformation-obstruction/resolution data.
    DeformationObstruction,
}

/// One finite graded vector/module space.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GradedSpaceV1 {
    /// Homological/cohomological degree.
    pub degree: i16,
    /// Finite basis dimension; zero spaces remain explicit.
    pub dimension: u32,
    /// Quantity kind carried by basis coefficients.
    pub quantity: DerivedQuantityKindIdV1,
}

/// One differential in a finite complex.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ComplexDifferentialV1 {
    /// Source degree.
    pub from_degree: i16,
    /// Target degree; must be source + 1.
    pub to_degree: i16,
    /// Exact/interval linear-map artifact.
    pub map: DerivedLinearMapIdV1,
    /// Witness for the relevant adjacent-composition square-zero obligation.
    pub square_zero_witness: DerivedWitnessIdV1,
}

/// Explicit finite resolution/truncation envelope.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FiniteResolutionV1 {
    /// Resolution artifact.
    pub id: DerivedResolutionIdV1,
    /// Lowest retained degree.
    pub min_degree: i16,
    /// Highest retained degree.
    pub max_degree: i16,
    /// Maximum admitted basis dimension in one degree. This declaration must
    /// not exceed [`DERIVED_GEOMETRY_HARD_MAX_DIMENSION_V1`].
    pub max_basis_dimension: u32,
    /// Zero means no series truncation; positive values require remainder.
    pub truncation_order: u32,
    /// Enclosure of discarded terms, when truncated.
    pub remainder: Option<DerivedWitnessIdV1>,
}

/// Finite tangent/cotangent/deformation complex.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct FiniteDerivedComplexV1 {
    /// Stable complex identity.
    pub id: DerivedComplexIdV1,
    /// Chart at which the complex is defined.
    pub chart: ConfigurationChartIdV1,
    /// Nominal complex role.
    pub role: DerivedComplexRoleV1,
    /// Explicit finite graded spaces.
    pub spaces: Vec<GradedSpaceV1>,
    /// Explicit degree-raising differentials.
    pub differentials: Vec<ComplexDifferentialV1>,
    /// Finite truncation/resolution.
    pub resolution: FiniteResolutionV1,
    /// Evaluation/assembly computability.
    pub computability: FiniteComputabilityV1,
}

/// Canonical local singular/regular mechanism class.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DerivedLocalModelClassV1 {
    /// Regular complete-intersection linkage/mechanism chart.
    RegularCompleteIntersection,
    /// Presentation with redundant equations retained explicitly.
    RedundantPresentation,
    /// Plane-cusp type local singularity.
    Cusp,
    /// Node/crossing type local singularity.
    Node,
    /// Relative-boundary half-space model.
    Boundary,
    /// Active inequality/contact corner.
    ContactCorner,
    /// General finite derived local model under the declared resolution.
    GeneralFiniteDerived,
}

/// Scope of presentation/equivalence authority. RD.1b supplies actual maps.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PresentationScopeV1 {
    /// Literal presentation only; no coordinate or generator equivalence.
    Literal {
        /// Explicit no-equivalence claim.
        no_claim: DerivedNoClaimIdV1,
    },
    /// Equivalence is asserted only at one retained finite resolution.
    FixedResolution {
        /// Exact resolution.
        resolution: DerivedResolutionIdV1,
        /// Scope witness.
        witness: DerivedWitnessIdV1,
    },
    /// External theorem metadata exists, but admission does not verify it.
    ExternallyChecked {
        /// Retained equivalence theorem artifact.
        witness: DerivedWitnessIdV1,
    },
}

/// One local derived model with typed references into all object families.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DerivedLocalModelV1 {
    /// Stable local-model identity.
    pub id: DerivedLocalModelIdV1,
    /// Containing configuration chart.
    pub chart: ConfigurationChartIdV1,
    /// Regular/singular local class.
    pub class: DerivedLocalModelClassV1,
    /// Equality germs participating in the presentation.
    pub equalities: Vec<EqualityConstraintIdV1>,
    /// Explicitly active inequalities only.
    pub active_inequalities: Vec<InequalityConstraintIdV1>,
    /// Explicitly active contacts only.
    pub active_contacts: Vec<ContactConstraintIdV1>,
    /// Constitutive metadata used for physical interpretation, not geometry.
    pub constitutive_data: Vec<ConstitutiveDatumIdV1>,
    /// Tangent complex.
    pub tangent_complex: DerivedComplexIdV1,
    /// Cotangent/dual complex.
    pub cotangent_complex: DerivedComplexIdV1,
    /// Deformation-obstruction complex/resolution.
    pub deformation_complex: DerivedComplexIdV1,
    /// Declared virtual dimension, which may differ from tangent nullity at a
    /// singular or redundant presentation.
    pub virtual_dimension: i32,
    /// Exact locality of this model.
    pub locality: LocalityScopeV1,
    /// Current presentation/equivalence scope.
    pub presentation: PresentationScopeV1,
}

/// One finite stratum.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StratumSpecV1 {
    /// Stable stratum identity.
    pub id: StratumIdV1,
    /// Containing chart.
    pub chart: ConfigurationChartIdV1,
    /// Local model governing this stratum.
    pub local_model: DerivedLocalModelIdV1,
    /// Finite geometric dimension.
    pub dimension: u32,
    /// Active inequalities on this stratum.
    pub active_inequalities: Vec<InequalityConstraintIdV1>,
    /// Active contacts on this stratum.
    pub active_contacts: Vec<ContactConstraintIdV1>,
    /// Optional relative-boundary classification.
    pub relative_boundary: Option<RelativeBoundaryIdV1>,
    /// Stratum regularity.
    pub regularity: RegularityClassV1,
    /// Stratum compactness.
    pub compactness: CompactnessV1,
}

/// One incidence relation in the finite stratum poset.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StratumIncidenceV1 {
    /// Lower-dimensional incident stratum.
    pub lower: StratumIdV1,
    /// Higher-dimensional containing stratum.
    pub upper: StratumIdV1,
    /// Exact dimension difference.
    pub codimension: u32,
    /// Incidence/frontier witness.
    pub witness: DerivedWitnessIdV1,
}

/// Local-link topology status.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum LocalLinkTopologyV1 {
    /// Finite simplicial/CW link with retained construction witness.
    FiniteComplex {
        /// Resolution/complex artifact.
        resolution: DerivedResolutionIdV1,
        /// Link-construction witness.
        witness: DerivedWitnessIdV1,
    },
    /// Topology unresolved; retained as an honest no-claim.
    Unknown {
        /// Explicit no-claim artifact.
        no_claim: DerivedNoClaimIdV1,
    },
}

/// Compact local link for one incident stratum pair.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct LocalLinkV1 {
    /// Stable link identity.
    pub id: LocalLinkIdV1,
    /// Lower stratum whose link is represented.
    pub stratum: StratumIdV1,
    /// Higher ambient stratum.
    pub ambient_stratum: StratumIdV1,
    /// Link dimension; admission requires `dimension + 1 == codimension`.
    pub dimension: u32,
    /// Compact-link witness.
    pub compactness_witness: DerivedWitnessIdV1,
    /// Topology state.
    pub topology: LocalLinkTopologyV1,
}

/// Declared stratification theorem level. Admission checks only finite
/// incidence and typing; it does not prove Whitney/Thom conditions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StratificationClassV1 {
    /// Finite incidence stratification, no Whitney claim.
    FiniteIncidence,
    /// Retained Whitney-A theorem metadata.
    WhitneyA {
        /// External theorem witness.
        witness: DerivedWitnessIdV1,
    },
    /// Retained Whitney-B theorem metadata.
    WhitneyB {
        /// External theorem witness.
        witness: DerivedWitnessIdV1,
    },
    /// Retained Thom condition metadata.
    Thom {
        /// External theorem witness.
        witness: DerivedWitnessIdV1,
    },
}

/// Complete finite stratification and local-link payload.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct StratificationV1 {
    /// Stable stratification identity.
    pub id: StratificationIdV1,
    /// Claimed theorem level (metadata only in RD.1a).
    pub class: StratificationClassV1,
    /// Finite strata.
    pub strata: Vec<StratumSpecV1>,
    /// Finite incidence poset edges.
    pub incidences: Vec<StratumIncidenceV1>,
    /// Finite compact local links.
    pub local_links: Vec<LocalLinkV1>,
}

/// Scope of a retained external proof/check record.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DerivedProofScopeV1 {
    /// Whole RD.1a structural object.
    Object,
    /// One local model.
    LocalModel(DerivedLocalModelIdV1),
    /// Whole stratification.
    Stratification(StratificationIdV1),
}

/// Proof-state metadata. `ExternallyChecked` is a retained assertion, not an
/// admitted authority capability; RD.1c independently checks it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DerivedProofStateV1 {
    /// Structural admission only.
    StructuralNoClaim {
        /// Explicit theorem no-claim.
        no_claim: DerivedNoClaimIdV1,
    },
    /// External theorem/check receipt retained with exact identities.
    ExternallyChecked {
        /// Theorem card.
        theorem: DerivedTheoremIdV1,
        /// Checker implementation.
        checker: DerivedCheckerIdV1,
        /// External receipt bytes/content.
        receipt: DerivedWitnessIdV1,
        /// Exact claim scope.
        scope: DerivedProofScopeV1,
    },
}

/// Versioned raw RD.1a object IR.
#[derive(Debug, Clone, PartialEq)]
pub struct DerivedGeometryIrV1 {
    /// Decoded schema version.
    pub schema_version: u32,
    /// Physical/machine subject.
    pub subject: DerivedSubjectIdV1,
    /// Immutable source-model version.
    pub model_version: DerivedModelVersionIdV1,
    /// Mathematical category.
    pub category: GeometricCategoryV1,
    /// Coefficient and order semantics.
    pub coefficients: CoefficientSystemV1,
    /// Global frame convention.
    pub frame: DerivedFrameIdV1,
    /// Global unit-system convention.
    pub unit_system: DerivedUnitSystemIdV1,
    /// Global locality/compactness scope.
    pub locality: LocalityScopeV1,
    /// Global compactness evidence.
    pub compactness: CompactnessV1,
    /// Configuration charts.
    pub charts: Vec<ConfigurationChartV1>,
    /// Equality germs.
    pub equalities: Vec<EqualityConstraintGermV1>,
    /// Inequality germs.
    pub inequalities: Vec<InequalityConstraintGermV1>,
    /// Relative boundaries.
    pub boundaries: Vec<RelativeBoundaryV1>,
    /// Unilateral/contact relations.
    pub contacts: Vec<ContactConstraintV1>,
    /// Physical constitutive metadata.
    pub constitutive_data: Vec<ConstitutiveDatumV1>,
    /// Finite tangent/cotangent/deformation complexes.
    pub complexes: Vec<FiniteDerivedComplexV1>,
    /// Finite derived local models.
    pub local_models: Vec<DerivedLocalModelV1>,
    /// Finite stratification and links.
    pub stratification: StratificationV1,
    /// Honest proof state.
    pub proof_state: DerivedProofStateV1,
}

/// Explicit resource budget for admission and identity construction.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DerivedAdmissionBudgetV1 {
    /// Maximum items in each top-level or nested collection.
    pub max_objects: usize,
    /// Maximum sum of all graded-space dimensions.
    pub max_total_rank: u64,
    /// Maximum canonical frame bytes.
    pub max_canonical_bytes: u64,
    /// Maximum bytes in one canonical field/item.
    pub max_field_bytes: u64,
}

impl DerivedAdmissionBudgetV1 {
    /// Conservative default below every hard ceiling.
    pub const STANDARD: Self = Self {
        max_objects: 1024,
        max_total_rank: 1 << 20,
        max_canonical_bytes: 1 << 22,
        max_field_bytes: 1 << 21,
    };

    fn canonical_limits(self) -> CanonicalLimits {
        CanonicalLimits::new(
            self.max_canonical_bytes,
            self.max_field_bytes,
            32,
            self.max_objects as u64,
            4096,
        )
    }
}

/// Stable collection/object families used in refusal reports.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DerivedObjectKindV1 {
    /// Configuration charts.
    Chart,
    /// Equality germs.
    Equality,
    /// Inequality germs.
    Inequality,
    /// Relative boundaries.
    Boundary,
    /// Contacts.
    Contact,
    /// Constitutive data.
    Constitutive,
    /// Finite complexes.
    Complex,
    /// Local models.
    LocalModel,
    /// Strata.
    Stratum,
    /// Incidences.
    Incidence,
    /// Local links.
    LocalLink,
    /// Proof/category/global context.
    Global,
}

/// Structured fail-closed RD.1a admission issue.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DerivedAdmissionIssueV1 {
    /// Unsupported decoded schema version.
    UnsupportedSchemaVersion {
        /// Supplied version.
        found: u32,
        /// Sole supported version.
        supported: u32,
    },
    /// Invalid caller resource budget.
    InvalidBudget {
        /// Stable budget field.
        field: &'static str,
    },
    /// A collection exceeded its budget before sorting/allocation work.
    ResourceLimit {
        /// Collection/resource family.
        kind: DerivedObjectKindV1,
        /// Requested amount or the first fail-fast lower bound proving excess.
        requested: u64,
        /// Budget/hard limit.
        limit: u64,
    },
    /// Required collection is empty.
    EmptyCollection {
        /// Collection family.
        kind: DerivedObjectKindV1,
    },
    /// Required typed identity is the all-zero sentinel.
    MissingIdentity {
        /// Object family.
        kind: DerivedObjectKindV1,
        /// Stable field.
        field: &'static str,
    },
    /// Stable ID is duplicated within its nominal family.
    DuplicateIdentity {
        /// Object family.
        kind: DerivedObjectKindV1,
    },
    /// Dimension/rank/truncation data are invalid or inconsistent.
    InvalidDimension {
        /// Object family.
        kind: DerivedObjectKindV1,
        /// Stable field.
        field: &'static str,
    },
    /// Chart/object frame differs without an RD.1b transition map.
    MixedFrame {
        /// Object family.
        kind: DerivedObjectKindV1,
    },
    /// Unit systems differ without an RD.1b unit map.
    MixedUnitSystem {
        /// Object family.
        kind: DerivedObjectKindV1,
    },
    /// Unit scale is not finite and strictly positive.
    InvalidUnitScale {
        /// Object family.
        kind: DerivedObjectKindV1,
    },
    /// Mathematical category and coefficient system disagree.
    CategoryCoefficientMismatch,
    /// Inequality/contact/normal-cone semantics require ordered real scalars.
    OrderedSemanticsRequiresReal,
    /// Locality is unbounded or infinite-dimensional.
    UnsupportedLocality {
        /// Object family.
        kind: DerivedObjectKindV1,
    },
    /// Compactness is assumed, unknown, or unbounded.
    UnsupportedCompactness {
        /// Object family.
        kind: DerivedObjectKindV1,
    },
    /// Regularity is unknown or inconsistent with the encoding.
    UnsupportedRegularity {
        /// Object family.
        kind: DerivedObjectKindV1,
    },
    /// Computability is opaque or infinite.
    UnsupportedComputability {
        /// Object family.
        kind: DerivedObjectKindV1,
    },
    /// Function encoding is opaque or categorically incompatible.
    UnsupportedFunctionEncoding {
        /// Object family.
        kind: DerivedObjectKindV1,
    },
    /// Typed reference is absent or points to the wrong family.
    MissingReference {
        /// Object family owning the reference.
        kind: DerivedObjectKindV1,
        /// Stable reference field.
        field: &'static str,
    },
    /// A typed reference resolves, but crosses chart ownership without an
    /// admitted RD.1b transition map.
    CrossChartReference {
        /// Object family owning the reference.
        kind: DerivedObjectKindV1,
        /// Stable reference field.
        field: &'static str,
    },
    /// A local model calls a candidate/inactive constraint active.
    ActiveSetMismatch {
        /// Inequality or contact family.
        kind: DerivedObjectKindV1,
    },
    /// A local model bound a complex under the wrong role or chart.
    ComplexRoleMismatch {
        /// Expected role field.
        field: &'static str,
    },
    /// Finite-complex structure is malformed.
    InvalidComplex {
        /// Stable reason field.
        field: &'static str,
    },
    /// Stratum incidence/frontier data are malformed.
    InvalidStratification {
        /// Stable reason field.
        field: &'static str,
    },
    /// Local-link dimensions or incidence binding are invalid.
    InvalidLocalLink {
        /// Stable reason field.
        field: &'static str,
    },
    /// Proof metadata names an absent local object or zero artifact.
    InvalidProofState,
    /// Bounded cooperative cancellation observed; no token was published.
    Cancelled {
        /// Stable validation stage.
        stage: &'static str,
        /// Fully processed deterministic work units in that stage.
        completed: usize,
    },
    /// Canonical identity construction failed.
    Identity(CanonicalError),
}

/// Deterministic fail-closed refusal report.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DerivedAdmissionReportV1 {
    issues: Vec<DerivedAdmissionIssueV1>,
}

impl DerivedAdmissionReportV1 {
    fn one(issue: DerivedAdmissionIssueV1) -> Self {
        Self {
            issues: vec![issue],
        }
    }

    fn new(issues: Vec<DerivedAdmissionIssueV1>) -> Self {
        Self { issues }
    }

    /// Deterministically ordered admission issues.
    #[must_use]
    pub fn issues(&self) -> &[DerivedAdmissionIssueV1] {
        &self.issues
    }
}

impl fmt::Display for DerivedAdmissionReportV1 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "derived machine geometry refused with {} issue(s)",
            self.issues.len()
        )
    }
}

impl core::error::Error for DerivedAdmissionReportV1 {}

/// Sealed, canonical, structurally admitted RD.1a object.
#[derive(Debug, Clone, PartialEq)]
pub struct AdmittedDerivedGeometryV1 {
    ir: DerivedGeometryIrV1,
    budget: DerivedAdmissionBudgetV1,
    receipt: IdentityReceipt<DerivedGeometryIdV1>,
}

impl AdmittedDerivedGeometryV1 {
    /// Canonically ordered object IR. This is structural data, not detachable
    /// theorem authority.
    #[must_use]
    pub const fn ir(&self) -> &DerivedGeometryIrV1 {
        &self.ir
    }

    /// Explicit admission budget retained beside the object.
    #[must_use]
    pub const fn budget(&self) -> DerivedAdmissionBudgetV1 {
        self.budget
    }

    /// Typed semantic identity.
    #[must_use]
    pub const fn id(&self) -> DerivedGeometryIdV1 {
        self.receipt.id()
    }

    /// Content-addressed canonical preimage and construction limits.
    #[must_use]
    pub const fn admission_receipt(&self) -> IdentityReceipt<DerivedGeometryIdV1> {
        self.receipt
    }
}

/// Validate, canonicalize, and content-address one finite RD.1a object.
///
/// Collection counts and aggregate rank are checked before sorting; dimensional
/// and canonical-byte ceilings are enforced in validation and identity
/// construction. Preflight, canonicalization, validation, and canonical item
/// loops poll the supplied execution context; each uninterruptible sort or
/// nested item loop remains bounded by the admitted collection ceiling. Canonical
/// identity construction adapts the same context to the streaming encoder. No
/// partial admitted value escapes on error or cancellation.
///
/// # Errors
/// Returns [`DerivedAdmissionReportV1`] for any unsupported category/scope,
/// dimensional or typed-reference defect, malformed finite complex,
/// stratification/link defect, budget violation, cancellation, or identity
/// failure.
#[must_use = "derived geometry admission must be handled before theorem use"]
pub fn admit_derived_geometry_v1(
    mut ir: DerivedGeometryIrV1,
    budget: DerivedAdmissionBudgetV1,
    cx: &Cx<'_>,
) -> Result<AdmittedDerivedGeometryV1, DerivedAdmissionReportV1> {
    validate_budget(budget)?;
    preflight_collections(&ir, budget, cx)?;
    canonicalize_ir(&mut ir, cx)?;

    let mut issues = Vec::new();
    if ir.schema_version != DERIVED_GEOMETRY_SCHEMA_VERSION_V1 {
        issues.push(DerivedAdmissionIssueV1::UnsupportedSchemaVersion {
            found: ir.schema_version,
            supported: DERIVED_GEOMETRY_SCHEMA_VERSION_V1,
        });
    }
    let issues_before_duplicate_check = issues.len();
    collect_duplicate_identity_issues(&ir, cx, &mut issues)?;
    if issues.len() != issues_before_duplicate_check {
        return Err(DerivedAdmissionReportV1::new(issues));
    }
    validate_global_context(&ir, &mut issues);
    validate_charts(&ir, cx, &mut issues)?;
    validate_constraints(&ir, cx, &mut issues)?;
    validate_complexes(&ir, cx, &mut issues)?;
    validate_local_models(&ir, cx, &mut issues)?;
    validate_stratification(&ir, cx, &mut issues)?;
    validate_proof_state(&ir, &mut issues);
    if !issues.is_empty() {
        return Err(DerivedAdmissionReportV1::new(issues));
    }

    checkpoint(cx, "identity", 0)?;
    let receipt = derived_geometry_receipt(&ir, budget, cx)
        .map_err(|error| DerivedAdmissionReportV1::one(DerivedAdmissionIssueV1::Identity(error)))?;
    checkpoint(cx, "publish", 0)?;
    Ok(AdmittedDerivedGeometryV1 {
        ir,
        budget,
        receipt,
    })
}

fn checkpoint(
    cx: &Cx<'_>,
    stage: &'static str,
    completed: usize,
) -> Result<(), DerivedAdmissionReportV1> {
    cx.checkpoint().map_err(|_| {
        DerivedAdmissionReportV1::one(DerivedAdmissionIssueV1::Cancelled { stage, completed })
    })
}

fn validate_budget(budget: DerivedAdmissionBudgetV1) -> Result<(), DerivedAdmissionReportV1> {
    let invalid =
        if budget.max_objects == 0 || budget.max_objects > DERIVED_GEOMETRY_HARD_MAX_OBJECTS_V1 {
            Some("max_objects")
        } else if budget.max_total_rank == 0
            || budget.max_total_rank > DERIVED_GEOMETRY_HARD_MAX_TOTAL_RANK_V1
        {
            Some("max_total_rank")
        } else if budget.max_canonical_bytes == 0
            || budget.max_canonical_bytes > DERIVED_GEOMETRY_HARD_MAX_CANONICAL_BYTES_V1
        {
            Some("max_canonical_bytes")
        } else if budget.max_field_bytes == 0
            || budget.max_field_bytes > DERIVED_GEOMETRY_HARD_MAX_FIELD_BYTES_V1
            || budget.max_field_bytes > budget.max_canonical_bytes
        {
            Some("max_field_bytes")
        } else {
            None
        };
    invalid.map_or(Ok(()), |field| {
        Err(DerivedAdmissionReportV1::one(
            DerivedAdmissionIssueV1::InvalidBudget { field },
        ))
    })
}

fn enforce_count(
    kind: DerivedObjectKindV1,
    found: usize,
    budget: DerivedAdmissionBudgetV1,
) -> Result<(), DerivedAdmissionReportV1> {
    if found > budget.max_objects {
        return Err(DerivedAdmissionReportV1::one(
            DerivedAdmissionIssueV1::ResourceLimit {
                kind,
                requested: budget.max_objects as u64 + 1,
                limit: budget.max_objects as u64,
            },
        ));
    }
    Ok(())
}

fn preflight_collections(
    ir: &DerivedGeometryIrV1,
    budget: DerivedAdmissionBudgetV1,
    cx: &Cx<'_>,
) -> Result<(), DerivedAdmissionReportV1> {
    checkpoint(cx, "preflight", 0)?;
    for (kind, count) in [
        (DerivedObjectKindV1::Chart, ir.charts.len()),
        (DerivedObjectKindV1::Equality, ir.equalities.len()),
        (DerivedObjectKindV1::Inequality, ir.inequalities.len()),
        (DerivedObjectKindV1::Boundary, ir.boundaries.len()),
        (DerivedObjectKindV1::Contact, ir.contacts.len()),
        (
            DerivedObjectKindV1::Constitutive,
            ir.constitutive_data.len(),
        ),
        (DerivedObjectKindV1::Complex, ir.complexes.len()),
        (DerivedObjectKindV1::LocalModel, ir.local_models.len()),
        (DerivedObjectKindV1::Stratum, ir.stratification.strata.len()),
        (
            DerivedObjectKindV1::Incidence,
            ir.stratification.incidences.len(),
        ),
        (
            DerivedObjectKindV1::LocalLink,
            ir.stratification.local_links.len(),
        ),
    ] {
        enforce_count(kind, count, budget)?;
    }

    let mut total_rank = 0u64;
    let mut completed_spaces = 0usize;
    for (completed, complex) in ir.complexes.iter().enumerate() {
        checkpoint(cx, "preflight", completed)?;
        enforce_count(DerivedObjectKindV1::Complex, complex.spaces.len(), budget)?;
        enforce_count(
            DerivedObjectKindV1::Complex,
            complex.differentials.len(),
            budget,
        )?;
        for space in &complex.spaces {
            if completed_spaces.is_multiple_of(DERIVED_GEOMETRY_CANCELLATION_STRIDE_V1) {
                checkpoint(cx, "preflight-rank", completed_spaces)?;
            }
            total_rank = total_rank
                .checked_add(u64::from(space.dimension))
                .ok_or_else(|| {
                    DerivedAdmissionReportV1::one(DerivedAdmissionIssueV1::ResourceLimit {
                        kind: DerivedObjectKindV1::Complex,
                        requested: u64::MAX,
                        limit: budget.max_total_rank,
                    })
                })?;
            completed_spaces += 1;
            if total_rank > budget.max_total_rank {
                return Err(DerivedAdmissionReportV1::one(
                    DerivedAdmissionIssueV1::ResourceLimit {
                        kind: DerivedObjectKindV1::Complex,
                        requested: budget.max_total_rank + 1,
                        limit: budget.max_total_rank,
                    },
                ));
            }
        }
    }
    for (completed, model) in ir.local_models.iter().enumerate() {
        checkpoint(cx, "preflight", completed)?;
        for count in [
            model.equalities.len(),
            model.active_inequalities.len(),
            model.active_contacts.len(),
            model.constitutive_data.len(),
        ] {
            enforce_count(DerivedObjectKindV1::LocalModel, count, budget)?;
        }
    }
    for (completed, stratum) in ir.stratification.strata.iter().enumerate() {
        checkpoint(cx, "preflight", completed)?;
        enforce_count(
            DerivedObjectKindV1::Stratum,
            stratum.active_inequalities.len(),
            budget,
        )?;
        enforce_count(
            DerivedObjectKindV1::Stratum,
            stratum.active_contacts.len(),
            budget,
        )?;
    }
    Ok(())
}

fn canonicalize_ir(
    ir: &mut DerivedGeometryIrV1,
    cx: &Cx<'_>,
) -> Result<(), DerivedAdmissionReportV1> {
    canonicalize_ir_with(ir, |completed| checkpoint(cx, "canonicalize", completed))
}

fn canonicalize_ir_with(
    ir: &mut DerivedGeometryIrV1,
    mut poll: impl FnMut(usize) -> Result<(), DerivedAdmissionReportV1>,
) -> Result<(), DerivedAdmissionReportV1> {
    let mut completed = 0usize;
    poll(completed)?;

    ir.charts.sort_unstable_by_key(|chart| chart.id);
    completed += 1;
    poll(completed)?;
    ir.equalities
        .sort_unstable_by_key(|constraint| constraint.id);
    completed += 1;
    poll(completed)?;
    ir.inequalities
        .sort_unstable_by_key(|constraint| constraint.id);
    completed += 1;
    poll(completed)?;
    ir.boundaries.sort_unstable_by_key(|boundary| boundary.id);
    completed += 1;
    poll(completed)?;
    ir.contacts.sort_unstable_by_key(|contact| contact.id);
    completed += 1;
    poll(completed)?;
    ir.constitutive_data.sort_unstable_by_key(|datum| datum.id);
    completed += 1;
    poll(completed)?;
    for complex in &mut ir.complexes {
        complex.spaces.sort_unstable_by_key(|space| space.degree);
        completed += 1;
        poll(completed)?;
        complex
            .differentials
            .sort_unstable_by_key(|map| (map.from_degree, map.to_degree, map.map));
        completed += 1;
        poll(completed)?;
    }
    ir.complexes.sort_unstable_by_key(|complex| complex.id);
    completed += 1;
    poll(completed)?;
    for model in &mut ir.local_models {
        model.equalities.sort_unstable();
        completed += 1;
        poll(completed)?;
        model.active_inequalities.sort_unstable();
        completed += 1;
        poll(completed)?;
        model.active_contacts.sort_unstable();
        completed += 1;
        poll(completed)?;
        model.constitutive_data.sort_unstable();
        completed += 1;
        poll(completed)?;
    }
    ir.local_models.sort_unstable_by_key(|model| model.id);
    completed += 1;
    poll(completed)?;
    for stratum in &mut ir.stratification.strata {
        stratum.active_inequalities.sort_unstable();
        completed += 1;
        poll(completed)?;
        stratum.active_contacts.sort_unstable();
        completed += 1;
        poll(completed)?;
    }
    ir.stratification
        .strata
        .sort_unstable_by_key(|stratum| stratum.id);
    completed += 1;
    poll(completed)?;
    ir.stratification
        .incidences
        .sort_unstable_by_key(|incidence| (incidence.lower, incidence.upper));
    completed += 1;
    poll(completed)?;
    ir.stratification
        .local_links
        .sort_unstable_by_key(|link| link.id);
    completed += 1;
    poll(completed)
}

fn validate_global_context(ir: &DerivedGeometryIrV1, issues: &mut Vec<DerivedAdmissionIssueV1>) {
    for (missing, field) in [
        (ir.subject.is_zero(), "subject"),
        (ir.model_version.is_zero(), "model_version"),
        (ir.frame.is_zero(), "frame"),
        (ir.unit_system.is_zero(), "unit_system"),
        (ir.stratification.id.is_zero(), "stratification"),
    ] {
        if missing {
            issues.push(DerivedAdmissionIssueV1::MissingIdentity {
                kind: DerivedObjectKindV1::Global,
                field,
            });
        }
    }
    match ir.coefficients {
        CoefficientSystemV1::DyadicIntervalReal { precision_bits: 0 } => {
            issues.push(DerivedAdmissionIssueV1::CategoryCoefficientMismatch);
        }
        CoefficientSystemV1::RationalReal
        | CoefficientSystemV1::AlgebraicReal
        | CoefficientSystemV1::DyadicIntervalReal { .. }
        | CoefficientSystemV1::RationalComplex
        | CoefficientSystemV1::AlgebraicComplex => {}
    }
    if matches!(
        ir.category,
        GeometricCategoryV1::Semialgebraic
            | GeometricCategoryV1::RestrictedAnalytic
            | GeometricCategoryV1::Subanalytic { .. }
    ) && !ir.coefficients.is_ordered_real()
    {
        issues.push(DerivedAdmissionIssueV1::CategoryCoefficientMismatch);
    }
    if let GeometricCategoryV1::Subanalytic { construction } = ir.category
        && construction.is_zero()
    {
        issues.push(DerivedAdmissionIssueV1::MissingIdentity {
            kind: DerivedObjectKindV1::Global,
            field: "subanalytic_construction",
        });
    }
    check_locality(ir.locality, &ir.charts, DerivedObjectKindV1::Global, issues);
    check_compactness(ir.compactness, DerivedObjectKindV1::Global, issues);
    if ir.charts.is_empty() {
        issues.push(DerivedAdmissionIssueV1::EmptyCollection {
            kind: DerivedObjectKindV1::Chart,
        });
    }
    if ir.local_models.is_empty() {
        issues.push(DerivedAdmissionIssueV1::EmptyCollection {
            kind: DerivedObjectKindV1::LocalModel,
        });
    }
    if ir.stratification.strata.is_empty() {
        issues.push(DerivedAdmissionIssueV1::EmptyCollection {
            kind: DerivedObjectKindV1::Stratum,
        });
    }
}

fn has_duplicate<T, K: PartialEq>(items: &[T], key: impl Fn(&T) -> K) -> bool {
    items
        .windows(2)
        .any(|window| key(&window[0]) == key(&window[1]))
}

fn collect_duplicate_identity_issues(
    ir: &DerivedGeometryIrV1,
    cx: &Cx<'_>,
    issues: &mut Vec<DerivedAdmissionIssueV1>,
) -> Result<(), DerivedAdmissionReportV1> {
    scan_duplicate_family(
        cx,
        0,
        &ir.charts,
        |chart| chart.id,
        DerivedObjectKindV1::Chart,
        issues,
    )?;
    scan_duplicate_family(
        cx,
        1,
        &ir.equalities,
        |constraint| constraint.id,
        DerivedObjectKindV1::Equality,
        issues,
    )?;
    scan_duplicate_family(
        cx,
        2,
        &ir.inequalities,
        |constraint| constraint.id,
        DerivedObjectKindV1::Inequality,
        issues,
    )?;
    scan_duplicate_family(
        cx,
        3,
        &ir.boundaries,
        |boundary| boundary.id,
        DerivedObjectKindV1::Boundary,
        issues,
    )?;
    scan_duplicate_family(
        cx,
        4,
        &ir.contacts,
        |contact| contact.id,
        DerivedObjectKindV1::Contact,
        issues,
    )?;
    scan_duplicate_family(
        cx,
        5,
        &ir.constitutive_data,
        |datum| datum.id,
        DerivedObjectKindV1::Constitutive,
        issues,
    )?;
    scan_duplicate_family(
        cx,
        6,
        &ir.complexes,
        |complex| complex.id,
        DerivedObjectKindV1::Complex,
        issues,
    )?;
    scan_duplicate_family(
        cx,
        7,
        &ir.local_models,
        |model| model.id,
        DerivedObjectKindV1::LocalModel,
        issues,
    )?;
    scan_duplicate_family(
        cx,
        8,
        &ir.stratification.strata,
        |stratum| stratum.id,
        DerivedObjectKindV1::Stratum,
        issues,
    )?;
    scan_duplicate_family(
        cx,
        9,
        &ir.stratification.incidences,
        |edge| (edge.lower, edge.upper),
        DerivedObjectKindV1::Incidence,
        issues,
    )?;
    scan_duplicate_family(
        cx,
        10,
        &ir.stratification.local_links,
        |link| link.id,
        DerivedObjectKindV1::LocalLink,
        issues,
    )
}

fn scan_duplicate_family<T, K: PartialEq>(
    cx: &Cx<'_>,
    completed: usize,
    items: &[T],
    key: impl Fn(&T) -> K,
    kind: DerivedObjectKindV1,
    issues: &mut Vec<DerivedAdmissionIssueV1>,
) -> Result<(), DerivedAdmissionReportV1> {
    checkpoint(cx, "duplicate-identities", completed)?;
    if has_duplicate(items, key) {
        issues.push(DerivedAdmissionIssueV1::DuplicateIdentity { kind });
    }
    Ok(())
}

fn unique_by_key<T, K: Copy + Ord>(items: &[T], target: K, key: impl Fn(&T) -> K) -> Option<&T> {
    let first = items.partition_point(|item| key(item) < target);
    let value = items.get(first)?;
    if key(value) != target || items.get(first + 1).is_some_and(|next| key(next) == target) {
        return None;
    }
    Some(value)
}

fn has_chart(ir: &DerivedGeometryIrV1, id: ConfigurationChartIdV1) -> bool {
    unique_by_key(&ir.charts, id, |chart| chart.id).is_some()
}

fn inequality(
    ir: &DerivedGeometryIrV1,
    id: InequalityConstraintIdV1,
) -> Option<&InequalityConstraintGermV1> {
    unique_by_key(&ir.inequalities, id, |constraint| constraint.id)
}

fn boundary(ir: &DerivedGeometryIrV1, id: RelativeBoundaryIdV1) -> Option<&RelativeBoundaryV1> {
    unique_by_key(&ir.boundaries, id, |boundary| boundary.id)
}

fn contact(ir: &DerivedGeometryIrV1, id: ContactConstraintIdV1) -> Option<&ContactConstraintV1> {
    unique_by_key(&ir.contacts, id, |contact| contact.id)
}

fn constitutive(
    ir: &DerivedGeometryIrV1,
    id: ConstitutiveDatumIdV1,
) -> Option<&ConstitutiveDatumV1> {
    unique_by_key(&ir.constitutive_data, id, |datum| datum.id)
}

fn complex(ir: &DerivedGeometryIrV1, id: DerivedComplexIdV1) -> Option<&FiniteDerivedComplexV1> {
    unique_by_key(&ir.complexes, id, |complex| complex.id)
}

fn local_model(
    ir: &DerivedGeometryIrV1,
    id: DerivedLocalModelIdV1,
) -> Option<&DerivedLocalModelV1> {
    unique_by_key(&ir.local_models, id, |model| model.id)
}

fn stratum(ir: &DerivedGeometryIrV1, id: StratumIdV1) -> Option<&StratumSpecV1> {
    unique_by_key(&ir.stratification.strata, id, |stratum| stratum.id)
}

fn check_unit(
    binding: UnitBindingV1,
    expected: DerivedUnitSystemIdV1,
    kind: DerivedObjectKindV1,
    issues: &mut Vec<DerivedAdmissionIssueV1>,
) {
    if binding.system != expected {
        issues.push(DerivedAdmissionIssueV1::MixedUnitSystem { kind });
    }
    if binding.system.is_zero() || binding.quantity.is_zero() {
        issues.push(DerivedAdmissionIssueV1::MissingIdentity {
            kind,
            field: "unit_binding",
        });
    }
    if !binding.scale_to_canonical.is_finite() || binding.scale_to_canonical <= 0.0 {
        issues.push(DerivedAdmissionIssueV1::InvalidUnitScale { kind });
    }
}

fn check_locality(
    locality: LocalityScopeV1,
    charts: &[ConfigurationChartV1],
    kind: DerivedObjectKindV1,
    issues: &mut Vec<DerivedAdmissionIssueV1>,
) {
    match locality {
        LocalityScopeV1::GermAt { chart, point } => {
            if unique_by_key(charts, chart, |candidate| candidate.id).is_none() {
                issues.push(DerivedAdmissionIssueV1::MissingReference {
                    kind,
                    field: "locality_chart",
                });
            }
            if point.is_zero() {
                issues.push(DerivedAdmissionIssueV1::MissingIdentity {
                    kind,
                    field: "locality_point",
                });
            }
        }
        LocalityScopeV1::CompactNeighborhood { chart, witness } => {
            if unique_by_key(charts, chart, |candidate| candidate.id).is_none() {
                issues.push(DerivedAdmissionIssueV1::MissingReference {
                    kind,
                    field: "locality_chart",
                });
            }
            if witness.is_zero() {
                issues.push(DerivedAdmissionIssueV1::MissingIdentity {
                    kind,
                    field: "locality_witness",
                });
            }
        }
        LocalityScopeV1::GlobalCompact { witness } => {
            if witness.is_zero() {
                issues.push(DerivedAdmissionIssueV1::MissingIdentity {
                    kind,
                    field: "locality_witness",
                });
            }
        }
        LocalityScopeV1::GlobalUnbounded | LocalityScopeV1::InfiniteDimensional => {
            issues.push(DerivedAdmissionIssueV1::UnsupportedLocality { kind });
        }
    }
}

fn check_compactness(
    compactness: CompactnessV1,
    kind: DerivedObjectKindV1,
    issues: &mut Vec<DerivedAdmissionIssueV1>,
) {
    match compactness {
        CompactnessV1::Proved { witness } | CompactnessV1::RelativelyCompact { witness } => {
            if witness.is_zero() {
                issues.push(DerivedAdmissionIssueV1::MissingIdentity {
                    kind,
                    field: "compactness_witness",
                });
            }
        }
        CompactnessV1::Assumed { no_claim } => {
            if no_claim.is_zero() {
                issues.push(DerivedAdmissionIssueV1::MissingIdentity {
                    kind,
                    field: "compactness_no_claim",
                });
            }
            issues.push(DerivedAdmissionIssueV1::UnsupportedCompactness { kind });
        }
        CompactnessV1::Unbounded | CompactnessV1::Unknown => {
            issues.push(DerivedAdmissionIssueV1::UnsupportedCompactness { kind });
        }
    }
}

fn check_regularity(
    regularity: RegularityClassV1,
    kind: DerivedObjectKindV1,
    issues: &mut Vec<DerivedAdmissionIssueV1>,
) {
    if matches!(
        regularity,
        RegularityClassV1::Unknown | RegularityClassV1::Differentiable { order: 0 }
    ) {
        issues.push(DerivedAdmissionIssueV1::UnsupportedRegularity { kind });
    }
}

fn check_computability(
    computability: FiniteComputabilityV1,
    kind: DerivedObjectKindV1,
    issues: &mut Vec<DerivedAdmissionIssueV1>,
) {
    match computability {
        FiniteComputabilityV1::ExactFinite { kernel } => {
            if kernel.is_zero() {
                issues.push(DerivedAdmissionIssueV1::MissingIdentity {
                    kind,
                    field: "computability_kernel",
                });
            }
        }
        FiniteComputabilityV1::IntervalFinite { enclosure } => {
            if enclosure.is_zero() {
                issues.push(DerivedAdmissionIssueV1::MissingIdentity {
                    kind,
                    field: "computability_enclosure",
                });
            }
        }
        FiniteComputabilityV1::TruncatedFinite {
            resolution,
            remainder,
        } => {
            if resolution.is_zero() || remainder.is_zero() {
                issues.push(DerivedAdmissionIssueV1::MissingIdentity {
                    kind,
                    field: "computability_truncation",
                });
            }
        }
        FiniteComputabilityV1::ExternalOpaque { no_claim } => {
            if no_claim.is_zero() {
                issues.push(DerivedAdmissionIssueV1::MissingIdentity {
                    kind,
                    field: "computability_no_claim",
                });
            }
            issues.push(DerivedAdmissionIssueV1::UnsupportedComputability { kind });
        }
        FiniteComputabilityV1::Infinite => {
            issues.push(DerivedAdmissionIssueV1::UnsupportedComputability { kind });
        }
    }
}

fn check_function(
    function: LocalFunctionEncodingV1,
    regularity: Option<RegularityClassV1>,
    category: GeometricCategoryV1,
    chart_dimension: u32,
    kind: DerivedObjectKindV1,
    issues: &mut Vec<DerivedAdmissionIssueV1>,
) {
    match function {
        LocalFunctionEncodingV1::Polynomial {
            polynomial,
            variables,
            ..
        } => {
            if polynomial.is_zero() {
                issues.push(DerivedAdmissionIssueV1::MissingIdentity {
                    kind,
                    field: "polynomial",
                });
            }
            if variables == 0 || variables != chart_dimension {
                issues.push(DerivedAdmissionIssueV1::InvalidDimension {
                    kind,
                    field: "function_variables",
                });
            }
            if regularity.is_some_and(|value| value != RegularityClassV1::Polynomial) {
                issues.push(DerivedAdmissionIssueV1::UnsupportedRegularity { kind });
            }
        }
        LocalFunctionEncodingV1::RestrictedAnalytic {
            program,
            primitives,
            derivative_order,
        } => {
            if program.is_zero() || primitives.is_zero() {
                issues.push(DerivedAdmissionIssueV1::MissingIdentity {
                    kind,
                    field: "analytic_program",
                });
            }
            if derivative_order == 0 {
                issues.push(DerivedAdmissionIssueV1::UnsupportedRegularity { kind });
            }
            if regularity.is_some_and(|value| value != RegularityClassV1::Analytic) {
                issues.push(DerivedAdmissionIssueV1::UnsupportedRegularity { kind });
            }
            if !matches!(
                category,
                GeometricCategoryV1::RestrictedAnalytic | GeometricCategoryV1::Subanalytic { .. }
            ) {
                issues.push(DerivedAdmissionIssueV1::UnsupportedFunctionEncoding { kind });
            }
        }
        LocalFunctionEncodingV1::OpaqueExternal { no_claim } => {
            if no_claim.is_zero() {
                issues.push(DerivedAdmissionIssueV1::MissingIdentity {
                    kind,
                    field: "function_no_claim",
                });
            }
            issues.push(DerivedAdmissionIssueV1::UnsupportedFunctionEncoding { kind });
        }
    }
    if matches!(
        category,
        GeometricCategoryV1::Semialgebraic | GeometricCategoryV1::Algebraic
    ) && !matches!(function, LocalFunctionEncodingV1::Polynomial { .. })
    {
        issues.push(DerivedAdmissionIssueV1::UnsupportedFunctionEncoding { kind });
    }
}

fn chart_dimension(ir: &DerivedGeometryIrV1, id: ConfigurationChartIdV1) -> Option<u32> {
    unique_by_key(&ir.charts, id, |chart| chart.id).map(|chart| chart.coordinate_dimension)
}

fn chart_class_matches(category: GeometricCategoryV1, class: ConfigurationChartClassV1) -> bool {
    match category {
        GeometricCategoryV1::Semialgebraic => matches!(
            class,
            ConfigurationChartClassV1::Semialgebraic | ConfigurationChartClassV1::Stratified
        ),
        GeometricCategoryV1::Algebraic => matches!(
            class,
            ConfigurationChartClassV1::Algebraic | ConfigurationChartClassV1::Stratified
        ),
        GeometricCategoryV1::RestrictedAnalytic | GeometricCategoryV1::Subanalytic { .. } => {
            matches!(
                class,
                ConfigurationChartClassV1::RestrictedAnalytic
                    | ConfigurationChartClassV1::Stratified
                    | ConfigurationChartClassV1::SmoothManifold
            )
        }
    }
}

fn validate_charts(
    ir: &DerivedGeometryIrV1,
    cx: &Cx<'_>,
    issues: &mut Vec<DerivedAdmissionIssueV1>,
) -> Result<(), DerivedAdmissionReportV1> {
    if has_duplicate(&ir.charts, |chart| chart.id) {
        issues.push(DerivedAdmissionIssueV1::DuplicateIdentity {
            kind: DerivedObjectKindV1::Chart,
        });
    }
    for (completed, chart) in ir.charts.iter().enumerate() {
        checkpoint(cx, "charts", completed)?;
        if chart.id.is_zero() {
            issues.push(DerivedAdmissionIssueV1::MissingIdentity {
                kind: DerivedObjectKindV1::Chart,
                field: "id",
            });
        }
        if chart.coordinate_dimension == 0
            || chart.ambient_dimension == 0
            || chart.coordinate_dimension > chart.ambient_dimension
            || chart.ambient_dimension > DERIVED_GEOMETRY_HARD_MAX_DIMENSION_V1
        {
            issues.push(DerivedAdmissionIssueV1::InvalidDimension {
                kind: DerivedObjectKindV1::Chart,
                field: "coordinate_or_ambient",
            });
        }
        if !chart_class_matches(ir.category, chart.class) {
            issues.push(DerivedAdmissionIssueV1::UnsupportedFunctionEncoding {
                kind: DerivedObjectKindV1::Chart,
            });
        }
        if chart.frame != ir.frame || chart.frame.is_zero() {
            issues.push(DerivedAdmissionIssueV1::MixedFrame {
                kind: DerivedObjectKindV1::Chart,
            });
        }
        check_unit(
            chart.coordinates,
            ir.unit_system,
            DerivedObjectKindV1::Chart,
            issues,
        );
        check_locality(
            chart.locality,
            &ir.charts,
            DerivedObjectKindV1::Chart,
            issues,
        );
        match chart.locality {
            LocalityScopeV1::GermAt { chart: owner, .. }
            | LocalityScopeV1::CompactNeighborhood { chart: owner, .. }
                if owner != chart.id =>
            {
                if has_chart(ir, owner) {
                    issues.push(DerivedAdmissionIssueV1::CrossChartReference {
                        kind: DerivedObjectKindV1::Chart,
                        field: "self_locality_chart",
                    });
                }
            }
            LocalityScopeV1::GermAt { .. }
            | LocalityScopeV1::CompactNeighborhood { .. }
            | LocalityScopeV1::GlobalCompact { .. }
            | LocalityScopeV1::GlobalUnbounded
            | LocalityScopeV1::InfiniteDimensional => {}
        }
        check_compactness(chart.compactness, DerivedObjectKindV1::Chart, issues);
        check_regularity(chart.regularity, DerivedObjectKindV1::Chart, issues);
        check_computability(chart.computability, DerivedObjectKindV1::Chart, issues);
    }
    Ok(())
}

fn check_active_state(
    state: ActiveSetStateV1,
    kind: DerivedObjectKindV1,
    issues: &mut Vec<DerivedAdmissionIssueV1>,
) {
    let id_missing = match state {
        ActiveSetStateV1::Inactive { witness } | ActiveSetStateV1::Active { witness } => {
            witness.is_zero()
        }
        ActiveSetStateV1::Candidate { no_claim } => no_claim.is_zero(),
    };
    if id_missing {
        issues.push(DerivedAdmissionIssueV1::MissingIdentity {
            kind,
            field: "active_set_witness",
        });
    }
}

fn check_normal_cone(
    cone: NormalConeClassV1,
    state: ActiveSetStateV1,
    kind: DerivedObjectKindV1,
    issues: &mut Vec<DerivedAdmissionIssueV1>,
) {
    if matches!(cone, NormalConeClassV1::Polyhedral { generators: 0 })
        || (matches!(state, ActiveSetStateV1::Active { .. })
            && matches!(cone, NormalConeClassV1::Unknown))
    {
        issues.push(DerivedAdmissionIssueV1::InvalidDimension {
            kind,
            field: "normal_cone",
        });
    }
}

fn validate_constraints(
    ir: &DerivedGeometryIrV1,
    cx: &Cx<'_>,
    issues: &mut Vec<DerivedAdmissionIssueV1>,
) -> Result<(), DerivedAdmissionReportV1> {
    for (duplicate, kind) in [
        (
            has_duplicate(&ir.equalities, |constraint| constraint.id),
            DerivedObjectKindV1::Equality,
        ),
        (
            has_duplicate(&ir.inequalities, |constraint| constraint.id),
            DerivedObjectKindV1::Inequality,
        ),
        (
            has_duplicate(&ir.boundaries, |boundary| boundary.id),
            DerivedObjectKindV1::Boundary,
        ),
        (
            has_duplicate(&ir.contacts, |contact| contact.id),
            DerivedObjectKindV1::Contact,
        ),
        (
            has_duplicate(&ir.constitutive_data, |datum| datum.id),
            DerivedObjectKindV1::Constitutive,
        ),
    ] {
        if duplicate {
            issues.push(DerivedAdmissionIssueV1::DuplicateIdentity { kind });
        }
    }

    for (completed, equality) in ir.equalities.iter().enumerate() {
        checkpoint(cx, "equalities", completed)?;
        if equality.id.is_zero() {
            issues.push(DerivedAdmissionIssueV1::MissingIdentity {
                kind: DerivedObjectKindV1::Equality,
                field: "id",
            });
        }
        let Some(dimension) = chart_dimension(ir, equality.chart) else {
            issues.push(DerivedAdmissionIssueV1::MissingReference {
                kind: DerivedObjectKindV1::Equality,
                field: "chart",
            });
            continue;
        };
        if equality.codomain_dimension == 0
            || equality.codomain_dimension > DERIVED_GEOMETRY_HARD_MAX_DIMENSION_V1
        {
            issues.push(DerivedAdmissionIssueV1::InvalidDimension {
                kind: DerivedObjectKindV1::Equality,
                field: "codomain_dimension",
            });
        }
        check_function(
            equality.equation,
            Some(equality.regularity),
            ir.category,
            dimension,
            DerivedObjectKindV1::Equality,
            issues,
        );
        check_regularity(equality.regularity, DerivedObjectKindV1::Equality, issues);
        check_unit(
            equality.units,
            ir.unit_system,
            DerivedObjectKindV1::Equality,
            issues,
        );
        check_computability(
            equality.computability,
            DerivedObjectKindV1::Equality,
            issues,
        );
    }

    for (completed, inequality) in ir.inequalities.iter().enumerate() {
        checkpoint(cx, "inequalities", completed)?;
        if inequality.id.is_zero() {
            issues.push(DerivedAdmissionIssueV1::MissingIdentity {
                kind: DerivedObjectKindV1::Inequality,
                field: "id",
            });
        }
        let Some(dimension) = chart_dimension(ir, inequality.chart) else {
            issues.push(DerivedAdmissionIssueV1::MissingReference {
                kind: DerivedObjectKindV1::Inequality,
                field: "chart",
            });
            continue;
        };
        check_function(
            inequality.function,
            None,
            ir.category,
            dimension,
            DerivedObjectKindV1::Inequality,
            issues,
        );
        check_active_state(inequality.state, DerivedObjectKindV1::Inequality, issues);
        check_normal_cone(
            inequality.normal_cone,
            inequality.state,
            DerivedObjectKindV1::Inequality,
            issues,
        );
        check_unit(
            inequality.units,
            ir.unit_system,
            DerivedObjectKindV1::Inequality,
            issues,
        );
        check_computability(
            inequality.computability,
            DerivedObjectKindV1::Inequality,
            issues,
        );
    }
    if (!ir.inequalities.is_empty() || !ir.contacts.is_empty())
        && !ir.coefficients.is_ordered_real()
    {
        issues.push(DerivedAdmissionIssueV1::OrderedSemanticsRequiresReal);
    }

    for (completed, boundary) in ir.boundaries.iter().enumerate() {
        checkpoint(cx, "boundaries", completed)?;
        if boundary.id.is_zero() || boundary.witness.is_zero() {
            issues.push(DerivedAdmissionIssueV1::MissingIdentity {
                kind: DerivedObjectKindV1::Boundary,
                field: "id_or_witness",
            });
        }
        if !has_chart(ir, boundary.chart) {
            issues.push(DerivedAdmissionIssueV1::MissingReference {
                kind: DerivedObjectKindV1::Boundary,
                field: "chart",
            });
        }
        let parent = stratum(ir, boundary.parent);
        let boundary_stratum = stratum(ir, boundary.boundary);
        if boundary.parent == boundary.boundary {
            issues.push(DerivedAdmissionIssueV1::MissingReference {
                kind: DerivedObjectKindV1::Boundary,
                field: "strata",
            });
        } else {
            for (stratum, field) in [(parent, "parent"), (boundary_stratum, "boundary")] {
                match stratum {
                    None => issues.push(DerivedAdmissionIssueV1::MissingReference {
                        kind: DerivedObjectKindV1::Boundary,
                        field,
                    }),
                    Some(value) if value.chart != boundary.chart => {
                        issues.push(DerivedAdmissionIssueV1::CrossChartReference {
                            kind: DerivedObjectKindV1::Boundary,
                            field,
                        });
                    }
                    Some(_) => {}
                }
            }
        }
        if let BoundaryOrientationV1::Unoriented { no_claim } = boundary.orientation
            && no_claim.is_zero()
        {
            issues.push(DerivedAdmissionIssueV1::MissingIdentity {
                kind: DerivedObjectKindV1::Boundary,
                field: "orientation_no_claim",
            });
        }
        check_unit(
            boundary.units,
            ir.unit_system,
            DerivedObjectKindV1::Boundary,
            issues,
        );
    }

    for (completed, contact) in ir.contacts.iter().enumerate() {
        checkpoint(cx, "contacts", completed)?;
        if contact.id.is_zero() {
            issues.push(DerivedAdmissionIssueV1::MissingIdentity {
                kind: DerivedObjectKindV1::Contact,
                field: "id",
            });
        }
        let Some(dimension) = chart_dimension(ir, contact.chart) else {
            issues.push(DerivedAdmissionIssueV1::MissingReference {
                kind: DerivedObjectKindV1::Contact,
                field: "chart",
            });
            continue;
        };
        if contact.side_a == contact.side_b {
            issues.push(DerivedAdmissionIssueV1::MissingReference {
                kind: DerivedObjectKindV1::Contact,
                field: "contact_sides",
            });
        }
        for (side, field) in [
            (boundary(ir, contact.side_a), "side_a"),
            (boundary(ir, contact.side_b), "side_b"),
        ] {
            match side {
                None => issues.push(DerivedAdmissionIssueV1::MissingReference {
                    kind: DerivedObjectKindV1::Contact,
                    field,
                }),
                Some(value) if value.chart != contact.chart => {
                    issues.push(DerivedAdmissionIssueV1::CrossChartReference {
                        kind: DerivedObjectKindV1::Contact,
                        field,
                    });
                }
                Some(_) => {}
            }
        }
        check_function(
            contact.gap,
            None,
            ir.category,
            dimension,
            DerivedObjectKindV1::Contact,
            issues,
        );
        check_active_state(contact.state, DerivedObjectKindV1::Contact, issues);
        check_normal_cone(
            contact.normal_cone,
            contact.state,
            DerivedObjectKindV1::Contact,
            issues,
        );
        match contact.law {
            ContactLawV1::Coulomb {
                friction_coefficient,
            } if !friction_coefficient.is_finite() || friction_coefficient.is_sign_negative() => {
                issues.push(DerivedAdmissionIssueV1::InvalidDimension {
                    kind: DerivedObjectKindV1::Contact,
                    field: "friction_coefficient",
                });
            }
            ContactLawV1::SetValued { graph } if graph.is_zero() => {
                issues.push(DerivedAdmissionIssueV1::MissingIdentity {
                    kind: DerivedObjectKindV1::Contact,
                    field: "contact_graph",
                });
            }
            ContactLawV1::Frictionless
            | ContactLawV1::Coulomb { .. }
            | ContactLawV1::SetValued { .. } => {}
        }
        check_unit(
            contact.units,
            ir.unit_system,
            DerivedObjectKindV1::Contact,
            issues,
        );
        check_computability(contact.computability, DerivedObjectKindV1::Contact, issues);
    }

    for (completed, datum) in ir.constitutive_data.iter().enumerate() {
        checkpoint(cx, "constitutive", completed)?;
        if datum.id.is_zero() {
            issues.push(DerivedAdmissionIssueV1::MissingIdentity {
                kind: DerivedObjectKindV1::Constitutive,
                field: "id",
            });
        }
        let Some(dimension) = chart_dimension(ir, datum.chart) else {
            issues.push(DerivedAdmissionIssueV1::MissingReference {
                kind: DerivedObjectKindV1::Constitutive,
                field: "chart",
            });
            continue;
        };
        if datum.state_dimension == 0
            || datum.state_dimension > DERIVED_GEOMETRY_HARD_MAX_DIMENSION_V1
        {
            issues.push(DerivedAdmissionIssueV1::InvalidDimension {
                kind: DerivedObjectKindV1::Constitutive,
                field: "state_dimension",
            });
        }
        check_function(
            datum.law,
            None,
            ir.category,
            dimension,
            DerivedObjectKindV1::Constitutive,
            issues,
        );
        check_unit(
            datum.units,
            ir.unit_system,
            DerivedObjectKindV1::Constitutive,
            issues,
        );
        check_computability(
            datum.computability,
            DerivedObjectKindV1::Constitutive,
            issues,
        );
    }
    Ok(())
}

fn validate_complexes(
    ir: &DerivedGeometryIrV1,
    cx: &Cx<'_>,
    issues: &mut Vec<DerivedAdmissionIssueV1>,
) -> Result<(), DerivedAdmissionReportV1> {
    if has_duplicate(&ir.complexes, |complex| complex.id) {
        issues.push(DerivedAdmissionIssueV1::DuplicateIdentity {
            kind: DerivedObjectKindV1::Complex,
        });
    }
    for (completed, complex) in ir.complexes.iter().enumerate() {
        checkpoint(cx, "complexes", completed)?;
        if complex.id.is_zero() || complex.resolution.id.is_zero() {
            issues.push(DerivedAdmissionIssueV1::MissingIdentity {
                kind: DerivedObjectKindV1::Complex,
                field: "id_or_resolution",
            });
        }
        if !has_chart(ir, complex.chart) {
            issues.push(DerivedAdmissionIssueV1::MissingReference {
                kind: DerivedObjectKindV1::Complex,
                field: "chart",
            });
        }
        if complex.spaces.is_empty() {
            issues.push(DerivedAdmissionIssueV1::EmptyCollection {
                kind: DerivedObjectKindV1::Complex,
            });
            continue;
        }
        if has_duplicate(&complex.spaces, |space| space.degree) {
            issues.push(DerivedAdmissionIssueV1::InvalidComplex {
                field: "duplicate_degree",
            });
        }
        if has_duplicate(&complex.differentials, |map| {
            (map.from_degree, map.to_degree)
        }) {
            issues.push(DerivedAdmissionIssueV1::InvalidComplex {
                field: "duplicate_differential",
            });
        }
        if complex
            .spaces
            .windows(2)
            .any(|spaces| spaces[0].degree.checked_add(1) != Some(spaces[1].degree))
        {
            issues.push(DerivedAdmissionIssueV1::InvalidComplex {
                field: "noncontiguous_degree",
            });
        }
        if complex.differentials.len() != complex.spaces.len().saturating_sub(1) {
            issues.push(DerivedAdmissionIssueV1::InvalidComplex {
                field: "differential_coverage",
            });
        }
        let first_degree = complex.spaces.first().map(|space| space.degree);
        let last_degree = complex.spaces.last().map(|space| space.degree);
        let maximum_dimension = complex
            .spaces
            .iter()
            .map(|space| space.dimension)
            .max()
            .unwrap_or(0);
        if first_degree != Some(complex.resolution.min_degree)
            || last_degree != Some(complex.resolution.max_degree)
            || complex.resolution.min_degree > complex.resolution.max_degree
            || complex.resolution.max_basis_dimension > DERIVED_GEOMETRY_HARD_MAX_DIMENSION_V1
            || complex.resolution.max_basis_dimension < maximum_dimension
        {
            issues.push(DerivedAdmissionIssueV1::InvalidComplex {
                field: "resolution_bounds",
            });
        }
        if (complex.resolution.truncation_order == 0) != complex.resolution.remainder.is_none()
            || complex
                .resolution
                .remainder
                .is_some_and(|witness| witness.is_zero())
        {
            issues.push(DerivedAdmissionIssueV1::InvalidComplex {
                field: "truncation_remainder",
            });
        }
        for space in &complex.spaces {
            if space.dimension > DERIVED_GEOMETRY_HARD_MAX_DIMENSION_V1 || space.quantity.is_zero()
            {
                issues.push(DerivedAdmissionIssueV1::InvalidDimension {
                    kind: DerivedObjectKindV1::Complex,
                    field: "graded_space",
                });
            }
        }
        for differential in &complex.differentials {
            let adjacent = differential
                .from_degree
                .checked_add(1)
                .is_some_and(|next| next == differential.to_degree);
            let source_exists = unique_by_key(&complex.spaces, differential.from_degree, |space| {
                space.degree
            })
            .is_some();
            let target_exists = unique_by_key(&complex.spaces, differential.to_degree, |space| {
                space.degree
            })
            .is_some();
            if !adjacent || !source_exists || !target_exists {
                issues.push(DerivedAdmissionIssueV1::InvalidComplex {
                    field: "differential_degrees",
                });
            }
            if differential.map.is_zero() || differential.square_zero_witness.is_zero() {
                issues.push(DerivedAdmissionIssueV1::MissingIdentity {
                    kind: DerivedObjectKindV1::Complex,
                    field: "differential_map_or_witness",
                });
            }
        }
        check_computability(complex.computability, DerivedObjectKindV1::Complex, issues);
    }
    Ok(())
}

fn validate_model_complex(
    ir: &DerivedGeometryIrV1,
    model: &DerivedLocalModelV1,
    id: DerivedComplexIdV1,
    expected: DerivedComplexRoleV1,
    field: &'static str,
    issues: &mut Vec<DerivedAdmissionIssueV1>,
) {
    match complex(ir, id) {
        Some(value) if value.chart != model.chart => {
            issues.push(DerivedAdmissionIssueV1::CrossChartReference {
                kind: DerivedObjectKindV1::LocalModel,
                field,
            });
        }
        Some(value) if value.role != expected => {
            issues.push(DerivedAdmissionIssueV1::ComplexRoleMismatch { field });
        }
        Some(_) => {}
        None => issues.push(DerivedAdmissionIssueV1::MissingReference {
            kind: DerivedObjectKindV1::LocalModel,
            field,
        }),
    }
}

fn validate_presentation(
    presentation: PresentationScopeV1,
    issues: &mut Vec<DerivedAdmissionIssueV1>,
) {
    let missing = match presentation {
        PresentationScopeV1::Literal { no_claim } => no_claim.is_zero(),
        PresentationScopeV1::FixedResolution {
            resolution,
            witness,
        } => resolution.is_zero() || witness.is_zero(),
        PresentationScopeV1::ExternallyChecked { witness } => witness.is_zero(),
    };
    if missing {
        issues.push(DerivedAdmissionIssueV1::MissingIdentity {
            kind: DerivedObjectKindV1::LocalModel,
            field: "presentation_scope",
        });
    }
}

fn validate_local_models(
    ir: &DerivedGeometryIrV1,
    cx: &Cx<'_>,
    issues: &mut Vec<DerivedAdmissionIssueV1>,
) -> Result<(), DerivedAdmissionReportV1> {
    if has_duplicate(&ir.local_models, |model| model.id) {
        issues.push(DerivedAdmissionIssueV1::DuplicateIdentity {
            kind: DerivedObjectKindV1::LocalModel,
        });
    }
    for (completed, model) in ir.local_models.iter().enumerate() {
        checkpoint(cx, "local-models", completed)?;
        if model.id.is_zero() {
            issues.push(DerivedAdmissionIssueV1::MissingIdentity {
                kind: DerivedObjectKindV1::LocalModel,
                field: "id",
            });
        }
        let Some(chart_dimension) = chart_dimension(ir, model.chart) else {
            issues.push(DerivedAdmissionIssueV1::MissingReference {
                kind: DerivedObjectKindV1::LocalModel,
                field: "chart",
            });
            continue;
        };
        check_locality(
            model.locality,
            &ir.charts,
            DerivedObjectKindV1::LocalModel,
            issues,
        );
        match model.locality {
            LocalityScopeV1::GermAt { chart, .. }
            | LocalityScopeV1::CompactNeighborhood { chart, .. }
                if chart != model.chart =>
            {
                if has_chart(ir, chart) {
                    issues.push(DerivedAdmissionIssueV1::CrossChartReference {
                        kind: DerivedObjectKindV1::LocalModel,
                        field: "locality_chart",
                    });
                }
            }
            LocalityScopeV1::GermAt { .. }
            | LocalityScopeV1::CompactNeighborhood { .. }
            | LocalityScopeV1::GlobalCompact { .. }
            | LocalityScopeV1::GlobalUnbounded
            | LocalityScopeV1::InfiniteDimensional => {}
        }
        validate_presentation(model.presentation, issues);

        for (duplicate, field) in [
            (
                model.equalities.windows(2).any(|ids| ids[0] == ids[1]),
                "duplicate_equality",
            ),
            (
                model
                    .active_inequalities
                    .windows(2)
                    .any(|ids| ids[0] == ids[1]),
                "duplicate_inequality",
            ),
            (
                model.active_contacts.windows(2).any(|ids| ids[0] == ids[1]),
                "duplicate_contact",
            ),
            (
                model
                    .constitutive_data
                    .windows(2)
                    .any(|ids| ids[0] == ids[1]),
                "duplicate_constitutive",
            ),
        ] {
            if duplicate {
                issues.push(DerivedAdmissionIssueV1::InvalidDimension {
                    kind: DerivedObjectKindV1::LocalModel,
                    field,
                });
            }
        }

        let mut equality_codimension = 0u64;
        for equality_id in &model.equalities {
            match unique_by_key(&ir.equalities, *equality_id, |constraint| constraint.id) {
                Some(equality) if equality.chart != model.chart => {
                    issues.push(DerivedAdmissionIssueV1::CrossChartReference {
                        kind: DerivedObjectKindV1::LocalModel,
                        field: "equality",
                    });
                }
                Some(equality) => {
                    equality_codimension += u64::from(equality.codomain_dimension);
                }
                None => issues.push(DerivedAdmissionIssueV1::MissingReference {
                    kind: DerivedObjectKindV1::LocalModel,
                    field: "equality",
                }),
            }
        }
        for inequality_id in &model.active_inequalities {
            match inequality(ir, *inequality_id) {
                Some(value) if value.chart != model.chart => {
                    issues.push(DerivedAdmissionIssueV1::CrossChartReference {
                        kind: DerivedObjectKindV1::LocalModel,
                        field: "active_inequality",
                    });
                }
                Some(value) if matches!(value.state, ActiveSetStateV1::Active { .. }) => {}
                Some(_) => issues.push(DerivedAdmissionIssueV1::ActiveSetMismatch {
                    kind: DerivedObjectKindV1::Inequality,
                }),
                None => issues.push(DerivedAdmissionIssueV1::MissingReference {
                    kind: DerivedObjectKindV1::LocalModel,
                    field: "active_inequality",
                }),
            }
        }
        for contact_id in &model.active_contacts {
            match contact(ir, *contact_id) {
                Some(value) if value.chart != model.chart => {
                    issues.push(DerivedAdmissionIssueV1::CrossChartReference {
                        kind: DerivedObjectKindV1::LocalModel,
                        field: "active_contact",
                    });
                }
                Some(value) if matches!(value.state, ActiveSetStateV1::Active { .. }) => {}
                Some(_) => issues.push(DerivedAdmissionIssueV1::ActiveSetMismatch {
                    kind: DerivedObjectKindV1::Contact,
                }),
                None => issues.push(DerivedAdmissionIssueV1::MissingReference {
                    kind: DerivedObjectKindV1::LocalModel,
                    field: "active_contact",
                }),
            }
        }
        for datum_id in &model.constitutive_data {
            match constitutive(ir, *datum_id) {
                Some(datum) if datum.chart == model.chart => {}
                Some(_) => issues.push(DerivedAdmissionIssueV1::CrossChartReference {
                    kind: DerivedObjectKindV1::LocalModel,
                    field: "constitutive_data",
                }),
                None => issues.push(DerivedAdmissionIssueV1::MissingReference {
                    kind: DerivedObjectKindV1::LocalModel,
                    field: "constitutive_data",
                }),
            }
        }

        validate_model_complex(
            ir,
            model,
            model.tangent_complex,
            DerivedComplexRoleV1::Tangent,
            "tangent_complex",
            issues,
        );
        validate_model_complex(
            ir,
            model,
            model.cotangent_complex,
            DerivedComplexRoleV1::Cotangent,
            "cotangent_complex",
            issues,
        );
        validate_model_complex(
            ir,
            model,
            model.deformation_complex,
            DerivedComplexRoleV1::DeformationObstruction,
            "deformation_complex",
            issues,
        );

        let active_codimension = model
            .active_inequalities
            .len()
            .saturating_add(model.active_contacts.len()) as u64;
        if model.class == DerivedLocalModelClassV1::RegularCompleteIntersection {
            let expected = i64::from(chart_dimension)
                - i64::try_from(equality_codimension).unwrap_or(i64::MAX)
                - i64::try_from(active_codimension).unwrap_or(i64::MAX);
            if expected != i64::from(model.virtual_dimension) || expected < 0 {
                issues.push(DerivedAdmissionIssueV1::InvalidDimension {
                    kind: DerivedObjectKindV1::LocalModel,
                    field: "regular_virtual_dimension",
                });
            }
        }
        if model.virtual_dimension.unsigned_abs() > DERIVED_GEOMETRY_HARD_MAX_DIMENSION_V1 {
            issues.push(DerivedAdmissionIssueV1::InvalidDimension {
                kind: DerivedObjectKindV1::LocalModel,
                field: "virtual_dimension",
            });
        }
        match model.class {
            DerivedLocalModelClassV1::RedundantPresentation
            | DerivedLocalModelClassV1::Cusp
            | DerivedLocalModelClassV1::Node
                if model.equalities.is_empty() =>
            {
                issues.push(DerivedAdmissionIssueV1::EmptyCollection {
                    kind: DerivedObjectKindV1::Equality,
                });
            }
            DerivedLocalModelClassV1::ContactCorner
                if model.active_inequalities.is_empty() && model.active_contacts.is_empty() =>
            {
                issues.push(DerivedAdmissionIssueV1::EmptyCollection {
                    kind: DerivedObjectKindV1::Contact,
                });
            }
            DerivedLocalModelClassV1::RegularCompleteIntersection
            | DerivedLocalModelClassV1::RedundantPresentation
            | DerivedLocalModelClassV1::Cusp
            | DerivedLocalModelClassV1::Node
            | DerivedLocalModelClassV1::Boundary
            | DerivedLocalModelClassV1::ContactCorner
            | DerivedLocalModelClassV1::GeneralFiniteDerived => {}
        }
    }
    Ok(())
}

fn incidence(
    ir: &DerivedGeometryIrV1,
    lower: StratumIdV1,
    upper: StratumIdV1,
) -> Option<&StratumIncidenceV1> {
    unique_by_key(&ir.stratification.incidences, (lower, upper), |edge| {
        (edge.lower, edge.upper)
    })
}

fn validate_stratification(
    ir: &DerivedGeometryIrV1,
    cx: &Cx<'_>,
    issues: &mut Vec<DerivedAdmissionIssueV1>,
) -> Result<(), DerivedAdmissionReportV1> {
    if has_duplicate(&ir.stratification.strata, |stratum| stratum.id) {
        issues.push(DerivedAdmissionIssueV1::DuplicateIdentity {
            kind: DerivedObjectKindV1::Stratum,
        });
    }
    if has_duplicate(&ir.stratification.incidences, |edge| {
        (edge.lower, edge.upper)
    }) {
        issues.push(DerivedAdmissionIssueV1::DuplicateIdentity {
            kind: DerivedObjectKindV1::Incidence,
        });
    }
    if has_duplicate(&ir.stratification.local_links, |link| link.id) {
        issues.push(DerivedAdmissionIssueV1::DuplicateIdentity {
            kind: DerivedObjectKindV1::LocalLink,
        });
    }
    match ir.stratification.class {
        StratificationClassV1::FiniteIncidence => {}
        StratificationClassV1::WhitneyA { witness }
        | StratificationClassV1::WhitneyB { witness }
        | StratificationClassV1::Thom { witness }
            if witness.is_zero() =>
        {
            issues.push(DerivedAdmissionIssueV1::MissingIdentity {
                kind: DerivedObjectKindV1::Global,
                field: "stratification_theorem_witness",
            });
        }
        StratificationClassV1::WhitneyA { .. }
        | StratificationClassV1::WhitneyB { .. }
        | StratificationClassV1::Thom { .. } => {}
    }

    for (completed, stratum) in ir.stratification.strata.iter().enumerate() {
        checkpoint(cx, "strata", completed)?;
        if stratum.id.is_zero() {
            issues.push(DerivedAdmissionIssueV1::MissingIdentity {
                kind: DerivedObjectKindV1::Stratum,
                field: "id",
            });
        }
        let chart_dimension = chart_dimension(ir, stratum.chart);
        if chart_dimension.is_none() {
            issues.push(DerivedAdmissionIssueV1::MissingReference {
                kind: DerivedObjectKindV1::Stratum,
                field: "chart",
            });
        }
        if chart_dimension.is_some_and(|dimension| stratum.dimension > dimension) {
            issues.push(DerivedAdmissionIssueV1::InvalidDimension {
                kind: DerivedObjectKindV1::Stratum,
                field: "dimension",
            });
        }
        let model = local_model(ir, stratum.local_model);
        match model {
            None => issues.push(DerivedAdmissionIssueV1::MissingReference {
                kind: DerivedObjectKindV1::Stratum,
                field: "local_model",
            }),
            Some(value) if value.chart != stratum.chart => {
                issues.push(DerivedAdmissionIssueV1::CrossChartReference {
                    kind: DerivedObjectKindV1::Stratum,
                    field: "local_model",
                });
            }
            Some(_) => {}
        }
        if model.is_some_and(|value| {
            value.active_inequalities != stratum.active_inequalities
                || value.active_contacts != stratum.active_contacts
        }) {
            issues.push(DerivedAdmissionIssueV1::InvalidStratification {
                field: "active_set",
            });
        }
        if let Some(boundary_id) = stratum.relative_boundary {
            match boundary(ir, boundary_id) {
                None => issues.push(DerivedAdmissionIssueV1::MissingReference {
                    kind: DerivedObjectKindV1::Stratum,
                    field: "relative_boundary",
                }),
                Some(value) if value.boundary != stratum.id => {
                    issues.push(DerivedAdmissionIssueV1::MissingReference {
                        kind: DerivedObjectKindV1::Stratum,
                        field: "relative_boundary",
                    });
                }
                Some(value) if value.chart != stratum.chart => {
                    issues.push(DerivedAdmissionIssueV1::CrossChartReference {
                        kind: DerivedObjectKindV1::Stratum,
                        field: "relative_boundary",
                    });
                }
                Some(_) => {}
            }
        }
        check_regularity(stratum.regularity, DerivedObjectKindV1::Stratum, issues);
        check_compactness(stratum.compactness, DerivedObjectKindV1::Stratum, issues);
    }

    for (completed, edge) in ir.stratification.incidences.iter().enumerate() {
        checkpoint(cx, "incidences", completed)?;
        let lower = stratum(ir, edge.lower);
        let upper = stratum(ir, edge.upper);
        if edge.witness.is_zero() {
            issues.push(DerivedAdmissionIssueV1::MissingIdentity {
                kind: DerivedObjectKindV1::Incidence,
                field: "witness",
            });
        }
        if lower.is_none() {
            issues.push(DerivedAdmissionIssueV1::MissingReference {
                kind: DerivedObjectKindV1::Incidence,
                field: "lower",
            });
        }
        if upper.is_none() {
            issues.push(DerivedAdmissionIssueV1::MissingReference {
                kind: DerivedObjectKindV1::Incidence,
                field: "upper",
            });
        }
        if let (Some(lower), Some(upper)) = (lower, upper) {
            if lower.chart != upper.chart {
                issues.push(DerivedAdmissionIssueV1::CrossChartReference {
                    kind: DerivedObjectKindV1::Incidence,
                    field: "lower_upper",
                });
            }
            if edge.lower == edge.upper
                || lower.dimension >= upper.dimension
                || upper.dimension - lower.dimension != edge.codimension
                || edge.codimension == 0
            {
                issues.push(DerivedAdmissionIssueV1::InvalidStratification {
                    field: "incidence_dimension",
                });
            }
        }
    }

    for boundary in &ir.boundaries {
        if incidence(ir, boundary.boundary, boundary.parent).is_none() {
            issues.push(DerivedAdmissionIssueV1::InvalidStratification {
                field: "boundary_incidence",
            });
        }
    }

    for (completed, link) in ir.stratification.local_links.iter().enumerate() {
        checkpoint(cx, "local-links", completed)?;
        if link.id.is_zero() || link.compactness_witness.is_zero() {
            issues.push(DerivedAdmissionIssueV1::MissingIdentity {
                kind: DerivedObjectKindV1::LocalLink,
                field: "id_or_compactness",
            });
        }
        let edge = incidence(ir, link.stratum, link.ambient_stratum);
        if edge.is_none() {
            issues.push(DerivedAdmissionIssueV1::MissingReference {
                kind: DerivedObjectKindV1::LocalLink,
                field: "incidence",
            });
        } else if edge.is_some_and(|value| link.dimension.checked_add(1) != Some(value.codimension))
        {
            issues.push(DerivedAdmissionIssueV1::InvalidLocalLink { field: "dimension" });
        }
        match link.topology {
            LocalLinkTopologyV1::FiniteComplex {
                resolution,
                witness,
            } if resolution.is_zero() || witness.is_zero() => {
                issues.push(DerivedAdmissionIssueV1::MissingIdentity {
                    kind: DerivedObjectKindV1::LocalLink,
                    field: "topology",
                });
            }
            LocalLinkTopologyV1::Unknown { no_claim } if no_claim.is_zero() => {
                issues.push(DerivedAdmissionIssueV1::MissingIdentity {
                    kind: DerivedObjectKindV1::LocalLink,
                    field: "topology_no_claim",
                });
            }
            LocalLinkTopologyV1::FiniteComplex { .. } | LocalLinkTopologyV1::Unknown { .. } => {}
        }
    }
    Ok(())
}

fn validate_proof_state(ir: &DerivedGeometryIrV1, issues: &mut Vec<DerivedAdmissionIssueV1>) {
    match ir.proof_state {
        DerivedProofStateV1::StructuralNoClaim { no_claim } => {
            if no_claim.is_zero() {
                issues.push(DerivedAdmissionIssueV1::InvalidProofState);
            }
        }
        DerivedProofStateV1::ExternallyChecked {
            theorem,
            checker,
            receipt,
            scope,
        } => {
            let scope_exists = match scope {
                DerivedProofScopeV1::Object => true,
                DerivedProofScopeV1::LocalModel(id) => local_model(ir, id).is_some(),
                DerivedProofScopeV1::Stratification(id) => id == ir.stratification.id,
            };
            if theorem.is_zero() || checker.is_zero() || receipt.is_zero() || !scope_exists {
                issues.push(DerivedAdmissionIssueV1::InvalidProofState);
            }
        }
    }
}

#[derive(Default)]
struct Wire(Vec<u8>);

impl Wire {
    fn tag(&mut self, value: u8) {
        self.0.push(value);
    }

    fn u16(&mut self, value: u16) {
        self.0.extend_from_slice(&value.to_le_bytes());
    }

    fn i16(&mut self, value: i16) {
        self.0.extend_from_slice(&value.to_le_bytes());
    }

    fn u32(&mut self, value: u32) {
        self.0.extend_from_slice(&value.to_le_bytes());
    }

    fn i32(&mut self, value: i32) {
        self.0.extend_from_slice(&value.to_le_bytes());
    }

    fn u64(&mut self, value: u64) {
        self.0.extend_from_slice(&value.to_le_bytes());
    }

    fn f64(&mut self, value: f64) {
        self.u64(value.to_bits());
    }

    fn id<I: DigestBytes>(&mut self, value: I) {
        self.0.extend_from_slice(value.digest_bytes());
    }

    fn optional_id<I: DigestBytes + Copy>(&mut self, value: Option<I>) {
        match value {
            Some(value) => {
                self.tag(1);
                self.id(value);
            }
            None => self.tag(0),
        }
    }

    fn ids<I: DigestBytes + Copy>(&mut self, values: &[I]) {
        self.u64(values.len() as u64);
        for value in values {
            self.id(*value);
        }
    }

    fn frames(&mut self, values: &[Vec<u8>]) {
        self.u64(values.len() as u64);
        for value in values {
            self.u64(value.len() as u64);
            self.0.extend_from_slice(value);
        }
    }

    fn finish(self) -> Vec<u8> {
        self.0
    }
}

fn encode_category(ir: &DerivedGeometryIrV1) -> Vec<u8> {
    let mut wire = Wire::default();
    wire.tag(ir.category.tag());
    if let GeometricCategoryV1::Subanalytic { construction } = ir.category {
        wire.id(construction);
    }
    wire.tag(ir.coefficients.tag());
    if let CoefficientSystemV1::DyadicIntervalReal { precision_bits } = ir.coefficients {
        wire.u16(precision_bits);
    }
    wire.finish()
}

fn encode_locality(wire: &mut Wire, locality: LocalityScopeV1) {
    match locality {
        LocalityScopeV1::GermAt { chart, point } => {
            wire.tag(0);
            wire.id(chart);
            wire.id(point);
        }
        LocalityScopeV1::CompactNeighborhood { chart, witness } => {
            wire.tag(1);
            wire.id(chart);
            wire.id(witness);
        }
        LocalityScopeV1::GlobalCompact { witness } => {
            wire.tag(2);
            wire.id(witness);
        }
        LocalityScopeV1::GlobalUnbounded => wire.tag(3),
        LocalityScopeV1::InfiniteDimensional => wire.tag(4),
    }
}

fn encode_compactness(wire: &mut Wire, compactness: CompactnessV1) {
    match compactness {
        CompactnessV1::Proved { witness } => {
            wire.tag(0);
            wire.id(witness);
        }
        CompactnessV1::RelativelyCompact { witness } => {
            wire.tag(1);
            wire.id(witness);
        }
        CompactnessV1::Assumed { no_claim } => {
            wire.tag(2);
            wire.id(no_claim);
        }
        CompactnessV1::Unbounded => wire.tag(3),
        CompactnessV1::Unknown => wire.tag(4),
    }
}

fn encode_regularity(wire: &mut Wire, regularity: RegularityClassV1) {
    match regularity {
        RegularityClassV1::Polynomial => wire.tag(0),
        RegularityClassV1::Analytic => wire.tag(1),
        RegularityClassV1::Differentiable { order } => {
            wire.tag(2);
            wire.u16(order);
        }
        RegularityClassV1::Unknown => wire.tag(3),
    }
}

fn encode_computability(wire: &mut Wire, computability: FiniteComputabilityV1) {
    match computability {
        FiniteComputabilityV1::ExactFinite { kernel } => {
            wire.tag(0);
            wire.id(kernel);
        }
        FiniteComputabilityV1::IntervalFinite { enclosure } => {
            wire.tag(1);
            wire.id(enclosure);
        }
        FiniteComputabilityV1::TruncatedFinite {
            resolution,
            remainder,
        } => {
            wire.tag(2);
            wire.id(resolution);
            wire.id(remainder);
        }
        FiniteComputabilityV1::ExternalOpaque { no_claim } => {
            wire.tag(3);
            wire.id(no_claim);
        }
        FiniteComputabilityV1::Infinite => wire.tag(4),
    }
}

fn encode_function(wire: &mut Wire, function: LocalFunctionEncodingV1) {
    match function {
        LocalFunctionEncodingV1::Polynomial {
            polynomial,
            variables,
            degree,
        } => {
            wire.tag(0);
            wire.id(polynomial);
            wire.u32(variables);
            wire.u32(degree);
        }
        LocalFunctionEncodingV1::RestrictedAnalytic {
            program,
            primitives,
            derivative_order,
        } => {
            wire.tag(1);
            wire.id(program);
            wire.id(primitives);
            wire.u16(derivative_order);
        }
        LocalFunctionEncodingV1::OpaqueExternal { no_claim } => {
            wire.tag(2);
            wire.id(no_claim);
        }
    }
}

fn encode_unit(wire: &mut Wire, units: UnitBindingV1) {
    wire.id(units.system);
    wire.id(units.quantity);
    wire.f64(units.scale_to_canonical);
}

fn encode_global_context(ir: &DerivedGeometryIrV1) -> Vec<u8> {
    let mut wire = Wire::default();
    wire.id(ir.frame);
    wire.id(ir.unit_system);
    encode_locality(&mut wire, ir.locality);
    encode_compactness(&mut wire, ir.compactness);
    wire.finish()
}

fn encode_chart(chart: &ConfigurationChartV1) -> Vec<u8> {
    let mut wire = Wire::default();
    wire.id(chart.id);
    wire.tag(match chart.class {
        ConfigurationChartClassV1::SmoothManifold => 0,
        ConfigurationChartClassV1::Semialgebraic => 1,
        ConfigurationChartClassV1::Algebraic => 2,
        ConfigurationChartClassV1::RestrictedAnalytic => 3,
        ConfigurationChartClassV1::Stratified => 4,
    });
    wire.u32(chart.coordinate_dimension);
    wire.u32(chart.ambient_dimension);
    wire.id(chart.frame);
    encode_unit(&mut wire, chart.coordinates);
    encode_locality(&mut wire, chart.locality);
    encode_compactness(&mut wire, chart.compactness);
    encode_regularity(&mut wire, chart.regularity);
    encode_computability(&mut wire, chart.computability);
    wire.finish()
}

fn encode_equality(constraint: &EqualityConstraintGermV1) -> Vec<u8> {
    let mut wire = Wire::default();
    wire.id(constraint.id);
    wire.id(constraint.chart);
    wire.u32(constraint.codomain_dimension);
    encode_function(&mut wire, constraint.equation);
    encode_regularity(&mut wire, constraint.regularity);
    encode_unit(&mut wire, constraint.units);
    encode_computability(&mut wire, constraint.computability);
    wire.finish()
}

fn encode_active_state(wire: &mut Wire, state: ActiveSetStateV1) {
    match state {
        ActiveSetStateV1::Inactive { witness } => {
            wire.tag(0);
            wire.id(witness);
        }
        ActiveSetStateV1::Active { witness } => {
            wire.tag(1);
            wire.id(witness);
        }
        ActiveSetStateV1::Candidate { no_claim } => {
            wire.tag(2);
            wire.id(no_claim);
        }
    }
}

fn encode_normal_cone(wire: &mut Wire, cone: NormalConeClassV1) {
    match cone {
        NormalConeClassV1::Ray => wire.tag(0),
        NormalConeClassV1::Polyhedral { generators } => {
            wire.tag(1);
            wire.u32(generators);
        }
        NormalConeClassV1::SmoothDual => wire.tag(2),
        NormalConeClassV1::Unknown => wire.tag(3),
    }
}

fn encode_inequality(constraint: &InequalityConstraintGermV1) -> Vec<u8> {
    let mut wire = Wire::default();
    wire.id(constraint.id);
    wire.id(constraint.chart);
    wire.tag(match constraint.sense {
        InequalitySenseV1::NonNegative => 0,
        InequalitySenseV1::NonPositive => 1,
    });
    encode_function(&mut wire, constraint.function);
    encode_active_state(&mut wire, constraint.state);
    encode_normal_cone(&mut wire, constraint.normal_cone);
    encode_unit(&mut wire, constraint.units);
    encode_computability(&mut wire, constraint.computability);
    wire.finish()
}

fn encode_boundary(boundary: &RelativeBoundaryV1) -> Vec<u8> {
    let mut wire = Wire::default();
    wire.id(boundary.id);
    wire.id(boundary.chart);
    wire.id(boundary.parent);
    wire.id(boundary.boundary);
    match boundary.orientation {
        BoundaryOrientationV1::Outward => wire.tag(0),
        BoundaryOrientationV1::Inward => wire.tag(1),
        BoundaryOrientationV1::Unoriented { no_claim } => {
            wire.tag(2);
            wire.id(no_claim);
        }
    }
    wire.id(boundary.witness);
    encode_unit(&mut wire, boundary.units);
    wire.finish()
}

fn encode_contact(contact: &ContactConstraintV1) -> Vec<u8> {
    let mut wire = Wire::default();
    wire.id(contact.id);
    wire.id(contact.chart);
    wire.id(contact.side_a);
    wire.id(contact.side_b);
    encode_function(&mut wire, contact.gap);
    encode_active_state(&mut wire, contact.state);
    encode_normal_cone(&mut wire, contact.normal_cone);
    match contact.law {
        ContactLawV1::Frictionless => wire.tag(0),
        ContactLawV1::Coulomb {
            friction_coefficient,
        } => {
            wire.tag(1);
            wire.f64(friction_coefficient);
        }
        ContactLawV1::SetValued { graph } => {
            wire.tag(2);
            wire.id(graph);
        }
    }
    encode_unit(&mut wire, contact.units);
    encode_computability(&mut wire, contact.computability);
    wire.finish()
}

fn encode_constitutive(datum: &ConstitutiveDatumV1) -> Vec<u8> {
    let mut wire = Wire::default();
    wire.id(datum.id);
    wire.id(datum.chart);
    wire.tag(match datum.role {
        ConstitutiveRoleV1::Energy => 0,
        ConstitutiveRoleV1::Dissipation => 1,
        ConstitutiveRoleV1::GeneralRelation => 2,
    });
    wire.u32(datum.state_dimension);
    encode_function(&mut wire, datum.law);
    encode_unit(&mut wire, datum.units);
    encode_computability(&mut wire, datum.computability);
    wire.finish()
}

fn encode_complex(complex: &FiniteDerivedComplexV1) -> Vec<u8> {
    let mut wire = Wire::default();
    wire.id(complex.id);
    wire.id(complex.chart);
    wire.tag(match complex.role {
        DerivedComplexRoleV1::Tangent => 0,
        DerivedComplexRoleV1::Cotangent => 1,
        DerivedComplexRoleV1::DeformationObstruction => 2,
    });
    wire.u64(complex.spaces.len() as u64);
    for space in &complex.spaces {
        wire.i16(space.degree);
        wire.u32(space.dimension);
        wire.id(space.quantity);
    }
    wire.u64(complex.differentials.len() as u64);
    for differential in &complex.differentials {
        wire.i16(differential.from_degree);
        wire.i16(differential.to_degree);
        wire.id(differential.map);
        wire.id(differential.square_zero_witness);
    }
    wire.id(complex.resolution.id);
    wire.i16(complex.resolution.min_degree);
    wire.i16(complex.resolution.max_degree);
    wire.u32(complex.resolution.max_basis_dimension);
    wire.u32(complex.resolution.truncation_order);
    wire.optional_id(complex.resolution.remainder);
    encode_computability(&mut wire, complex.computability);
    wire.finish()
}

fn encode_presentation(wire: &mut Wire, presentation: PresentationScopeV1) {
    match presentation {
        PresentationScopeV1::Literal { no_claim } => {
            wire.tag(0);
            wire.id(no_claim);
        }
        PresentationScopeV1::FixedResolution {
            resolution,
            witness,
        } => {
            wire.tag(1);
            wire.id(resolution);
            wire.id(witness);
        }
        PresentationScopeV1::ExternallyChecked { witness } => {
            wire.tag(2);
            wire.id(witness);
        }
    }
}

fn encode_local_model(model: &DerivedLocalModelV1) -> Vec<u8> {
    let mut wire = Wire::default();
    wire.id(model.id);
    wire.id(model.chart);
    wire.tag(match model.class {
        DerivedLocalModelClassV1::RegularCompleteIntersection => 0,
        DerivedLocalModelClassV1::RedundantPresentation => 1,
        DerivedLocalModelClassV1::Cusp => 2,
        DerivedLocalModelClassV1::Node => 3,
        DerivedLocalModelClassV1::Boundary => 4,
        DerivedLocalModelClassV1::ContactCorner => 5,
        DerivedLocalModelClassV1::GeneralFiniteDerived => 6,
    });
    wire.ids(&model.equalities);
    wire.ids(&model.active_inequalities);
    wire.ids(&model.active_contacts);
    wire.ids(&model.constitutive_data);
    wire.id(model.tangent_complex);
    wire.id(model.cotangent_complex);
    wire.id(model.deformation_complex);
    wire.i32(model.virtual_dimension);
    encode_locality(&mut wire, model.locality);
    encode_presentation(&mut wire, model.presentation);
    wire.finish()
}

fn encode_stratum(stratum: &StratumSpecV1) -> Vec<u8> {
    let mut wire = Wire::default();
    wire.id(stratum.id);
    wire.id(stratum.chart);
    wire.id(stratum.local_model);
    wire.u32(stratum.dimension);
    wire.ids(&stratum.active_inequalities);
    wire.ids(&stratum.active_contacts);
    wire.optional_id(stratum.relative_boundary);
    encode_regularity(&mut wire, stratum.regularity);
    encode_compactness(&mut wire, stratum.compactness);
    wire.finish()
}

fn encode_incidence(edge: &StratumIncidenceV1) -> Vec<u8> {
    let mut wire = Wire::default();
    wire.id(edge.lower);
    wire.id(edge.upper);
    wire.u32(edge.codimension);
    wire.id(edge.witness);
    wire.finish()
}

fn encode_local_link(link: &LocalLinkV1) -> Vec<u8> {
    let mut wire = Wire::default();
    wire.id(link.id);
    wire.id(link.stratum);
    wire.id(link.ambient_stratum);
    wire.u32(link.dimension);
    wire.id(link.compactness_witness);
    match link.topology {
        LocalLinkTopologyV1::FiniteComplex {
            resolution,
            witness,
        } => {
            wire.tag(0);
            wire.id(resolution);
            wire.id(witness);
        }
        LocalLinkTopologyV1::Unknown { no_claim } => {
            wire.tag(1);
            wire.id(no_claim);
        }
    }
    wire.finish()
}

fn encode_stratification(
    stratification: &StratificationV1,
    cx: &Cx<'_>,
) -> Result<Vec<u8>, CanonicalError> {
    let mut wire = Wire::default();
    wire.id(stratification.id);
    match stratification.class {
        StratificationClassV1::FiniteIncidence => wire.tag(0),
        StratificationClassV1::WhitneyA { witness } => {
            wire.tag(1);
            wire.id(witness);
        }
        StratificationClassV1::WhitneyB { witness } => {
            wire.tag(2);
            wire.id(witness);
        }
        StratificationClassV1::Thom { witness } => {
            wire.tag(3);
            wire.id(witness);
        }
    }
    let strata = encode_collection(&stratification.strata, encode_stratum, cx)?;
    let incidences = encode_collection(&stratification.incidences, encode_incidence, cx)?;
    let links = encode_collection(&stratification.local_links, encode_local_link, cx)?;
    wire.frames(&strata);
    wire.frames(&incidences);
    wire.frames(&links);
    Ok(wire.finish())
}

fn encode_proof_state(proof: DerivedProofStateV1) -> Vec<u8> {
    let mut wire = Wire::default();
    match proof {
        DerivedProofStateV1::StructuralNoClaim { no_claim } => {
            wire.tag(0);
            wire.id(no_claim);
        }
        DerivedProofStateV1::ExternallyChecked {
            theorem,
            checker,
            receipt,
            scope,
        } => {
            wire.tag(1);
            wire.id(theorem);
            wire.id(checker);
            wire.id(receipt);
            match scope {
                DerivedProofScopeV1::Object => wire.tag(0),
                DerivedProofScopeV1::LocalModel(id) => {
                    wire.tag(1);
                    wire.id(id);
                }
                DerivedProofScopeV1::Stratification(id) => {
                    wire.tag(2);
                    wire.id(id);
                }
            }
        }
    }
    wire.finish()
}

fn encode_collection<T>(
    values: &[T],
    encode: fn(&T) -> Vec<u8>,
    cx: &Cx<'_>,
) -> Result<Vec<Vec<u8>>, CanonicalError> {
    let mut encoded = Vec::with_capacity(values.len());
    for value in values {
        if cx.checkpoint().is_err() {
            return Err(CanonicalError::Cancelled { absorbed_bytes: 0 });
        }
        encoded.push(encode(value));
    }
    Ok(encoded)
}

fn derived_geometry_receipt(
    ir: &DerivedGeometryIrV1,
    budget: DerivedAdmissionBudgetV1,
    cx: &Cx<'_>,
) -> Result<IdentityReceipt<DerivedGeometryIdV1>, CanonicalError> {
    let category = encode_category(ir);
    let global_context = encode_global_context(ir);
    let charts = encode_collection(&ir.charts, encode_chart, cx)?;
    let equalities = encode_collection(&ir.equalities, encode_equality, cx)?;
    let inequalities = encode_collection(&ir.inequalities, encode_inequality, cx)?;
    let boundaries = encode_collection(&ir.boundaries, encode_boundary, cx)?;
    let contacts = encode_collection(&ir.contacts, encode_contact, cx)?;
    let constitutive = encode_collection(&ir.constitutive_data, encode_constitutive, cx)?;
    let complexes = encode_collection(&ir.complexes, encode_complex, cx)?;
    let local_models = encode_collection(&ir.local_models, encode_local_model, cx)?;
    if cx.checkpoint().is_err() {
        return Err(CanonicalError::Cancelled { absorbed_bytes: 0 });
    }
    let stratification = encode_stratification(&ir.stratification, cx)?;
    let proof_state = encode_proof_state(ir.proof_state);

    CanonicalEncoder::<DerivedGeometryIdV1, _>::new(budget.canonical_limits(), || {
        cx.checkpoint().is_err()
    })?
    .bytes(Field::new(0, "subject"), ir.subject.as_bytes())?
    .bytes(Field::new(1, "model-version"), ir.model_version.as_bytes())?
    .bytes(Field::new(2, "category"), &category)?
    .bytes(Field::new(3, "global-context"), &global_context)?
    .canonical_set(
        Field::new(4, "charts"),
        charts.len() as u64,
        charts.iter().map(Vec::as_slice),
    )?
    .canonical_set(
        Field::new(5, "equalities"),
        equalities.len() as u64,
        equalities.iter().map(Vec::as_slice),
    )?
    .canonical_set(
        Field::new(6, "inequalities"),
        inequalities.len() as u64,
        inequalities.iter().map(Vec::as_slice),
    )?
    .canonical_set(
        Field::new(7, "boundaries"),
        boundaries.len() as u64,
        boundaries.iter().map(Vec::as_slice),
    )?
    .canonical_set(
        Field::new(8, "contacts"),
        contacts.len() as u64,
        contacts.iter().map(Vec::as_slice),
    )?
    .canonical_set(
        Field::new(9, "constitutive-data"),
        constitutive.len() as u64,
        constitutive.iter().map(Vec::as_slice),
    )?
    .canonical_set(
        Field::new(10, "complexes"),
        complexes.len() as u64,
        complexes.iter().map(Vec::as_slice),
    )?
    .canonical_set(
        Field::new(11, "local-models"),
        local_models.len() as u64,
        local_models.iter().map(Vec::as_slice),
    )?
    .bytes(Field::new(12, "stratification"), &stratification)?
    .bytes(Field::new(13, "proof-state"), &proof_state)?
    .finish()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn empty_ir() -> DerivedGeometryIrV1 {
        let zero = [0; 32];
        DerivedGeometryIrV1 {
            schema_version: DERIVED_GEOMETRY_SCHEMA_VERSION_V1,
            subject: DerivedSubjectIdV1::from_bytes(zero),
            model_version: DerivedModelVersionIdV1::from_bytes(zero),
            category: GeometricCategoryV1::Semialgebraic,
            coefficients: CoefficientSystemV1::RationalReal,
            frame: DerivedFrameIdV1::from_bytes(zero),
            unit_system: DerivedUnitSystemIdV1::from_bytes(zero),
            locality: LocalityScopeV1::GlobalUnbounded,
            compactness: CompactnessV1::Unknown,
            charts: Vec::new(),
            equalities: Vec::new(),
            inequalities: Vec::new(),
            boundaries: Vec::new(),
            contacts: Vec::new(),
            constitutive_data: Vec::new(),
            complexes: Vec::new(),
            local_models: Vec::new(),
            stratification: StratificationV1 {
                id: StratificationIdV1::from_bytes(zero),
                class: StratificationClassV1::FiniteIncidence,
                strata: Vec::new(),
                incidences: Vec::new(),
                local_links: Vec::new(),
            },
            proof_state: DerivedProofStateV1::StructuralNoClaim {
                no_claim: DerivedNoClaimIdV1::from_bytes(zero),
            },
        }
    }

    #[test]
    fn canonicalization_propagates_a_mid_sequence_cancellation() {
        let mut ir = empty_ir();
        let mut observed = Vec::new();
        let report = canonicalize_ir_with(&mut ir, |completed| {
            observed.push(completed);
            if completed == 3 {
                Err(DerivedAdmissionReportV1::one(
                    DerivedAdmissionIssueV1::Cancelled {
                        stage: "canonicalize",
                        completed,
                    },
                ))
            } else {
                Ok(())
            }
        })
        .expect_err("injected canonicalization cancellation must propagate");

        assert_eq!(observed, vec![0, 1, 2, 3]);
        assert_eq!(
            report.issues(),
            &[DerivedAdmissionIssueV1::Cancelled {
                stage: "canonicalize",
                completed: 3,
            }]
        );
    }
}
