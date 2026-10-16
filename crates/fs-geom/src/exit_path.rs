//! Exit-path approximation and constructibility theorem statements (RD.X1).
//!
//! An incidence poset is not silently identified with an exit-path category.
//! This module records the hypotheses under which poset, groupoid, simplicial,
//! finite higher-truncation, or full-higher statements are even well formed.
//! It derives a complete fallback lattice and retains adversarial examples.
//! Structural admission and identity are not theorem proof.

use core::fmt;

use fs_blake3::identity::{
    CanonicalEncoder, CanonicalError, CanonicalLimits, CanonicalSchema, EvidenceNodeId, Field,
    FieldSpec, IdentityReceipt, StrongIdentity, WireType,
};
use fs_exec::Cx;

use crate::derived::{
    AdmittedDerivedGeometryV1, CoefficientSystemV1, DerivedFrameIdV1, DerivedGeometryIdV1,
    DerivedModelVersionIdV1, DerivedUnitSystemIdV1, StratificationIdV1,
};

/// Current RD.X1 statement schema.
pub const EXIT_PATH_SCHEMA_VERSION_V1: u32 = 1;
/// Largest finite homotopy truncation admitted by v1.
pub const MAX_EXIT_PATH_TRUNCATION_V1: u8 = 16;
/// Hard cap on preregistered falsifiers.
pub const MAX_EXIT_PATH_FALSIFIERS_V1: usize = 128;

const EXIT_PATH_IDENTITY_LIMITS: CanonicalLimits =
    CanonicalLimits::new(1 << 20, 1 << 20, 16, 4096, 4096);

trait DigestBytes {
    fn digest_bytes(&self) -> &[u8; 32];

    fn is_zero(&self) -> bool {
        self.digest_bytes().iter().all(|byte| *byte == 0)
    }
}

macro_rules! opaque_exit_id {
    ($(#[$meta:meta])* $name:ident) => {
        $(#[$meta])*
        #[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
        pub struct $name([u8; 32]);

        impl $name {
            /// Construct from exact typed digest bytes.
            #[must_use]
            pub const fn from_bytes(bytes: [u8; 32]) -> Self {
                Self(bytes)
            }

            /// Exact typed digest bytes.
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

opaque_exit_id!(
    /// Exit-path convention bundle.
    ExitPathConventionIdV1
);
opaque_exit_id!(
    /// Exact stratified-path equivalence relation.
    StratifiedPathEquivalenceIdV1
);
opaque_exit_id!(
    /// Link catalog or link-homotopy artifact.
    ExitLinkCatalogIdV1
);
opaque_exit_id!(
    /// Stratum fundamental-group/groupoid catalog.
    StratumGroupoidCatalogIdV1
);
opaque_exit_id!(
    /// Constructible local-system catalog.
    LocalSystemCatalogIdV1
);
opaque_exit_id!(
    /// Retained simplicial or higher-coherence artifact.
    HigherCoherenceArtifactIdV1
);
opaque_exit_id!(
    /// Deterministic stratification refinement identity.
    StratifiedRefinementIdV1
);
opaque_exit_id!(
    /// Exact common-refinement map.
    RefinementMapIdV1
);
opaque_exit_id!(
    /// Theorem-card identity.
    ExitPathTheoremCardIdV1
);
opaque_exit_id!(
    /// External checker identity.
    ExitPathCheckerIdV1
);
opaque_exit_id!(
    /// Trusted-code-base declaration identity.
    ExitPathTcbIdV1
);
opaque_exit_id!(
    /// Exact retained hypothesis witness.
    ExitPathWitnessIdV1
);
opaque_exit_id!(
    /// Explicit no-claim artifact.
    ExitPathNoClaimIdV1
);
opaque_exit_id!(
    /// Preregistered counterexample identity.
    ExitPathFalsifierIdV1
);
opaque_exit_id!(
    /// Topological model appearing in a counterexample.
    ExitPathCountermodelIdV1
);

/// Domain-separated identity schema for one complete theorem-family snapshot.
pub enum ExitPathFamilySnapshotIdentitySchemaV1 {}

impl CanonicalSchema for ExitPathFamilySnapshotIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-geom.exit-path-family-snapshot.v1";
    const NAME: &'static str = "exit-path-approximation-constructibility-family-snapshot";
    const VERSION: u32 = EXIT_PATH_SCHEMA_VERSION_V1;
    const CONTEXT: &'static str = "derived model, stratification, direction and variance, link and monodromy fidelity, constructibility, properness, refinement, homotopy truncation, falsifiers, theorem state, TCB, and budget";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("subject", WireType::Bytes),
        FieldSpec::required("conventions", WireType::Bytes),
        FieldSpec::required("hypotheses", WireType::Bytes),
        FieldSpec::required("truncation-family", WireType::Bytes),
        FieldSpec::required("falsifiers", WireType::CanonicalSet),
        FieldSpec::required("theorem-state", WireType::Bytes),
        FieldSpec::required("tcb", WireType::Bytes),
        FieldSpec::required("budget", WireType::Bytes),
    ];
}

/// Typed identity of one complete RD.X1 statement/evidence/operation snapshot.
///
/// The lifecycle, TCB, falsifiers, and budgets are intentionally identity
/// bearing, so this is not a stable theorem-statement identity.
pub type ExitPathFamilySnapshotIdV1 = EvidenceNodeId<ExitPathFamilySnapshotIdentitySchemaV1>;

/// Admitted finite stratified-space class.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExitStratifiedSpaceClassV1 {
    /// Finite regular cell stratification.
    FiniteRegularCell,
    /// Finite conically stratified semialgebraic model.
    ConicalSemialgebraic,
    /// Finite conically stratified subanalytic model.
    ConicalSubanalytic,
    /// Unsupported nonconical, infinite, or unbounded presentation.
    Unsupported,
}

/// Direction of admissible stratified paths.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExitPathDirectionV1 {
    /// Paths may leave a stratum toward incident larger strata.
    Exit,
    /// Entrance-path convention, the categorical opposite of exit direction.
    Entrance,
}

/// Constructible (co)sheaf variance.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConstructibleVarianceV1 {
    /// Contravariant constructible sheaves.
    SheafContravariant,
    /// Covariant constructible cosheaves.
    CosheafCovariant,
}

/// Equivalence relation imposed on directed stratified paths.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StratifiedPathEquivalenceV1 {
    /// Endpoint-fixed stratified homotopy.
    EndpointFixed {
        /// Exact equivalence artifact.
        relation: StratifiedPathEquivalenceIdV1,
    },
    /// Thin stratified homotopy.
    Thin {
        /// Exact equivalence artifact.
        relation: StratifiedPathEquivalenceIdV1,
    },
    /// Higher coherent stratified homotopy through a finite degree.
    HigherThrough {
        /// Highest retained homotopy/coherence degree.
        degree: u8,
        /// Exact equivalence artifact.
        relation: StratifiedPathEquivalenceIdV1,
    },
    /// Full higher coherent stratified homotopy.
    FullHigher {
        /// Exact equivalence artifact.
        relation: StratifiedPathEquivalenceIdV1,
        /// Higher coherence.
        coherence: HigherCoherenceArtifactIdV1,
    },
    /// Equivalence is unresolved.
    Unknown {
        /// Explicit no-claim artifact.
        no_claim: ExitPathNoClaimIdV1,
    },
}

/// Local conical/link hypothesis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConicalLinkHypothesisV1 {
    /// Regular-cell links are contractible in every required slice.
    Contractible {
        /// Exact link catalog.
        links: ExitLinkCatalogIdV1,
        /// Contractibility witness.
        witness: ExitPathWitnessIdV1,
    },
    /// Link homotopy data is retained through a finite degree.
    RetainedThrough {
        /// Exact link catalog.
        links: ExitLinkCatalogIdV1,
        /// Highest retained degree.
        degree: u8,
        /// Checker witness.
        witness: ExitPathWitnessIdV1,
    },
    /// Full higher link object is retained.
    FullHigher {
        /// Exact link catalog.
        links: ExitLinkCatalogIdV1,
        /// Higher-coherence artifact.
        coherence: HigherCoherenceArtifactIdV1,
        /// Checker witness.
        witness: ExitPathWitnessIdV1,
    },
    /// Link fidelity is unresolved.
    Unknown {
        /// Explicit no-claim artifact.
        no_claim: ExitPathNoClaimIdV1,
    },
}

/// Stratum fundamental-group and local-system hypothesis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MonodromyHypothesisV1 {
    /// All relevant stratum fundamental groups and local monodromies are trivial.
    Trivial {
        /// Exact triviality witness.
        witness: ExitPathWitnessIdV1,
    },
    /// Within-stratum fundamental groupoids are retained.
    Groupoids {
        /// Exact groupoid catalog.
        groupoids: StratumGroupoidCatalogIdV1,
        /// Checker witness.
        witness: ExitPathWitnessIdV1,
    },
    /// Within-stratum groupoids and constructible local systems are retained
    /// through a degree.
    LocalSystemsThrough {
        /// Exact groupoid catalog.
        groupoids: StratumGroupoidCatalogIdV1,
        /// Local-system catalog.
        local_systems: LocalSystemCatalogIdV1,
        /// Highest retained coherence degree.
        degree: u8,
        /// Checker witness.
        witness: ExitPathWitnessIdV1,
    },
    /// Full higher monodromy/local-system data is retained.
    FullHigher {
        /// Exact groupoid catalog.
        groupoids: StratumGroupoidCatalogIdV1,
        /// Local-system catalog.
        local_systems: LocalSystemCatalogIdV1,
        /// Higher coherence.
        coherence: HigherCoherenceArtifactIdV1,
        /// Checker witness.
        witness: ExitPathWitnessIdV1,
    },
    /// Monodromy fidelity is unresolved.
    Unknown {
        /// Explicit no-claim artifact.
        no_claim: ExitPathNoClaimIdV1,
    },
}

/// Constructibility hypothesis on the chosen stratification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ConstructibilityHypothesisV1 {
    /// Locally constant on every stratum with retained restriction data.
    LocallyConstantOnStrata {
        /// Exact witness.
        witness: ExitPathWitnessIdV1,
    },
    /// Whitney/Thom-controlled constructibility.
    Controlled {
        /// Exact control witness.
        witness: ExitPathWitnessIdV1,
    },
    /// Constructibility is unresolved.
    Unknown {
        /// Explicit no-claim artifact.
        no_claim: ExitPathNoClaimIdV1,
    },
}

/// Compactness/properness hypothesis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExitPropernessHypothesisV1 {
    /// Compact stratified space.
    Compact {
        /// Compactness witness.
        witness: ExitPathWitnessIdV1,
    },
    /// Proper exit-path/descent map on a locally finite model.
    ProperLocallyFinite {
        /// Properness witness.
        witness: ExitPathWitnessIdV1,
    },
    /// Properness is unresolved.
    Unknown {
        /// Explicit no-claim artifact.
        no_claim: ExitPathNoClaimIdV1,
    },
}

/// Refinement-invariance hypothesis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RefinementHypothesisV1 {
    /// Exact identity refinement.
    Identity {
        /// Deterministic refinement identity.
        refinement: StratifiedRefinementIdV1,
    },
    /// Certified common refinement with forward and reverse comparison maps.
    CommonRefinement {
        /// Deterministic refinement identity.
        refinement: StratifiedRefinementIdV1,
        /// Forward map.
        forward: RefinementMapIdV1,
        /// Reverse map.
        reverse: RefinementMapIdV1,
        /// Equivalence witness.
        witness: ExitPathWitnessIdV1,
    },
    /// Refinement invariance is unresolved.
    Unknown {
        /// Explicit no-claim artifact.
        no_claim: ExitPathNoClaimIdV1,
    },
}

/// Global retained path/homotopy fidelity.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum HomotopyFidelityV1 {
    /// Incidence data only.
    IncidenceOnly,
    /// Homotopy and composition data through a finite degree.
    RetainedThrough {
        /// Highest retained degree.
        degree: u8,
        /// Exact coherence artifact.
        coherence: HigherCoherenceArtifactIdV1,
        /// Checker witness.
        witness: ExitPathWitnessIdV1,
    },
    /// Full higher coherent object retained.
    FullHigher {
        /// Exact coherence artifact.
        coherence: HigherCoherenceArtifactIdV1,
        /// Checker witness.
        witness: ExitPathWitnessIdV1,
    },
    /// Fidelity is unresolved.
    Unknown {
        /// Explicit no-claim artifact.
        no_claim: ExitPathNoClaimIdV1,
    },
}

/// Complete explicit hypothesis bundle.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ExitPathHypothesesV1 {
    /// Local conical/link data.
    pub links: ConicalLinkHypothesisV1,
    /// Stratum group and local-system data.
    pub monodromy: MonodromyHypothesisV1,
    /// Constructibility.
    pub constructibility: ConstructibilityHypothesisV1,
    /// Compactness/properness.
    pub properness: ExitPropernessHypothesisV1,
    /// Refinement invariance.
    pub refinement: RefinementHypothesisV1,
    /// Global path/homotopy fidelity.
    pub homotopy: HomotopyFidelityV1,
}

/// One target object in the maximal reduction lattice.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExitPathApproximationV1 {
    /// Incidence poset.
    IncidencePoset,
    /// Directed one-category enriched by within-stratum fundamental groupoids.
    StratumGroupoidEnrichedExitCategory,
    /// Simplicial category through a fixed simplex dimension.
    SimplicialCategory {
        /// Highest retained simplex dimension.
        max_simplex_dimension: u8,
    },
    /// General finite homotopy truncation.
    HigherTruncation {
        /// Highest retained homotopy degree.
        degree: u8,
    },
    /// Full exit/entrance path infinity-category statement.
    FullHigherCategory,
}

/// Why a reduction node remains Unknown.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExitPathUnknownReasonV1 {
    /// Base conical/constructibility/properness/refinement hypotheses are missing.
    CommonHypothesisMissing,
    /// Contractible or sufficiently retained link data is missing.
    LinkDataInsufficient,
    /// Fundamental-group/local-system fidelity is insufficient.
    MonodromyDataInsufficient,
    /// The declared path-equivalence relation is too weak for this degree.
    PathEquivalenceDataInsufficient,
    /// Global homotopy/coherence fidelity is insufficient.
    HomotopyDataInsufficient,
    /// Full higher coherence is not retained on every axis.
    FullHigherDataInsufficient,
}

/// Structural status of one theorem-lattice node.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExitPathNodeStateV1 {
    /// The implication is well formed under the declared hypotheses.
    ///
    /// This is a statement-admission status, not proof authority.
    SufficientStatement,
    /// The reduction cannot be stated as sufficient from the retained data.
    Unknown {
        /// Exact missing-data reason.
        reason: ExitPathUnknownReasonV1,
    },
    /// An unauthenticated counterexample record targets this exact node.
    RefutationRecorded {
        /// Recorded counterexample identity.
        falsifier: ExitPathFalsifierIdV1,
    },
}

/// One derived theorem-lattice node.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ExitPathTheoremNodeV1 {
    /// Candidate approximation.
    pub approximation: ExitPathApproximationV1,
    /// Structural sufficiency, precise Unknown, or node-scoped refutation.
    pub state: ExitPathNodeStateV1,
}

/// Preregistered counterexample class.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum ExitPathFalsifierKindV1 {
    /// Same incidence data but inequivalent links.
    SameIncidenceDifferentLink,
    /// Same incidence data but different fundamental-group/local-system monodromy.
    SameIncidenceDifferentMonodromy,
    /// Exit/entrance direction reversal changes variance/category.
    DirectionReversal,
    /// Deleting one named sufficient hypothesis breaks the conclusion.
    HypothesisDeletion,
}

/// Exact retained falsifier.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ExitPathFalsifierV1 {
    /// Falsifier identity.
    pub id: ExitPathFalsifierIdV1,
    /// Counterexample class.
    pub kind: ExitPathFalsifierKindV1,
    /// First topological/stratified model.
    pub left: ExitPathCountermodelIdV1,
    /// Second model or mutated twin.
    pub right: ExitPathCountermodelIdV1,
    /// Exact distinguishing witness.
    pub witness: ExitPathWitnessIdV1,
}

/// Theorem-family lifecycle state. V1 has no proved state.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExitPathTheoremStateV1 {
    /// Statement and falsifiers are preregistered.
    Preregistered {
        /// Exact theorem card.
        card: ExitPathTheoremCardIdV1,
    },
    /// Native/formal checking is in progress without authority.
    Candidate {
        /// Exact theorem card.
        card: ExitPathTheoremCardIdV1,
        /// Checker invocation witness.
        witness: ExitPathWitnessIdV1,
    },
    /// An unauthenticated counterexample record targets one exact node.
    RefutationRecorded {
        /// Exact targeted node; richer fallbacks remain independently eligible.
        approximation: ExitPathApproximationV1,
        /// Recorded counterexample.
        falsifier: ExitPathFalsifierIdV1,
    },
    /// No theorem-state claim.
    Unknown {
        /// Explicit no-claim artifact.
        no_claim: ExitPathNoClaimIdV1,
    },
}

/// Explicit trusted-code-base declaration.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ExitPathTcbV1 {
    /// TCB policy identity.
    pub tcb: ExitPathTcbIdV1,
    /// Independent checker implementation.
    pub checker: ExitPathCheckerIdV1,
    /// Exact theorem card.
    pub theorem_card: ExitPathTheoremCardIdV1,
}

/// Bounded structural analysis budget.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ExitPathBudgetV1 {
    /// Largest finite theorem node requested.
    pub max_truncation: u8,
    /// Maximum referenced artifact slots, counting distinct semantic roles.
    pub max_referenced_artifact_slots: u64,
    /// Maximum deterministic implication checks.
    pub max_implication_checks: u64,
    /// Positive wall-time declaration in seconds; not an execution receipt.
    pub declared_wall_seconds: f64,
}

/// Raw versioned RD.X1 theorem-family description.
#[derive(Debug, Clone, PartialEq)]
pub struct ExitPathFamilyIrV1 {
    schema_version: u32,
    /// Admitted derived geometry subject.
    pub geometry: DerivedGeometryIdV1,
    /// Exact underlying model version.
    pub model_version: DerivedModelVersionIdV1,
    /// Exact stratification.
    pub stratification: StratificationIdV1,
    /// State frame.
    pub frame: DerivedFrameIdV1,
    /// State units.
    pub units: DerivedUnitSystemIdV1,
    /// Finite admitted space class.
    pub space_class: ExitStratifiedSpaceClassV1,
    /// Exit or entrance direction.
    pub direction: ExitPathDirectionV1,
    /// Sheaf or cosheaf variance.
    pub variance: ConstructibleVarianceV1,
    /// Exact direction/endpoint convention bundle.
    pub convention: ExitPathConventionIdV1,
    /// Path equivalence.
    pub path_equivalence: StratifiedPathEquivalenceV1,
    /// Constructible (co)sheaf coefficient category, independent of the
    /// subject geometry's coordinate algebra.
    pub constructible_coefficients: CoefficientSystemV1,
    /// Explicit sufficient-hypothesis bundle.
    pub hypotheses: ExitPathHypothesesV1,
    /// Required counterexamples.
    pub falsifiers: Vec<ExitPathFalsifierV1>,
    /// Theorem lifecycle state.
    pub theorem_state: ExitPathTheoremStateV1,
    /// Trusted-code-base declaration.
    pub tcb: ExitPathTcbV1,
    /// Analysis budget and finite lattice ceiling.
    pub budget: ExitPathBudgetV1,
}

impl ExitPathFamilyIrV1 {
    /// Construct a current-version raw theorem family.
    #[allow(clippy::too_many_arguments)]
    #[must_use]
    pub fn new(
        geometry: DerivedGeometryIdV1,
        model_version: DerivedModelVersionIdV1,
        stratification: StratificationIdV1,
        frame: DerivedFrameIdV1,
        units: DerivedUnitSystemIdV1,
        space_class: ExitStratifiedSpaceClassV1,
        direction: ExitPathDirectionV1,
        variance: ConstructibleVarianceV1,
        convention: ExitPathConventionIdV1,
        path_equivalence: StratifiedPathEquivalenceV1,
        constructible_coefficients: CoefficientSystemV1,
        hypotheses: ExitPathHypothesesV1,
        falsifiers: Vec<ExitPathFalsifierV1>,
        theorem_state: ExitPathTheoremStateV1,
        tcb: ExitPathTcbV1,
        budget: ExitPathBudgetV1,
    ) -> Self {
        Self::with_schema_version(
            EXIT_PATH_SCHEMA_VERSION_V1,
            geometry,
            model_version,
            stratification,
            frame,
            units,
            space_class,
            direction,
            variance,
            convention,
            path_equivalence,
            constructible_coefficients,
            hypotheses,
            falsifiers,
            theorem_state,
            tcb,
            budget,
        )
    }

    /// Construct decoded versioned input. Unsupported versions fail closed.
    #[allow(clippy::too_many_arguments)]
    #[must_use]
    pub fn with_schema_version(
        schema_version: u32,
        geometry: DerivedGeometryIdV1,
        model_version: DerivedModelVersionIdV1,
        stratification: StratificationIdV1,
        frame: DerivedFrameIdV1,
        units: DerivedUnitSystemIdV1,
        space_class: ExitStratifiedSpaceClassV1,
        direction: ExitPathDirectionV1,
        variance: ConstructibleVarianceV1,
        convention: ExitPathConventionIdV1,
        path_equivalence: StratifiedPathEquivalenceV1,
        constructible_coefficients: CoefficientSystemV1,
        hypotheses: ExitPathHypothesesV1,
        falsifiers: Vec<ExitPathFalsifierV1>,
        theorem_state: ExitPathTheoremStateV1,
        tcb: ExitPathTcbV1,
        budget: ExitPathBudgetV1,
    ) -> Self {
        Self {
            schema_version,
            geometry,
            model_version,
            stratification,
            frame,
            units,
            space_class,
            direction,
            variance,
            convention,
            path_equivalence,
            constructible_coefficients,
            hypotheses,
            falsifiers,
            theorem_state,
            tcb,
            budget,
        }
    }

    /// Declared schema version.
    #[must_use]
    pub const fn schema_version(&self) -> u32 {
        self.schema_version
    }

    /// Number of externally identified artifact slots retained by this IR.
    ///
    /// Repeated digests in distinct semantic roles count separately because
    /// each role drives validation and canonical-identity work.
    #[must_use]
    pub fn required_referenced_artifact_slots(&self) -> u64 {
        referenced_artifact_slots(self)
    }
}

/// Invalid bounded numeric field.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExitPathFieldV1 {
    /// Path-equivalence truncation.
    PathEquivalenceDegree,
    /// Link truncation.
    LinkDegree,
    /// Monodromy/local-system truncation.
    MonodromyDegree,
    /// Global homotopy truncation.
    HomotopyDegree,
    /// Dyadic coefficient precision.
    CoefficientPrecision,
    /// Requested theorem-family truncation.
    MaxTruncation,
    /// Referenced-artifact-slot budget.
    ReferencedArtifactBudget,
    /// Implication-check budget.
    ImplicationBudget,
    /// Wall-time budget.
    WallTimeBudget,
}

/// Stable semantic group containing an absent all-zero identity sentinel.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExitPathIdentityFieldV1 {
    /// Geometry/model/stratification/frame/unit subject bundle.
    Subject,
    /// Exit/entrance convention artifact.
    Convention,
    /// Path-equivalence relation or coherence.
    PathEquivalence,
    /// Link catalog, coherence, witness, or no-claim artifact.
    LinkData,
    /// Groupoid/local-system/coherence/witness data.
    MonodromyData,
    /// Constructibility witness or no-claim artifact.
    Constructibility,
    /// Compactness/properness witness or no-claim artifact.
    Properness,
    /// Refinement, maps, witness, or no-claim artifact.
    Refinement,
    /// Homotopy coherence, witness, or no-claim artifact.
    Homotopy,
    /// Falsifier, countermodel, or distinguishing witness.
    Falsifier,
    /// Theorem card, candidate witness, refutation, or no-claim artifact.
    TheoremState,
    /// Trusted-code-base, checker, or theorem-card identity.
    TrustedCodeBase,
}

/// Stable field in the RD.1a admitted-subject binding.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExitPathSubjectBindingFieldV1 {
    /// Complete derived-geometry semantic identity.
    Geometry,
    /// Immutable source-model version.
    ModelVersion,
    /// Finite stratification identity.
    Stratification,
    /// Global coordinate frame.
    Frame,
    /// Global unit system.
    Units,
}

/// Deterministic structural refusal.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum ExitPathSemanticIssueV1 {
    /// Unsupported decoded schema.
    UnsupportedSchemaVersion {
        /// Supplied version.
        found: u32,
        /// Supported version.
        supported: u32,
    },
    /// Unsupported space presentation.
    UnsupportedSpaceClass,
    /// A bounded numeric field is zero, non-finite, or out of range.
    InvalidValue {
        /// Invalid field.
        field: ExitPathFieldV1,
    },
    /// A required typed identity is the all-zero sentinel.
    MissingIdentity {
        /// Semantic group containing the missing identity.
        field: ExitPathIdentityFieldV1,
    },
    /// Declared subject metadata does not belong to the supplied RD.1a object.
    SubjectBindingMismatch {
        /// Mismatched admitted-subject field.
        field: ExitPathSubjectBindingFieldV1,
    },
    /// Too many falsifiers were supplied before canonical work.
    TooManyFalsifiers {
        /// Supplied count.
        found: usize,
        /// Hard maximum.
        limit: usize,
    },
    /// The declared artifact budget cannot cover all referenced semantic slots.
    ReferencedArtifactBudgetExceeded {
        /// Slots required by the supplied IR.
        required: u64,
        /// Declared maximum.
        available: u64,
    },
    /// Required falsifier class is absent.
    MissingFalsifier {
        /// Missing class.
        kind: ExitPathFalsifierKindV1,
    },
    /// Falsifier identity is duplicated.
    DuplicateFalsifierId,
    /// Falsifier compares one model with itself.
    DegenerateFalsifier,
    /// Candidate/preregistered theorem card disagrees with the TCB card.
    TheoremCardMismatch,
    /// Refutation record names no retained falsifier.
    UnknownRecordedFalsifier,
    /// Refutation record names a node absent from the requested finite lattice.
    UnknownRecordedApproximation,
    /// Cooperative cancellation occurred before publication.
    Cancelled,
    /// Canonical identity construction failed.
    Identity(CanonicalError),
}

/// Complete deterministic refusal report.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ExitPathSemanticReportV1 {
    issues: Vec<ExitPathSemanticIssueV1>,
}

impl ExitPathSemanticReportV1 {
    fn new(issues: Vec<ExitPathSemanticIssueV1>) -> Self {
        Self { issues }
    }

    /// Ordered issues.
    #[must_use]
    pub fn issues(&self) -> &[ExitPathSemanticIssueV1] {
        &self.issues
    }
}

impl fmt::Display for ExitPathSemanticReportV1 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "exit-path theorem-family admission refused with {} issue(s)",
            self.issues.len()
        )
    }
}

impl core::error::Error for ExitPathSemanticReportV1 {}

/// Permanent v1 authority boundary.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ExitPathScientificAuthorityV1 {
    /// Statement admission and a sufficiency lattice are not theorem proof.
    ScientificCorrectnessNotProven,
}

/// Sealed theorem-family statement and its derived fallback lattice.
#[derive(Debug, Clone, PartialEq)]
pub struct ValidatedExitPathFamilyV1 {
    ir: ExitPathFamilyIrV1,
    receipt: IdentityReceipt<ExitPathFamilySnapshotIdV1>,
    lattice: Vec<ExitPathTheoremNodeV1>,
}

impl ValidatedExitPathFamilyV1 {
    /// Canonical admitted statement.
    #[must_use]
    pub const fn ir(&self) -> &ExitPathFamilyIrV1 {
        &self.ir
    }

    /// Typed identity of this complete statement/evidence/operation snapshot.
    #[must_use]
    pub const fn snapshot_id(&self) -> ExitPathFamilySnapshotIdV1 {
        self.receipt.id()
    }

    /// Exact identity/preimage receipt.
    #[must_use]
    pub const fn identity_receipt(&self) -> IdentityReceipt<ExitPathFamilySnapshotIdV1> {
        self.receipt
    }

    /// Poset-to-full-higher fallback lattice.
    #[must_use]
    pub fn theorem_lattice(&self) -> &[ExitPathTheoremNodeV1] {
        &self.lattice
    }

    /// Permanent no-authority boundary.
    #[must_use]
    pub const fn scientific_authority(&self) -> ExitPathScientificAuthorityV1 {
        ExitPathScientificAuthorityV1::ScientificCorrectnessNotProven
    }
}

/// Validate and canonically identify an RD.X1 theorem-family snapshot.
///
/// `admitted` is the authority-bearing RD.1a structural object. The redundant
/// geometry, model-version, stratification, frame, and unit IDs in `ir` must
/// match it exactly; a merely parsed [`DerivedGeometryIdV1`] is insufficient.
///
/// # Errors
/// Returns a deterministic report for schema, admitted-subject binding,
/// admitted class, truncation, budget, falsifier, theorem-state, cancellation,
/// or identity failures.
#[must_use = "the exit-path theorem-family admission result must be handled"]
pub fn validate_exit_path_family_v1(
    mut ir: ExitPathFamilyIrV1,
    admitted: &AdmittedDerivedGeometryV1,
    cx: &Cx<'_>,
) -> Result<ValidatedExitPathFamilyV1, ExitPathSemanticReportV1> {
    checkpoint(cx)?;
    if ir.falsifiers.len() > MAX_EXIT_PATH_FALSIFIERS_V1 {
        return Err(ExitPathSemanticReportV1::new(vec![
            ExitPathSemanticIssueV1::TooManyFalsifiers {
                found: ir.falsifiers.len(),
                limit: MAX_EXIT_PATH_FALSIFIERS_V1,
            },
        ]));
    }
    canonicalize_zero(&mut ir.budget.declared_wall_seconds);
    let mut issues = Vec::new();
    if ir.schema_version != EXIT_PATH_SCHEMA_VERSION_V1 {
        issues.push(ExitPathSemanticIssueV1::UnsupportedSchemaVersion {
            found: ir.schema_version,
            supported: EXIT_PATH_SCHEMA_VERSION_V1,
        });
    }
    if matches!(ir.space_class, ExitStratifiedSpaceClassV1::Unsupported) {
        issues.push(ExitPathSemanticIssueV1::UnsupportedSpaceClass);
    }
    validate_subject_binding(&ir, admitted, &mut issues);
    validate_identities(&ir, &mut issues);
    validate_degrees_and_budget(&ir, &mut issues);

    ir.falsifiers.sort_by_key(|falsifier| falsifier.id);
    if ir
        .falsifiers
        .windows(2)
        .any(|pair| pair[0].id == pair[1].id)
    {
        issues.push(ExitPathSemanticIssueV1::DuplicateFalsifierId);
    }
    if ir
        .falsifiers
        .iter()
        .any(|falsifier| falsifier.left == falsifier.right)
    {
        issues.push(ExitPathSemanticIssueV1::DegenerateFalsifier);
    }
    for required in [
        ExitPathFalsifierKindV1::SameIncidenceDifferentLink,
        ExitPathFalsifierKindV1::SameIncidenceDifferentMonodromy,
        ExitPathFalsifierKindV1::DirectionReversal,
        ExitPathFalsifierKindV1::HypothesisDeletion,
    ] {
        if !ir
            .falsifiers
            .iter()
            .any(|falsifier| falsifier.kind == required)
        {
            issues.push(ExitPathSemanticIssueV1::MissingFalsifier { kind: required });
        }
    }
    validate_theorem_state(&ir, &mut issues);
    if !issues.is_empty() {
        return Err(ExitPathSemanticReportV1::new(issues));
    }
    checkpoint(cx)?;
    let lattice = derive_theorem_lattice(&ir, cx)?;
    let receipt = exit_path_receipt(&ir, cx).map_err(identity_report)?;
    Ok(ValidatedExitPathFamilyV1 {
        ir,
        receipt,
        lattice,
    })
}

fn checkpoint(cx: &Cx<'_>) -> Result<(), ExitPathSemanticReportV1> {
    cx.checkpoint()
        .map_err(|_| ExitPathSemanticReportV1::new(vec![ExitPathSemanticIssueV1::Cancelled]))
}

fn identity_report(error: CanonicalError) -> ExitPathSemanticReportV1 {
    let issue = if matches!(&error, CanonicalError::Cancelled { .. }) {
        ExitPathSemanticIssueV1::Cancelled
    } else {
        ExitPathSemanticIssueV1::Identity(error)
    };
    ExitPathSemanticReportV1::new(vec![issue])
}

fn canonicalize_zero(value: &mut f64) {
    if *value == 0.0 {
        *value = 0.0;
    }
}

fn is_zero(bytes: &[u8; 32]) -> bool {
    bytes.iter().all(|byte| *byte == 0)
}

fn validate_subject_binding(
    ir: &ExitPathFamilyIrV1,
    admitted: &AdmittedDerivedGeometryV1,
    issues: &mut Vec<ExitPathSemanticIssueV1>,
) {
    let admitted_ir = admitted.ir();
    for (mismatch, field) in [
        (
            ir.geometry != admitted.id(),
            ExitPathSubjectBindingFieldV1::Geometry,
        ),
        (
            ir.model_version != admitted_ir.model_version,
            ExitPathSubjectBindingFieldV1::ModelVersion,
        ),
        (
            ir.stratification != admitted_ir.stratification.id,
            ExitPathSubjectBindingFieldV1::Stratification,
        ),
        (
            ir.frame != admitted_ir.frame,
            ExitPathSubjectBindingFieldV1::Frame,
        ),
        (
            ir.units != admitted_ir.unit_system,
            ExitPathSubjectBindingFieldV1::Units,
        ),
    ] {
        if mismatch {
            issues.push(ExitPathSemanticIssueV1::SubjectBindingMismatch { field });
        }
    }
}

fn validate_identities(ir: &ExitPathFamilyIrV1, issues: &mut Vec<ExitPathSemanticIssueV1>) {
    if is_zero(ir.geometry.as_bytes())
        || is_zero(ir.model_version.as_bytes())
        || is_zero(ir.stratification.as_bytes())
        || is_zero(ir.frame.as_bytes())
        || is_zero(ir.units.as_bytes())
    {
        issues.push(ExitPathSemanticIssueV1::MissingIdentity {
            field: ExitPathIdentityFieldV1::Subject,
        });
    }
    if ir.convention.is_zero() {
        issues.push(ExitPathSemanticIssueV1::MissingIdentity {
            field: ExitPathIdentityFieldV1::Convention,
        });
    }
    if path_equivalence_identity_missing(ir.path_equivalence) {
        issues.push(ExitPathSemanticIssueV1::MissingIdentity {
            field: ExitPathIdentityFieldV1::PathEquivalence,
        });
    }
    for (missing, field) in [
        (
            link_identity_missing(ir.hypotheses.links),
            ExitPathIdentityFieldV1::LinkData,
        ),
        (
            monodromy_identity_missing(ir.hypotheses.monodromy),
            ExitPathIdentityFieldV1::MonodromyData,
        ),
        (
            constructibility_identity_missing(ir.hypotheses.constructibility),
            ExitPathIdentityFieldV1::Constructibility,
        ),
        (
            properness_identity_missing(ir.hypotheses.properness),
            ExitPathIdentityFieldV1::Properness,
        ),
        (
            refinement_identity_missing(ir.hypotheses.refinement),
            ExitPathIdentityFieldV1::Refinement,
        ),
        (
            homotopy_identity_missing(ir.hypotheses.homotopy),
            ExitPathIdentityFieldV1::Homotopy,
        ),
    ] {
        if missing {
            issues.push(ExitPathSemanticIssueV1::MissingIdentity { field });
        }
    }
    if ir.falsifiers.iter().any(|falsifier| {
        falsifier.id.is_zero()
            || falsifier.left.is_zero()
            || falsifier.right.is_zero()
            || falsifier.witness.is_zero()
    }) {
        issues.push(ExitPathSemanticIssueV1::MissingIdentity {
            field: ExitPathIdentityFieldV1::Falsifier,
        });
    }
    if theorem_state_identity_missing(ir.theorem_state) {
        issues.push(ExitPathSemanticIssueV1::MissingIdentity {
            field: ExitPathIdentityFieldV1::TheoremState,
        });
    }
    if ir.tcb.tcb.is_zero() || ir.tcb.checker.is_zero() || ir.tcb.theorem_card.is_zero() {
        issues.push(ExitPathSemanticIssueV1::MissingIdentity {
            field: ExitPathIdentityFieldV1::TrustedCodeBase,
        });
    }
}

fn path_equivalence_identity_missing(value: StratifiedPathEquivalenceV1) -> bool {
    match value {
        StratifiedPathEquivalenceV1::EndpointFixed { relation }
        | StratifiedPathEquivalenceV1::Thin { relation }
        | StratifiedPathEquivalenceV1::HigherThrough { relation, .. } => relation.is_zero(),
        StratifiedPathEquivalenceV1::FullHigher {
            relation,
            coherence,
        } => relation.is_zero() || coherence.is_zero(),
        StratifiedPathEquivalenceV1::Unknown { no_claim } => no_claim.is_zero(),
    }
}

fn link_identity_missing(value: ConicalLinkHypothesisV1) -> bool {
    match value {
        ConicalLinkHypothesisV1::Contractible { links, witness }
        | ConicalLinkHypothesisV1::RetainedThrough { links, witness, .. } => {
            links.is_zero() || witness.is_zero()
        }
        ConicalLinkHypothesisV1::FullHigher {
            links,
            coherence,
            witness,
        } => links.is_zero() || coherence.is_zero() || witness.is_zero(),
        ConicalLinkHypothesisV1::Unknown { no_claim } => no_claim.is_zero(),
    }
}

fn monodromy_identity_missing(value: MonodromyHypothesisV1) -> bool {
    match value {
        MonodromyHypothesisV1::Trivial { witness } => witness.is_zero(),
        MonodromyHypothesisV1::Groupoids { groupoids, witness } => {
            groupoids.is_zero() || witness.is_zero()
        }
        MonodromyHypothesisV1::LocalSystemsThrough {
            groupoids,
            local_systems,
            witness,
            ..
        } => groupoids.is_zero() || local_systems.is_zero() || witness.is_zero(),
        MonodromyHypothesisV1::FullHigher {
            groupoids,
            local_systems,
            coherence,
            witness,
        } => {
            groupoids.is_zero()
                || local_systems.is_zero()
                || coherence.is_zero()
                || witness.is_zero()
        }
        MonodromyHypothesisV1::Unknown { no_claim } => no_claim.is_zero(),
    }
}

fn constructibility_identity_missing(value: ConstructibilityHypothesisV1) -> bool {
    match value {
        ConstructibilityHypothesisV1::LocallyConstantOnStrata { witness }
        | ConstructibilityHypothesisV1::Controlled { witness } => witness.is_zero(),
        ConstructibilityHypothesisV1::Unknown { no_claim } => no_claim.is_zero(),
    }
}

fn properness_identity_missing(value: ExitPropernessHypothesisV1) -> bool {
    match value {
        ExitPropernessHypothesisV1::Compact { witness }
        | ExitPropernessHypothesisV1::ProperLocallyFinite { witness } => witness.is_zero(),
        ExitPropernessHypothesisV1::Unknown { no_claim } => no_claim.is_zero(),
    }
}

fn refinement_identity_missing(value: RefinementHypothesisV1) -> bool {
    match value {
        RefinementHypothesisV1::Identity { refinement } => refinement.is_zero(),
        RefinementHypothesisV1::CommonRefinement {
            refinement,
            forward,
            reverse,
            witness,
        } => refinement.is_zero() || forward.is_zero() || reverse.is_zero() || witness.is_zero(),
        RefinementHypothesisV1::Unknown { no_claim } => no_claim.is_zero(),
    }
}

fn homotopy_identity_missing(value: HomotopyFidelityV1) -> bool {
    match value {
        HomotopyFidelityV1::IncidenceOnly => false,
        HomotopyFidelityV1::RetainedThrough {
            coherence, witness, ..
        }
        | HomotopyFidelityV1::FullHigher { coherence, witness } => {
            coherence.is_zero() || witness.is_zero()
        }
        HomotopyFidelityV1::Unknown { no_claim } => no_claim.is_zero(),
    }
}

fn theorem_state_identity_missing(value: ExitPathTheoremStateV1) -> bool {
    match value {
        ExitPathTheoremStateV1::Preregistered { card } => card.is_zero(),
        ExitPathTheoremStateV1::Candidate { card, witness } => card.is_zero() || witness.is_zero(),
        ExitPathTheoremStateV1::RefutationRecorded { falsifier, .. } => falsifier.is_zero(),
        ExitPathTheoremStateV1::Unknown { no_claim } => no_claim.is_zero(),
    }
}

fn validate_degree(degree: u8, field: ExitPathFieldV1, issues: &mut Vec<ExitPathSemanticIssueV1>) {
    if degree > MAX_EXIT_PATH_TRUNCATION_V1 {
        issues.push(ExitPathSemanticIssueV1::InvalidValue { field });
    }
}

fn validate_degrees_and_budget(ir: &ExitPathFamilyIrV1, issues: &mut Vec<ExitPathSemanticIssueV1>) {
    match ir.path_equivalence {
        StratifiedPathEquivalenceV1::HigherThrough { degree, .. } => {
            validate_degree(degree, ExitPathFieldV1::PathEquivalenceDegree, issues);
        }
        StratifiedPathEquivalenceV1::EndpointFixed { .. }
        | StratifiedPathEquivalenceV1::Thin { .. }
        | StratifiedPathEquivalenceV1::FullHigher { .. }
        | StratifiedPathEquivalenceV1::Unknown { .. } => {}
    }
    match ir.hypotheses.links {
        ConicalLinkHypothesisV1::RetainedThrough { degree, .. } => {
            validate_degree(degree, ExitPathFieldV1::LinkDegree, issues);
        }
        ConicalLinkHypothesisV1::Contractible { .. }
        | ConicalLinkHypothesisV1::FullHigher { .. }
        | ConicalLinkHypothesisV1::Unknown { .. } => {}
    }
    if let MonodromyHypothesisV1::LocalSystemsThrough { degree, .. } = ir.hypotheses.monodromy {
        validate_degree(degree, ExitPathFieldV1::MonodromyDegree, issues);
    }
    if let HomotopyFidelityV1::RetainedThrough { degree, .. } = ir.hypotheses.homotopy {
        validate_degree(degree, ExitPathFieldV1::HomotopyDegree, issues);
    }
    if let CoefficientSystemV1::DyadicIntervalReal { precision_bits: 0 } =
        ir.constructible_coefficients
    {
        issues.push(ExitPathSemanticIssueV1::InvalidValue {
            field: ExitPathFieldV1::CoefficientPrecision,
        });
    }
    if ir.budget.max_truncation > MAX_EXIT_PATH_TRUNCATION_V1 {
        issues.push(ExitPathSemanticIssueV1::InvalidValue {
            field: ExitPathFieldV1::MaxTruncation,
        });
    }
    if ir.budget.max_referenced_artifact_slots == 0 {
        issues.push(ExitPathSemanticIssueV1::InvalidValue {
            field: ExitPathFieldV1::ReferencedArtifactBudget,
        });
    } else {
        let required = referenced_artifact_slots(ir);
        if ir.budget.max_referenced_artifact_slots < required {
            issues.push(ExitPathSemanticIssueV1::ReferencedArtifactBudgetExceeded {
                required,
                available: ir.budget.max_referenced_artifact_slots,
            });
        }
    }
    if ir.budget.max_implication_checks == 0 {
        issues.push(ExitPathSemanticIssueV1::InvalidValue {
            field: ExitPathFieldV1::ImplicationBudget,
        });
    } else if ir.budget.max_implication_checks < lattice_node_count(ir.budget.max_truncation) {
        issues.push(ExitPathSemanticIssueV1::InvalidValue {
            field: ExitPathFieldV1::ImplicationBudget,
        });
    }
    if !(ir.budget.declared_wall_seconds.is_finite() && ir.budget.declared_wall_seconds > 0.0) {
        issues.push(ExitPathSemanticIssueV1::InvalidValue {
            field: ExitPathFieldV1::WallTimeBudget,
        });
    }
}

fn referenced_artifact_slots(ir: &ExitPathFamilyIrV1) -> u64 {
    let path_equivalence = match ir.path_equivalence {
        StratifiedPathEquivalenceV1::FullHigher { .. } => 2,
        StratifiedPathEquivalenceV1::EndpointFixed { .. }
        | StratifiedPathEquivalenceV1::Thin { .. }
        | StratifiedPathEquivalenceV1::HigherThrough { .. }
        | StratifiedPathEquivalenceV1::Unknown { .. } => 1,
    };
    let links = match ir.hypotheses.links {
        ConicalLinkHypothesisV1::Contractible { .. }
        | ConicalLinkHypothesisV1::RetainedThrough { .. } => 2,
        ConicalLinkHypothesisV1::FullHigher { .. } => 3,
        ConicalLinkHypothesisV1::Unknown { .. } => 1,
    };
    let monodromy = match ir.hypotheses.monodromy {
        MonodromyHypothesisV1::Trivial { .. } | MonodromyHypothesisV1::Unknown { .. } => 1,
        MonodromyHypothesisV1::Groupoids { .. } => 2,
        MonodromyHypothesisV1::LocalSystemsThrough { .. } => 3,
        MonodromyHypothesisV1::FullHigher { .. } => 4,
    };
    let refinement = match ir.hypotheses.refinement {
        RefinementHypothesisV1::Identity { .. } | RefinementHypothesisV1::Unknown { .. } => 1,
        RefinementHypothesisV1::CommonRefinement { .. } => 4,
    };
    let homotopy = match ir.hypotheses.homotopy {
        HomotopyFidelityV1::IncidenceOnly => 0,
        HomotopyFidelityV1::RetainedThrough { .. } | HomotopyFidelityV1::FullHigher { .. } => 2,
        HomotopyFidelityV1::Unknown { .. } => 1,
    };
    let theorem_state = match ir.theorem_state {
        ExitPathTheoremStateV1::Candidate { .. } => 2,
        ExitPathTheoremStateV1::Preregistered { .. }
        | ExitPathTheoremStateV1::RefutationRecorded { .. }
        | ExitPathTheoremStateV1::Unknown { .. } => 1,
    };

    // Five subject IDs, one convention, one constructibility witness/no-claim,
    // one properness witness/no-claim, three TCB IDs, and four artifact slots
    // per bounded falsifier. The hard falsifier cap makes every conversion and
    // multiplication below exact in u64.
    11 + path_equivalence
        + links
        + monodromy
        + refinement
        + homotopy
        + theorem_state
        + 4 * ir.falsifiers.len() as u64
}

fn validate_theorem_state(ir: &ExitPathFamilyIrV1, issues: &mut Vec<ExitPathSemanticIssueV1>) {
    match ir.theorem_state {
        ExitPathTheoremStateV1::Preregistered { card }
        | ExitPathTheoremStateV1::Candidate { card, .. } => {
            if card != ir.tcb.theorem_card {
                issues.push(ExitPathSemanticIssueV1::TheoremCardMismatch);
            }
        }
        ExitPathTheoremStateV1::RefutationRecorded {
            approximation,
            falsifier,
        } => {
            if !ir
                .falsifiers
                .iter()
                .any(|candidate| candidate.id == falsifier)
            {
                issues.push(ExitPathSemanticIssueV1::UnknownRecordedFalsifier);
            }
            if !approximation_requested(approximation, ir.budget.max_truncation) {
                issues.push(ExitPathSemanticIssueV1::UnknownRecordedApproximation);
            }
        }
        ExitPathTheoremStateV1::Unknown { .. } => {}
    }
}

fn approximation_requested(approximation: ExitPathApproximationV1, max_truncation: u8) -> bool {
    match approximation {
        ExitPathApproximationV1::IncidencePoset | ExitPathApproximationV1::FullHigherCategory => {
            true
        }
        ExitPathApproximationV1::StratumGroupoidEnrichedExitCategory => max_truncation >= 1,
        ExitPathApproximationV1::SimplicialCategory {
            max_simplex_dimension,
        } => max_truncation >= 2 && max_simplex_dimension == 2,
        ExitPathApproximationV1::HigherTruncation { degree } => {
            degree >= 2 && degree <= max_truncation
        }
    }
}

fn common_hypotheses_hold(ir: &ExitPathFamilyIrV1) -> bool {
    !matches!(ir.space_class, ExitStratifiedSpaceClassV1::Unsupported)
        && !matches!(
            ir.path_equivalence,
            StratifiedPathEquivalenceV1::Unknown { .. }
        )
        && !matches!(
            ir.hypotheses.constructibility,
            ConstructibilityHypothesisV1::Unknown { .. }
        )
        && !matches!(
            ir.hypotheses.properness,
            ExitPropernessHypothesisV1::Unknown { .. }
        )
        && !matches!(
            ir.hypotheses.refinement,
            RefinementHypothesisV1::Unknown { .. }
        )
}

fn finite_link_support(links: ConicalLinkHypothesisV1, degree: u8) -> bool {
    match links {
        ConicalLinkHypothesisV1::Contractible { .. }
        | ConicalLinkHypothesisV1::FullHigher { .. } => true,
        ConicalLinkHypothesisV1::RetainedThrough {
            degree: retained, ..
        } => retained >= degree,
        ConicalLinkHypothesisV1::Unknown { .. } => false,
    }
}

fn finite_monodromy_support(monodromy: MonodromyHypothesisV1, degree: u8) -> bool {
    match monodromy {
        MonodromyHypothesisV1::Trivial { .. } | MonodromyHypothesisV1::Groupoids { .. } => {
            degree <= 1
        }
        MonodromyHypothesisV1::FullHigher { .. } => true,
        MonodromyHypothesisV1::LocalSystemsThrough {
            degree: retained, ..
        } => retained >= degree,
        MonodromyHypothesisV1::Unknown { .. } => false,
    }
}

fn finite_homotopy_support(homotopy: HomotopyFidelityV1, degree: u8) -> bool {
    match homotopy {
        HomotopyFidelityV1::IncidenceOnly => degree == 0,
        HomotopyFidelityV1::RetainedThrough {
            degree: retained, ..
        } => retained >= degree,
        HomotopyFidelityV1::FullHigher { .. } => true,
        HomotopyFidelityV1::Unknown { .. } => false,
    }
}

fn finite_path_equivalence_support(equivalence: StratifiedPathEquivalenceV1, degree: u8) -> bool {
    match equivalence {
        StratifiedPathEquivalenceV1::EndpointFixed { .. }
        | StratifiedPathEquivalenceV1::Thin { .. } => degree <= 1,
        StratifiedPathEquivalenceV1::HigherThrough {
            degree: retained, ..
        } => retained >= degree,
        StratifiedPathEquivalenceV1::FullHigher { .. } => true,
        StratifiedPathEquivalenceV1::Unknown { .. } => false,
    }
}

fn node_for(
    ir: &ExitPathFamilyIrV1,
    approximation: ExitPathApproximationV1,
) -> ExitPathTheoremNodeV1 {
    let state = if !common_hypotheses_hold(ir) {
        ExitPathNodeStateV1::Unknown {
            reason: ExitPathUnknownReasonV1::CommonHypothesisMissing,
        }
    } else {
        match approximation {
            ExitPathApproximationV1::IncidencePoset => {
                if !matches!(
                    ir.hypotheses.links,
                    ConicalLinkHypothesisV1::Contractible { .. }
                ) {
                    ExitPathNodeStateV1::Unknown {
                        reason: ExitPathUnknownReasonV1::LinkDataInsufficient,
                    }
                } else if !matches!(
                    ir.hypotheses.monodromy,
                    MonodromyHypothesisV1::Trivial { .. }
                ) {
                    ExitPathNodeStateV1::Unknown {
                        reason: ExitPathUnknownReasonV1::MonodromyDataInsufficient,
                    }
                } else if !finite_homotopy_support(ir.hypotheses.homotopy, 0) {
                    ExitPathNodeStateV1::Unknown {
                        reason: ExitPathUnknownReasonV1::HomotopyDataInsufficient,
                    }
                } else {
                    ExitPathNodeStateV1::SufficientStatement
                }
            }
            ExitPathApproximationV1::StratumGroupoidEnrichedExitCategory => {
                finite_node_state(ir, 1)
            }
            ExitPathApproximationV1::SimplicialCategory {
                max_simplex_dimension,
            } => finite_node_state(ir, max_simplex_dimension),
            ExitPathApproximationV1::HigherTruncation { degree } => finite_node_state(ir, degree),
            ExitPathApproximationV1::FullHigherCategory => {
                let full_equivalence = matches!(
                    ir.path_equivalence,
                    StratifiedPathEquivalenceV1::FullHigher { .. }
                );
                let full_links = matches!(
                    ir.hypotheses.links,
                    ConicalLinkHypothesisV1::Contractible { .. }
                        | ConicalLinkHypothesisV1::FullHigher { .. }
                );
                let full_monodromy = matches!(
                    ir.hypotheses.monodromy,
                    MonodromyHypothesisV1::FullHigher { .. }
                );
                let full_homotopy = matches!(
                    ir.hypotheses.homotopy,
                    HomotopyFidelityV1::FullHigher { .. }
                );
                if full_equivalence && full_links && full_monodromy && full_homotopy {
                    ExitPathNodeStateV1::SufficientStatement
                } else {
                    ExitPathNodeStateV1::Unknown {
                        reason: ExitPathUnknownReasonV1::FullHigherDataInsufficient,
                    }
                }
            }
        }
    };
    let state = match ir.theorem_state {
        ExitPathTheoremStateV1::RefutationRecorded {
            approximation: refuted,
            falsifier,
        } if refuted == approximation => ExitPathNodeStateV1::RefutationRecorded { falsifier },
        ExitPathTheoremStateV1::Preregistered { .. }
        | ExitPathTheoremStateV1::Candidate { .. }
        | ExitPathTheoremStateV1::RefutationRecorded { .. }
        | ExitPathTheoremStateV1::Unknown { .. } => state,
    };
    ExitPathTheoremNodeV1 {
        approximation,
        state,
    }
}

fn finite_node_state(ir: &ExitPathFamilyIrV1, degree: u8) -> ExitPathNodeStateV1 {
    if !finite_path_equivalence_support(ir.path_equivalence, degree) {
        ExitPathNodeStateV1::Unknown {
            reason: ExitPathUnknownReasonV1::PathEquivalenceDataInsufficient,
        }
    } else if !finite_link_support(ir.hypotheses.links, degree) {
        ExitPathNodeStateV1::Unknown {
            reason: ExitPathUnknownReasonV1::LinkDataInsufficient,
        }
    } else if !finite_monodromy_support(ir.hypotheses.monodromy, degree) {
        ExitPathNodeStateV1::Unknown {
            reason: ExitPathUnknownReasonV1::MonodromyDataInsufficient,
        }
    } else if !finite_homotopy_support(ir.hypotheses.homotopy, degree) {
        ExitPathNodeStateV1::Unknown {
            reason: ExitPathUnknownReasonV1::HomotopyDataInsufficient,
        }
    } else {
        ExitPathNodeStateV1::SufficientStatement
    }
}

fn derive_theorem_lattice(
    ir: &ExitPathFamilyIrV1,
    cx: &Cx<'_>,
) -> Result<Vec<ExitPathTheoremNodeV1>, ExitPathSemanticReportV1> {
    let mut lattice = Vec::with_capacity(usize::from(ir.budget.max_truncation) + 4);
    push_lattice_node(
        &mut lattice,
        ir,
        ExitPathApproximationV1::IncidencePoset,
        cx,
    )?;
    if ir.budget.max_truncation >= 1 {
        push_lattice_node(
            &mut lattice,
            ir,
            ExitPathApproximationV1::StratumGroupoidEnrichedExitCategory,
            cx,
        )?;
    }
    if ir.budget.max_truncation >= 2 {
        push_lattice_node(
            &mut lattice,
            ir,
            ExitPathApproximationV1::SimplicialCategory {
                max_simplex_dimension: 2,
            },
            cx,
        )?;
    }
    for degree in 2..=ir.budget.max_truncation {
        push_lattice_node(
            &mut lattice,
            ir,
            ExitPathApproximationV1::HigherTruncation { degree },
            cx,
        )?;
    }
    push_lattice_node(
        &mut lattice,
        ir,
        ExitPathApproximationV1::FullHigherCategory,
        cx,
    )?;
    Ok(lattice)
}

fn push_lattice_node(
    lattice: &mut Vec<ExitPathTheoremNodeV1>,
    ir: &ExitPathFamilyIrV1,
    approximation: ExitPathApproximationV1,
    cx: &Cx<'_>,
) -> Result<(), ExitPathSemanticReportV1> {
    checkpoint(cx)?;
    lattice.push(node_for(ir, approximation));
    Ok(())
}

fn lattice_node_count(max_truncation: u8) -> u64 {
    match max_truncation {
        0 => 2,
        1 => 3,
        degree => u64::from(degree) + 3,
    }
}

fn push_u16(out: &mut Vec<u8>, value: u16) {
    out.extend_from_slice(&value.to_le_bytes());
}

fn push_u64(out: &mut Vec<u8>, value: u64) {
    out.extend_from_slice(&value.to_le_bytes());
}

fn push_f64(out: &mut Vec<u8>, value: f64) {
    let bits = if value == 0.0 {
        0.0_f64.to_bits()
    } else {
        value.to_bits()
    };
    push_u64(out, bits);
}

fn push_digest<I: DigestBytes>(out: &mut Vec<u8>, id: I) {
    out.extend_from_slice(id.digest_bytes());
}

fn subject_bytes(ir: &ExitPathFamilyIrV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(200);
    out.extend_from_slice(ir.geometry.as_bytes());
    out.extend_from_slice(ir.model_version.as_bytes());
    out.extend_from_slice(ir.stratification.as_bytes());
    out.extend_from_slice(ir.frame.as_bytes());
    out.extend_from_slice(ir.units.as_bytes());
    out.push(match ir.space_class {
        ExitStratifiedSpaceClassV1::FiniteRegularCell => 0,
        ExitStratifiedSpaceClassV1::ConicalSemialgebraic => 1,
        ExitStratifiedSpaceClassV1::ConicalSubanalytic => 2,
        ExitStratifiedSpaceClassV1::Unsupported => 3,
    });
    out
}

fn constructible_coefficient_bytes(coefficients: CoefficientSystemV1, out: &mut Vec<u8>) {
    match coefficients {
        CoefficientSystemV1::RationalReal => out.push(0),
        CoefficientSystemV1::AlgebraicReal => out.push(1),
        CoefficientSystemV1::DyadicIntervalReal { precision_bits } => {
            out.push(2);
            push_u16(out, precision_bits);
        }
        CoefficientSystemV1::RationalComplex => out.push(3),
        CoefficientSystemV1::AlgebraicComplex => out.push(4),
    }
}

fn conventions_bytes(ir: &ExitPathFamilyIrV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(160);
    out.push(match ir.direction {
        ExitPathDirectionV1::Exit => 0,
        ExitPathDirectionV1::Entrance => 1,
    });
    out.push(match ir.variance {
        ConstructibleVarianceV1::SheafContravariant => 0,
        ConstructibleVarianceV1::CosheafCovariant => 1,
    });
    push_digest(&mut out, ir.convention);
    match ir.path_equivalence {
        StratifiedPathEquivalenceV1::EndpointFixed { relation } => {
            out.push(0);
            push_digest(&mut out, relation);
        }
        StratifiedPathEquivalenceV1::Thin { relation } => {
            out.push(1);
            push_digest(&mut out, relation);
        }
        StratifiedPathEquivalenceV1::HigherThrough { degree, relation } => {
            out.push(2);
            out.push(degree);
            push_digest(&mut out, relation);
        }
        StratifiedPathEquivalenceV1::FullHigher {
            relation,
            coherence,
        } => {
            out.push(3);
            push_digest(&mut out, relation);
            push_digest(&mut out, coherence);
        }
        StratifiedPathEquivalenceV1::Unknown { no_claim } => {
            out.push(4);
            push_digest(&mut out, no_claim);
        }
    }
    constructible_coefficient_bytes(ir.constructible_coefficients, &mut out);
    out
}

fn push_links(out: &mut Vec<u8>, links: ConicalLinkHypothesisV1) {
    match links {
        ConicalLinkHypothesisV1::Contractible { links, witness } => {
            out.push(0);
            push_digest(out, links);
            push_digest(out, witness);
        }
        ConicalLinkHypothesisV1::RetainedThrough {
            links,
            degree,
            witness,
        } => {
            out.push(1);
            push_digest(out, links);
            out.push(degree);
            push_digest(out, witness);
        }
        ConicalLinkHypothesisV1::FullHigher {
            links,
            coherence,
            witness,
        } => {
            out.push(2);
            push_digest(out, links);
            push_digest(out, coherence);
            push_digest(out, witness);
        }
        ConicalLinkHypothesisV1::Unknown { no_claim } => {
            out.push(3);
            push_digest(out, no_claim);
        }
    }
}

fn push_monodromy(out: &mut Vec<u8>, monodromy: MonodromyHypothesisV1) {
    match monodromy {
        MonodromyHypothesisV1::Trivial { witness } => {
            out.push(0);
            push_digest(out, witness);
        }
        MonodromyHypothesisV1::Groupoids { groupoids, witness } => {
            out.push(1);
            push_digest(out, groupoids);
            push_digest(out, witness);
        }
        MonodromyHypothesisV1::LocalSystemsThrough {
            groupoids,
            local_systems,
            degree,
            witness,
        } => {
            out.push(2);
            push_digest(out, groupoids);
            push_digest(out, local_systems);
            out.push(degree);
            push_digest(out, witness);
        }
        MonodromyHypothesisV1::FullHigher {
            groupoids,
            local_systems,
            coherence,
            witness,
        } => {
            out.push(3);
            push_digest(out, groupoids);
            push_digest(out, local_systems);
            push_digest(out, coherence);
            push_digest(out, witness);
        }
        MonodromyHypothesisV1::Unknown { no_claim } => {
            out.push(4);
            push_digest(out, no_claim);
        }
    }
}

fn push_constructibility(out: &mut Vec<u8>, value: ConstructibilityHypothesisV1) {
    match value {
        ConstructibilityHypothesisV1::LocallyConstantOnStrata { witness } => {
            out.push(0);
            push_digest(out, witness);
        }
        ConstructibilityHypothesisV1::Controlled { witness } => {
            out.push(1);
            push_digest(out, witness);
        }
        ConstructibilityHypothesisV1::Unknown { no_claim } => {
            out.push(2);
            push_digest(out, no_claim);
        }
    }
}

fn push_properness(out: &mut Vec<u8>, value: ExitPropernessHypothesisV1) {
    match value {
        ExitPropernessHypothesisV1::Compact { witness } => {
            out.push(0);
            push_digest(out, witness);
        }
        ExitPropernessHypothesisV1::ProperLocallyFinite { witness } => {
            out.push(1);
            push_digest(out, witness);
        }
        ExitPropernessHypothesisV1::Unknown { no_claim } => {
            out.push(2);
            push_digest(out, no_claim);
        }
    }
}

fn push_refinement(out: &mut Vec<u8>, value: RefinementHypothesisV1) {
    match value {
        RefinementHypothesisV1::Identity { refinement } => {
            out.push(0);
            push_digest(out, refinement);
        }
        RefinementHypothesisV1::CommonRefinement {
            refinement,
            forward,
            reverse,
            witness,
        } => {
            out.push(1);
            push_digest(out, refinement);
            push_digest(out, forward);
            push_digest(out, reverse);
            push_digest(out, witness);
        }
        RefinementHypothesisV1::Unknown { no_claim } => {
            out.push(2);
            push_digest(out, no_claim);
        }
    }
}

fn push_homotopy(out: &mut Vec<u8>, value: HomotopyFidelityV1) {
    match value {
        HomotopyFidelityV1::IncidenceOnly => out.push(0),
        HomotopyFidelityV1::RetainedThrough {
            degree,
            coherence,
            witness,
        } => {
            out.push(1);
            out.push(degree);
            push_digest(out, coherence);
            push_digest(out, witness);
        }
        HomotopyFidelityV1::FullHigher { coherence, witness } => {
            out.push(2);
            push_digest(out, coherence);
            push_digest(out, witness);
        }
        HomotopyFidelityV1::Unknown { no_claim } => {
            out.push(3);
            push_digest(out, no_claim);
        }
    }
}

fn hypotheses_bytes(hypotheses: ExitPathHypothesesV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(512);
    push_links(&mut out, hypotheses.links);
    push_monodromy(&mut out, hypotheses.monodromy);
    push_constructibility(&mut out, hypotheses.constructibility);
    push_properness(&mut out, hypotheses.properness);
    push_refinement(&mut out, hypotheses.refinement);
    push_homotopy(&mut out, hypotheses.homotopy);
    out
}

fn truncation_family_bytes(max_truncation: u8) -> Vec<u8> {
    vec![
        max_truncation,
        1,
        u8::from(max_truncation >= 1),
        u8::from(max_truncation >= 2),
        1,
    ]
}

fn falsifier_bytes(falsifier: &ExitPathFalsifierV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(132);
    push_digest(&mut out, falsifier.id);
    out.push(match falsifier.kind {
        ExitPathFalsifierKindV1::SameIncidenceDifferentLink => 0,
        ExitPathFalsifierKindV1::SameIncidenceDifferentMonodromy => 1,
        ExitPathFalsifierKindV1::DirectionReversal => 2,
        ExitPathFalsifierKindV1::HypothesisDeletion => 3,
    });
    push_digest(&mut out, falsifier.left);
    push_digest(&mut out, falsifier.right);
    push_digest(&mut out, falsifier.witness);
    out
}

fn theorem_state_bytes(state: ExitPathTheoremStateV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(72);
    match state {
        ExitPathTheoremStateV1::Preregistered { card } => {
            out.push(0);
            push_digest(&mut out, card);
        }
        ExitPathTheoremStateV1::Candidate { card, witness } => {
            out.push(1);
            push_digest(&mut out, card);
            push_digest(&mut out, witness);
        }
        ExitPathTheoremStateV1::RefutationRecorded {
            approximation,
            falsifier,
        } => {
            out.push(2);
            push_approximation(&mut out, approximation);
            push_digest(&mut out, falsifier);
        }
        ExitPathTheoremStateV1::Unknown { no_claim } => {
            out.push(3);
            push_digest(&mut out, no_claim);
        }
    }
    out
}

fn push_approximation(out: &mut Vec<u8>, approximation: ExitPathApproximationV1) {
    match approximation {
        ExitPathApproximationV1::IncidencePoset => out.push(0),
        ExitPathApproximationV1::StratumGroupoidEnrichedExitCategory => out.push(1),
        ExitPathApproximationV1::SimplicialCategory {
            max_simplex_dimension,
        } => {
            out.push(2);
            out.push(max_simplex_dimension);
        }
        ExitPathApproximationV1::HigherTruncation { degree } => {
            out.push(3);
            out.push(degree);
        }
        ExitPathApproximationV1::FullHigherCategory => out.push(4),
    }
}

fn tcb_bytes(tcb: ExitPathTcbV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(96);
    push_digest(&mut out, tcb.tcb);
    push_digest(&mut out, tcb.checker);
    push_digest(&mut out, tcb.theorem_card);
    out
}

fn budget_bytes(budget: ExitPathBudgetV1) -> Vec<u8> {
    let mut out = Vec::with_capacity(32);
    out.push(budget.max_truncation);
    push_u64(&mut out, budget.max_referenced_artifact_slots);
    push_u64(&mut out, budget.max_implication_checks);
    push_f64(&mut out, budget.declared_wall_seconds);
    out
}

fn exit_path_receipt(
    ir: &ExitPathFamilyIrV1,
    cx: &Cx<'_>,
) -> Result<IdentityReceipt<ExitPathFamilySnapshotIdV1>, CanonicalError> {
    let subject = subject_bytes(ir);
    let conventions = conventions_bytes(ir);
    let hypotheses = hypotheses_bytes(ir.hypotheses);
    let truncation = truncation_family_bytes(ir.budget.max_truncation);
    let falsifiers: Vec<_> = ir.falsifiers.iter().map(falsifier_bytes).collect();
    let theorem_state = theorem_state_bytes(ir.theorem_state);
    let tcb = tcb_bytes(ir.tcb);
    let budget = budget_bytes(ir.budget);
    CanonicalEncoder::<ExitPathFamilySnapshotIdV1, _>::new(EXIT_PATH_IDENTITY_LIMITS, || {
        cx.is_cancel_requested()
    })?
    .bytes(Field::new(0, "subject"), &subject)?
    .bytes(Field::new(1, "conventions"), &conventions)?
    .bytes(Field::new(2, "hypotheses"), &hypotheses)?
    .bytes(Field::new(3, "truncation-family"), &truncation)?
    .canonical_set(
        Field::new(4, "falsifiers"),
        falsifiers.len() as u64,
        falsifiers.iter().map(Vec::as_slice),
    )?
    .bytes(Field::new(5, "theorem-state"), &theorem_state)?
    .bytes(Field::new(6, "tcb"), &tcb)?
    .bytes(Field::new(7, "budget"), &budget)?
    .finish()
}
