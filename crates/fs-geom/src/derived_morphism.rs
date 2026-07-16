//! Typed structural morphisms between admitted RD.1a geometries (RD.1b).
//!
//! This RD.1b spine admits category identities, generic strict maps, typed
//! declared chart maps, finite-complex rank refinements, and whole-object
//! inclusion declarations; checks structural evidence restriction/corestriction;
//! and composes ordered typed primitive paths with content-addressed lineage. A
//! separate stratum-scoped category admits component declarations only between
//! exact `(geometry, stratification, stratum)` objects and deliberately exposes
//! no whole-geometry evidence transport. A finite assembly candidate can bind
//! exactly one direct sealed component for every source stratum without granting
//! global-map authority. Another standalone token seals declared spans from two
//! admitted legs without folding correspondences into directed-map composition.
//! A final standalone token binds a fixed-resolution quasi-isomorphism
//! *candidate* to an exact refinement path and exact local presentations,
//! without granting theorem authority. This module deliberately cannot mint a
//! non-identity equivalence: a witness digest is data, not a proof of an inverse,
//! quasi-isomorphism, refinement theorem, or physical crosswalk.

use core::fmt;

use fs_blake3::identity::{
    CanonicalEncoder, CanonicalError, CanonicalLimits, CanonicalSchema, EvidenceNodeId, Field,
    FieldSpec, IdentityReceipt, StrongIdentity, WireType,
};
use fs_evidence::ColorRank;
use fs_exec::Cx;

use crate::derived::{
    AdmittedDerivedGeometryV1, CoefficientSystemV1, ConfigurationChartIdV1, ConfigurationChartV1,
    DerivedCheckerIdV1, DerivedComplexIdV1, DerivedComplexRoleV1, DerivedFrameIdV1,
    DerivedGeometryIdV1, DerivedLocalModelIdV1, DerivedLocalModelV1, DerivedModelVersionIdV1,
    DerivedNoClaimIdV1, DerivedResolutionIdV1, DerivedSubjectIdV1, DerivedTheoremIdV1,
    DerivedUnitSystemIdV1, DerivedWitnessIdV1, FiniteDerivedComplexV1, GeometricCategoryV1,
    PresentationScopeV1, StratificationIdV1, StratumIdV1,
};

/// Current schema for structural RD.1b morphism receipts.
pub const DERIVED_MORPHISM_SCHEMA_VERSION_V1: u32 = 1;
/// Current schema for standalone stratum-scoped morphism receipts.
pub const DERIVED_STRATUM_MORPHISM_SCHEMA_VERSION_V1: u32 = 1;
/// Current schema for exhaustive finite stratified-map assembly candidates.
pub const DERIVED_EXHAUSTIVE_STRATIFIED_MAP_CANDIDATE_SCHEMA_VERSION_V1: u32 = 1;
/// Current schema for standalone declared span-correspondence receipts.
pub const DERIVED_SPAN_CORRESPONDENCE_SCHEMA_VERSION_V1: u32 = 1;
/// Current schema for fixed-resolution quasi-isomorphism candidate receipts.
pub const DERIVED_FIXED_RESOLUTION_QUASI_ISOMORPHISM_CANDIDATE_SCHEMA_VERSION_V1: u32 = 1;
/// Maximum primitive nonidentity factors retained in one flattened composition.
pub const DERIVED_MORPHISM_MAX_FACTORS_V1: usize = 1024;
const DERIVED_MORPHISM_CANCELLATION_STRIDE_V1: usize = 64;
const DERIVED_MORPHISM_IDENTITY_LIMITS_V1: CanonicalLimits =
    CanonicalLimits::new(1 << 17, 1 << 16, 8, 1 << 11, 4096);
const DERIVED_STRATUM_MORPHISM_IDENTITY_LIMITS_V1: CanonicalLimits =
    CanonicalLimits::new(1 << 17, 1 << 16, 9, 1 << 11, 4096);
const DERIVED_EXHAUSTIVE_STRATIFIED_MAP_CANDIDATE_IDENTITY_LIMITS_V1: CanonicalLimits =
    CanonicalLimits::new(1 << 17, 1 << 16, 8, 1 << 11, 4096);
const DERIVED_FIXED_RESOLUTION_QUASI_ISOMORPHISM_CANDIDATE_IDENTITY_LIMITS_V1: CanonicalLimits =
    CanonicalLimits::new(1 << 17, 1 << 16, 16, 1 << 11, 4096);

/// Domain-separated semantic identity for one admitted structural morphism.
pub enum DerivedMorphismIdentitySchemaV1 {}

impl CanonicalSchema for DerivedMorphismIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-geom.derived-morphism.v1";
    const NAME: &'static str = "derived-geometry-structural-morphism";
    const VERSION: u32 = DERIVED_MORPHISM_SCHEMA_VERSION_V1;
    // Frozen v1 identity material: change only with a schema-version bump.
    const CONTEXT: &'static str = "typed endpoints, strict map class, evidence variance, no-equivalence boundary, and ordered primitive lineage";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("source", WireType::Bytes),
        FieldSpec::required("target", WireType::Bytes),
        FieldSpec::required("class", WireType::Bytes),
        FieldSpec::required("evidence-transport", WireType::Bytes),
        FieldSpec::required("no-equivalence-claims", WireType::OrderedBytes),
        FieldSpec::required("primitive-lineage", WireType::OrderedBytes),
    ];
}

/// Typed identity of one admitted RD.1b structural morphism.
pub type DerivedMorphismIdV1 = EvidenceNodeId<DerivedMorphismIdentitySchemaV1>;

/// Domain-separated semantic identity for one admitted stratum-scoped morphism.
pub enum DerivedStratumMorphismIdentitySchemaV1 {}

impl CanonicalSchema for DerivedStratumMorphismIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-geom.derived-stratum-morphism.v1";
    const NAME: &'static str = "derived-geometry-stratum-scoped-morphism";
    const VERSION: u32 = DERIVED_STRATUM_MORPHISM_SCHEMA_VERSION_V1;
    const CONTEXT: &'static str = "exact geometry, stratification, and stratum endpoints, nominal component declarations, no-authority boundary, and ordered primitive lineage";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("source-geometry", WireType::Bytes),
        FieldSpec::required("source-stratification", WireType::Bytes),
        FieldSpec::required("source-stratum", WireType::Bytes),
        FieldSpec::required("target-geometry", WireType::Bytes),
        FieldSpec::required("target-stratification", WireType::Bytes),
        FieldSpec::required("target-stratum", WireType::Bytes),
        FieldSpec::required("class", WireType::Bytes),
        FieldSpec::required("no-authority-claims", WireType::OrderedBytes),
        FieldSpec::required("primitive-lineage", WireType::OrderedBytes),
    ];
}

/// Typed identity of one admitted stratum-scoped morphism.
pub type DerivedStratumMorphismIdV1 = EvidenceNodeId<DerivedStratumMorphismIdentitySchemaV1>;

/// Domain-separated identity for one exhaustive finite component assembly.
pub enum DerivedExhaustiveStratifiedMapCandidateIdentitySchemaV1 {}

impl CanonicalSchema for DerivedExhaustiveStratifiedMapCandidateIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-geom.exhaustive-stratified-map-candidate.v1";
    const NAME: &'static str = "exhaustive-finite-stratified-map-assembly-candidate";
    const VERSION: u32 = DERIVED_EXHAUSTIVE_STRATIFIED_MAP_CANDIDATE_SCHEMA_VERSION_V1;
    const CONTEXT: &'static str = "exact endpoint geometries and stratifications, one explicit direct sealed component binding per source stratum in canonical order, nominal assembly and global constructibility declarations, and an explicit no-authority boundary";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("source-geometry", WireType::Bytes),
        FieldSpec::required("source-stratification", WireType::Bytes),
        FieldSpec::required("target-geometry", WireType::Bytes),
        FieldSpec::required("target-stratification", WireType::Bytes),
        FieldSpec::required("components", WireType::OrderedBytes),
        FieldSpec::required("nominal-assembly", WireType::Bytes),
        FieldSpec::required("nominal-constructibility", WireType::Bytes),
        FieldSpec::required("no-authority", WireType::Bytes),
    ];
}

/// Typed identity of one exhaustive finite stratified-map assembly candidate.
pub type DerivedExhaustiveStratifiedMapCandidateIdV1 =
    EvidenceNodeId<DerivedExhaustiveStratifiedMapCandidateIdentitySchemaV1>;

/// Domain-separated semantic identity for one admitted declared span.
pub enum DerivedSpanCorrespondenceIdentitySchemaV1 {}

impl CanonicalSchema for DerivedSpanCorrespondenceIdentitySchemaV1 {
    const DOMAIN: &'static str = "org.frankensim.fs-geom.derived-span-correspondence.v1";
    const NAME: &'static str = "derived-geometry-declared-span-correspondence";
    const VERSION: u32 = DERIVED_SPAN_CORRESPONDENCE_SCHEMA_VERSION_V1;
    const CONTEXT: &'static str =
        "exact source, common apex, exact target, ordered admitted legs, and no-authority boundary";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("source", WireType::Bytes),
        FieldSpec::required("apex", WireType::Bytes),
        FieldSpec::required("target", WireType::Bytes),
        FieldSpec::required("left-leg", WireType::Bytes),
        FieldSpec::required("right-leg", WireType::Bytes),
        FieldSpec::required("no-claim", WireType::Bytes),
    ];
}

/// Typed identity of one admitted standalone span `source <- apex -> target`.
pub type DerivedSpanCorrespondenceIdV1 = EvidenceNodeId<DerivedSpanCorrespondenceIdentitySchemaV1>;

/// Domain-separated semantic identity for one structurally admitted
/// fixed-resolution quasi-isomorphism candidate.
pub enum DerivedFixedResolutionQuasiIsomorphismCandidateIdentitySchemaV1 {}

impl CanonicalSchema for DerivedFixedResolutionQuasiIsomorphismCandidateIdentitySchemaV1 {
    const DOMAIN: &'static str =
        "org.frankensim.fs-geom.fixed-resolution-quasi-isomorphism-candidate.v1";
    const NAME: &'static str = "fixed-resolution-quasi-isomorphism-candidate";
    const VERSION: u32 = DERIVED_FIXED_RESOLUTION_QUASI_ISOMORPHISM_CANDIDATE_SCHEMA_VERSION_V1;
    const CONTEXT: &'static str = "exact refinement path, exact endpoint local models, complexes, resolutions, and roles, retained fixed-resolution scope witnesses, nominal external check metadata, and no-authority boundary";
    const FIELDS: &'static [FieldSpec] = &[
        FieldSpec::required("refinement-path", WireType::Bytes),
        FieldSpec::required("source-geometry", WireType::Bytes),
        FieldSpec::required("target-geometry", WireType::Bytes),
        FieldSpec::required("source-local-model", WireType::Bytes),
        FieldSpec::required("target-local-model", WireType::Bytes),
        FieldSpec::required("complex-role", WireType::Bytes),
        FieldSpec::required("source-complex", WireType::Bytes),
        FieldSpec::required("target-complex", WireType::Bytes),
        FieldSpec::required("source-resolution", WireType::Bytes),
        FieldSpec::required("target-resolution", WireType::Bytes),
        FieldSpec::required("source-scope-witness", WireType::Bytes),
        FieldSpec::required("target-scope-witness", WireType::Bytes),
        FieldSpec::required("nominal-theorem", WireType::Bytes),
        FieldSpec::required("nominal-checker", WireType::Bytes),
        FieldSpec::required("nominal-check-receipt", WireType::Bytes),
        FieldSpec::required("no-authority", WireType::Bytes),
    ];
}

/// Typed identity of one fixed-resolution quasi-isomorphism candidate.
pub type DerivedFixedResolutionQuasiIsomorphismCandidateIdV1 =
    EvidenceNodeId<DerivedFixedResolutionQuasiIsomorphismCandidateIdentitySchemaV1>;

/// Nominal artifact implementing one declared chart map.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DerivedChartMapIdV1([u8; 32]);

impl DerivedChartMapIdV1 {
    /// Construct a nominal chart-map artifact identity from exact bytes.
    #[must_use]
    pub const fn from_bytes(bytes: [u8; 32]) -> Self {
        Self(bytes)
    }

    /// Borrow the exact identity bytes.
    #[must_use]
    pub const fn as_bytes(&self) -> &[u8; 32] {
        &self.0
    }
}

/// Nominal scope artifact for the overlap on which a chart map is declared.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DerivedChartOverlapIdV1([u8; 32]);

impl DerivedChartOverlapIdV1 {
    /// Construct a nominal overlap-scope identity from exact bytes.
    #[must_use]
    pub const fn from_bytes(bytes: [u8; 32]) -> Self {
        Self(bytes)
    }

    /// Borrow the exact identity bytes.
    #[must_use]
    pub const fn as_bytes(&self) -> &[u8; 32] {
        &self.0
    }
}

/// Nominal aggregate prolongation artifact for one finite-complex refinement.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DerivedComplexRefinementMapIdV1([u8; 32]);

impl DerivedComplexRefinementMapIdV1 {
    /// Construct a nominal refinement-map artifact identity from exact bytes.
    #[must_use]
    pub const fn from_bytes(bytes: [u8; 32]) -> Self {
        Self(bytes)
    }

    /// Borrow the exact identity bytes.
    #[must_use]
    pub const fn as_bytes(&self) -> &[u8; 32] {
        &self.0
    }
}

/// Nominal map-artifact identity for one declared whole-object inclusion.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DerivedInclusionMapIdV1([u8; 32]);

impl DerivedInclusionMapIdV1 {
    /// Construct a nominal inclusion-map artifact identity from exact bytes.
    #[must_use]
    pub const fn from_bytes(bytes: [u8; 32]) -> Self {
        Self(bytes)
    }

    /// Borrow the exact identity bytes.
    #[must_use]
    pub const fn as_bytes(&self) -> &[u8; 32] {
        &self.0
    }
}

/// Nominal map artifact for one declared stratum-map component.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DerivedStratumMapIdV1([u8; 32]);

impl DerivedStratumMapIdV1 {
    /// Construct a nominal stratum-map identity from exact bytes.
    #[must_use]
    pub const fn from_bytes(bytes: [u8; 32]) -> Self {
        Self(bytes)
    }

    /// Borrow the exact identity bytes.
    #[must_use]
    pub const fn as_bytes(&self) -> &[u8; 32] {
        &self.0
    }
}

/// Nominal constructibility declaration for one stratum-map component.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DerivedConstructibilityDeclarationIdV1([u8; 32]);

impl DerivedConstructibilityDeclarationIdV1 {
    /// Construct a nominal declaration identity from exact bytes.
    #[must_use]
    pub const fn from_bytes(bytes: [u8; 32]) -> Self {
        Self(bytes)
    }

    /// Borrow the exact identity bytes.
    #[must_use]
    pub const fn as_bytes(&self) -> &[u8; 32] {
        &self.0
    }
}

/// Nominal whole-family artifact for one finite stratified-map assembly.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DerivedStratifiedMapAssemblyIdV1([u8; 32]);

impl DerivedStratifiedMapAssemblyIdV1 {
    /// Construct a nominal assembly identity from exact bytes.
    #[must_use]
    pub const fn from_bytes(bytes: [u8; 32]) -> Self {
        Self(bytes)
    }

    /// Borrow the exact identity bytes.
    #[must_use]
    pub const fn as_bytes(&self) -> &[u8; 32] {
        &self.0
    }
}

/// Nominal declaration of global constructibility for one component family.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DerivedGlobalConstructibilityDeclarationIdV1([u8; 32]);

impl DerivedGlobalConstructibilityDeclarationIdV1 {
    /// Construct a nominal global declaration identity from exact bytes.
    #[must_use]
    pub const fn from_bytes(bytes: [u8; 32]) -> Self {
        Self(bytes)
    }

    /// Borrow the exact identity bytes.
    #[must_use]
    pub const fn as_bytes(&self) -> &[u8; 32] {
        &self.0
    }
}

/// One exact object in the standalone stratum-scoped category.
///
/// The geometry alone is not the object: both the finite stratification and
/// selected stratum are identity-bearing parts of the endpoint.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DerivedStratumObjectV1 {
    /// Exact admitted geometry that owns the stratification.
    pub geometry: DerivedGeometryIdV1,
    /// Exact finite stratification owned by `geometry`.
    pub stratification: StratificationIdV1,
    /// Exact stratum owned by `stratification`.
    pub stratum: StratumIdV1,
}

/// Caller-supplied primitive in the standalone stratum-scoped category.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DerivedStratumMorphismKindV1 {
    /// The categorical identity on one exact stratum object.
    Identity,
    /// A nominal component-map declaration with no whole-map authority.
    DeclaredComponent {
        /// Nominal implementation or mathematical map artifact.
        map: DerivedStratumMapIdV1,
        /// Nominal declaration that this component is constructible.
        constructibility: DerivedConstructibilityDeclarationIdV1,
    },
}

/// Explicit authority boundary for one stratum-scoped request.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DerivedStratumAuthorityBoundaryV1 {
    /// Available only for the exact identity on one stratum object.
    IdentityOnly,
    /// No execution, constructibility, whole-map, or evidence authority claimed.
    NoClaim {
        /// Retained no-authority artifact.
        artifact: DerivedNoClaimIdV1,
    },
}

/// Versioned primitive request in the standalone stratum-scoped category.
///
/// This IR has no whole-geometry evidence field and cannot be converted into
/// `DerivedMorphismIrV1` by this module.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DerivedStratumMorphismIrV1 {
    /// Decoded schema version.
    pub schema_version: u32,
    /// Exact source stratum object.
    pub source: DerivedStratumObjectV1,
    /// Exact target stratum object.
    pub target: DerivedStratumObjectV1,
    /// Identity or nominal component declaration.
    pub kind: DerivedStratumMorphismKindV1,
    /// Honest authority boundary.
    pub authority: DerivedStratumAuthorityBoundaryV1,
}

/// Admitted family in the standalone stratum-scoped category.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdmittedDerivedStratumMorphismClassV1 {
    /// Exact identity on one `(geometry, stratification, stratum)` object.
    Identity,
    /// One or more ordered declared component primitives.
    DeclaredPath,
}

/// One retained nominal stratum-map component.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DeclaredStratumMapPrimitiveV1 {
    /// Exact source object of this component.
    pub source: DerivedStratumObjectV1,
    /// Exact target object of this component.
    pub target: DerivedStratumObjectV1,
    /// Nominal component-map artifact.
    pub map: DerivedStratumMapIdV1,
    /// Nominal constructibility declaration; not authenticated here.
    pub constructibility: DerivedConstructibilityDeclarationIdV1,
}

/// One explicit source/target selector binding in a finite map assembly.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DerivedStratifiedMapComponentBindingV1 {
    /// Source stratum covered by this exact component.
    pub source_stratum: StratumIdV1,
    /// Target stratum selected by this exact component.
    pub target_stratum: StratumIdV1,
    /// Exact sealed direct component receipt.
    pub component: DerivedStratumMorphismIdV1,
}

/// Versioned exhaustive finite stratified-map assembly candidate.
///
/// `components` must follow the source geometry's canonical stratum order.
/// Each ID names one separately sealed direct stratum component. The candidate
/// establishes finite source coverage only; it is not a whole-geometry map.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DerivedExhaustiveStratifiedMapCandidateIrV1 {
    /// Decoded schema version.
    pub schema_version: u32,
    /// Exact admitted source geometry.
    pub source_geometry: DerivedGeometryIdV1,
    /// Exact finite source stratification.
    pub source_stratification: StratificationIdV1,
    /// Exact admitted target geometry.
    pub target_geometry: DerivedGeometryIdV1,
    /// Exact finite target stratification.
    pub target_stratification: StratificationIdV1,
    /// One direct sealed component binding per source stratum, in canonical order.
    pub components: Vec<DerivedStratifiedMapComponentBindingV1>,
    /// Nominal whole-family/assembly artifact; not resolved by RD.1b.
    pub nominal_assembly: DerivedStratifiedMapAssemblyIdV1,
    /// Nominal global constructibility declaration; not authenticated by RD.1b.
    pub nominal_constructibility: DerivedGlobalConstructibilityDeclarationIdV1,
    /// Explicit denial of global map, gluing, and theorem authority.
    pub no_authority: DerivedNoClaimIdV1,
}

/// Caller-supplied primitive map class.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DerivedMorphismKindV1 {
    /// The categorical identity on one exact admitted geometry.
    Identity,
    /// A strict map whose construction is retained but not theorem-promoted.
    Strict {
        /// Exact map/construction artifact.
        witness: DerivedWitnessIdV1,
    },
    /// A declared coordinate map on one nominal overlap, without invertibility.
    DeclaredChartMap {
        /// Chart owned by the source geometry.
        source_chart: ConfigurationChartIdV1,
        /// Chart owned by the target geometry.
        target_chart: ConfigurationChartIdV1,
        /// Exact declared overlap/scope artifact.
        overlap: DerivedChartOverlapIdV1,
        /// Exact forward coordinate-map artifact.
        map: DerivedChartMapIdV1,
    },
    /// A coarse-to-refined finite-complex declaration without chain-map authority.
    DeclaredComplexRefinement {
        /// Complex owned by the coarse/source geometry.
        source_complex: DerivedComplexIdV1,
        /// Complex owned by the refined/target geometry.
        target_complex: DerivedComplexIdV1,
        /// Exact resolution retained by the source complex.
        source_resolution: DerivedResolutionIdV1,
        /// Exact resolution retained by the target complex.
        target_resolution: DerivedResolutionIdV1,
        /// Nominal aggregate coarse-to-refined prolongation artifact.
        prolongation: DerivedComplexRefinementMapIdV1,
        /// Nominal differential-commutation witness; not authenticated here.
        commutation: DerivedWitnessIdV1,
    },
    /// A whole-object inclusion declaration without containment authority.
    DeclaredInclusion {
        /// Nominal declared source-to-target map artifact.
        map: DerivedInclusionMapIdV1,
        /// Nominal claimed-containment artifact; not authenticated here.
        containment: DerivedWitnessIdV1,
    },
}

/// Admitted map family. Composition flattens ordered typed primitive lineage.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdmittedDerivedMorphismClassV1 {
    /// Exact categorical identity.
    Identity,
    /// One or more ordered strict primitive maps.
    Strict,
    /// One or more ordered declared chart-map primitives.
    DeclaredChartMapPath,
    /// One or more ordered finite-complex refinement primitives.
    DeclaredComplexRefinementPath,
    /// One or more ordered whole-object inclusion declarations.
    DeclaredInclusionPath,
    /// An ordered path containing more than one nonidentity primitive family.
    HeterogeneousPath,
}

/// Exact chart endpoints retained for a homogeneous declared chart-map path.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DerivedChartPathV1 {
    /// First source chart in the path.
    pub source_chart: ConfigurationChartIdV1,
    /// Last target chart in the path.
    pub target_chart: ConfigurationChartIdV1,
}

/// One retained declared chart-map primitive.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DeclaredChartMapPrimitiveV1 {
    /// Exact admitted geometry owning the source chart.
    pub source_geometry: DerivedGeometryIdV1,
    /// Exact admitted geometry owning the target chart.
    pub target_geometry: DerivedGeometryIdV1,
    /// Chart owned by the primitive source geometry.
    pub source_chart: ConfigurationChartIdV1,
    /// Chart owned by the primitive target geometry.
    pub target_chart: ConfigurationChartIdV1,
    /// Nominal overlap/scope artifact.
    pub overlap: DerivedChartOverlapIdV1,
    /// Nominal forward coordinate-map artifact.
    pub map: DerivedChartMapIdV1,
}

/// One retained coarse-to-refined finite-complex declaration.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DeclaredComplexRefinementPrimitiveV1 {
    /// Exact admitted geometry owning the coarse complex.
    pub source_geometry: DerivedGeometryIdV1,
    /// Exact admitted geometry owning the refined complex.
    pub target_geometry: DerivedGeometryIdV1,
    /// Exact coarse/source complex.
    pub source_complex: DerivedComplexIdV1,
    /// Exact refined/target complex.
    pub target_complex: DerivedComplexIdV1,
    /// Exact source resolution.
    pub source_resolution: DerivedResolutionIdV1,
    /// Exact target resolution.
    pub target_resolution: DerivedResolutionIdV1,
    /// Nominal aggregate prolongation artifact.
    pub prolongation: DerivedComplexRefinementMapIdV1,
    /// Nominal commutation witness with zero theorem authority in v1.
    pub commutation: DerivedWitnessIdV1,
}

/// One retained whole-object inclusion declaration.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DeclaredInclusionPrimitiveV1 {
    /// Exact admitted source geometry declared to be included.
    pub source_geometry: DerivedGeometryIdV1,
    /// Exact admitted target geometry declared to contain the source.
    pub target_geometry: DerivedGeometryIdV1,
    /// Nominal declared inclusion-map artifact.
    pub map: DerivedInclusionMapIdV1,
    /// Nominal claimed-containment artifact with zero theorem authority in v1.
    pub containment: DerivedWitnessIdV1,
}

/// One typed nonidentity primitive retained in semantic path order.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdmittedDerivedPrimitiveV1 {
    /// One generic strict primitive.
    Strict {
        /// Exact admitted source geometry.
        source_geometry: DerivedGeometryIdV1,
        /// Exact admitted target geometry.
        target_geometry: DerivedGeometryIdV1,
        /// Exact map/construction artifact.
        witness: DerivedWitnessIdV1,
    },
    /// One declared chart-map primitive.
    DeclaredChartMap(DeclaredChartMapPrimitiveV1),
    /// One declared finite-complex refinement primitive.
    DeclaredComplexRefinement(DeclaredComplexRefinementPrimitiveV1),
    /// One declared whole-object inclusion primitive.
    DeclaredInclusion(DeclaredInclusionPrimitiveV1),
}

impl AdmittedDerivedPrimitiveV1 {
    /// Exact admitted source geometry of this primitive.
    #[must_use]
    pub const fn source_geometry(&self) -> DerivedGeometryIdV1 {
        match self {
            Self::Strict {
                source_geometry, ..
            } => *source_geometry,
            Self::DeclaredChartMap(primitive) => primitive.source_geometry,
            Self::DeclaredComplexRefinement(primitive) => primitive.source_geometry,
            Self::DeclaredInclusion(primitive) => primitive.source_geometry,
        }
    }

    /// Exact admitted target geometry of this primitive.
    #[must_use]
    pub const fn target_geometry(&self) -> DerivedGeometryIdV1 {
        match self {
            Self::Strict {
                target_geometry, ..
            } => *target_geometry,
            Self::DeclaredChartMap(primitive) => primitive.target_geometry,
            Self::DeclaredComplexRefinement(primitive) => primitive.target_geometry,
            Self::DeclaredInclusion(primitive) => primitive.target_geometry,
        }
    }
}

/// Direction in which evidence is transported along an object map `X -> Y`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DerivedEvidenceVarianceV1 {
    /// Identity evidence on one exact object.
    Identity,
    /// Sheaf-like restriction `E(Y) -> E(X)`.
    RestrictionContravariant,
    /// Cosheaf/balance-like corestriction `B(X) -> B(Y)`.
    BalanceCorestrictionCovariant,
}

/// Nominal caller-declared identity of an evidence artifact.
///
/// This identity binds structural transport lineage into a morphism receipt.
/// It does not authenticate a payload or grant evidence authority.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct DerivedEvidenceArtifactIdV1([u8; 32]);

impl DerivedEvidenceArtifactIdV1 {
    /// Construct a nominal evidence-artifact identity from exact bytes.
    #[must_use]
    pub const fn from_bytes(bytes: [u8; 32]) -> Self {
        Self(bytes)
    }

    /// Borrow the exact identity bytes.
    #[must_use]
    pub const fn as_bytes(&self) -> &[u8; 32] {
        &self.0
    }
}

/// Structural evidence transport attached to one object map.
///
/// Identity arrows are deliberately rank-neutral, so every admitted geometry
/// has one categorical identity. Nonidentity variants bind caller-declared
/// evidence artifact identities and ranks, but do not authenticate either
/// artifact.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DerivedEvidenceTransportV1 {
    /// Unique rank-neutral evidence transport on one identity arrow.
    Identity,
    /// Sheaf-like restriction `E(Y) -> E(X)` along `X -> Y`.
    RestrictionContravariant {
        /// Geometry from which evidence is consumed.
        input_geometry: DerivedGeometryIdV1,
        /// Geometry at which transported evidence is published.
        output_geometry: DerivedGeometryIdV1,
        /// Nominal consumed evidence artifact.
        input_evidence: DerivedEvidenceArtifactIdV1,
        /// Nominal published evidence artifact.
        output_evidence: DerivedEvidenceArtifactIdV1,
        /// Caller-declared authority rank consumed by the transport.
        input_rank: ColorRank,
        /// Caller-declared authority rank published by the transport.
        output_rank: ColorRank,
    },
    /// Cosheaf/balance-like corestriction `B(X) -> B(Y)` along `X -> Y`.
    BalanceCorestrictionCovariant {
        /// Geometry from which evidence is consumed.
        input_geometry: DerivedGeometryIdV1,
        /// Geometry at which transported evidence is published.
        output_geometry: DerivedGeometryIdV1,
        /// Nominal consumed evidence artifact.
        input_evidence: DerivedEvidenceArtifactIdV1,
        /// Nominal published evidence artifact.
        output_evidence: DerivedEvidenceArtifactIdV1,
        /// Caller-declared authority rank consumed by the transport.
        input_rank: ColorRank,
        /// Caller-declared authority rank published by the transport.
        output_rank: ColorRank,
    },
}

impl DerivedEvidenceTransportV1 {
    /// Variance of this structural evidence transport.
    #[must_use]
    pub const fn variance(self) -> DerivedEvidenceVarianceV1 {
        match self {
            Self::Identity => DerivedEvidenceVarianceV1::Identity,
            Self::RestrictionContravariant { .. } => {
                DerivedEvidenceVarianceV1::RestrictionContravariant
            }
            Self::BalanceCorestrictionCovariant { .. } => {
                DerivedEvidenceVarianceV1::BalanceCorestrictionCovariant
            }
        }
    }

    /// Caller-declared input/output ranks for a strict transport.
    #[must_use]
    pub const fn ranks(self) -> Option<(ColorRank, ColorRank)> {
        match self {
            Self::Identity => None,
            Self::RestrictionContravariant {
                input_rank,
                output_rank,
                ..
            }
            | Self::BalanceCorestrictionCovariant {
                input_rank,
                output_rank,
                ..
            } => Some((input_rank, output_rank)),
        }
    }
}

/// Explicit equivalence boundary carried by a primitive map request.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DerivedEquivalenceBoundaryV1 {
    /// Available only for an exact identity arrow.
    IdentityOnly,
    /// No inverse, equivalence, or physical correspondence is claimed.
    NoClaim {
        /// Retained no-claim artifact.
        artifact: DerivedNoClaimIdV1,
    },
}

/// Versioned primitive RD.1b morphism request.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DerivedMorphismIrV1 {
    /// Decoded schema version.
    pub schema_version: u32,
    /// Exact admitted source geometry.
    pub source: DerivedGeometryIdV1,
    /// Exact admitted target geometry.
    pub target: DerivedGeometryIdV1,
    /// Primitive map class.
    pub kind: DerivedMorphismKindV1,
    /// Directional evidence transport.
    pub evidence: DerivedEvidenceTransportV1,
    /// Honest equivalence/no-equivalence boundary.
    pub equivalence: DerivedEquivalenceBoundaryV1,
}

/// Versioned standalone span request `source <- apex -> target`.
///
/// `left_leg` and `right_leg` name already-admitted morphisms with orientations
/// `apex -> source` and `apex -> target`. The request has no direct
/// source-to-target evidence transport and makes no pullback or functionality
/// claim.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DerivedSpanCorrespondenceIrV1 {
    /// Decoded schema version.
    pub schema_version: u32,
    /// Exact left endpoint.
    pub source: DerivedGeometryIdV1,
    /// Exact common source of both admitted legs.
    pub apex: DerivedGeometryIdV1,
    /// Exact right endpoint.
    pub target: DerivedGeometryIdV1,
    /// Exact admitted `apex -> source` morphism.
    pub left_leg: DerivedMorphismIdV1,
    /// Exact admitted `apex -> target` morphism.
    pub right_leg: DerivedMorphismIdV1,
    /// Explicit no-correspondence-authority artifact.
    pub no_claim: DerivedNoClaimIdV1,
}

/// Versioned declaration of one fixed-resolution quasi-isomorphism candidate.
///
/// The exact supplied refinement path is the only directed map. The theorem,
/// checker, and check-receipt IDs are nominal metadata retained for RD.1c; RD.1b
/// neither dereferences them nor promotes the candidate to an equivalence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DerivedFixedResolutionQuasiIsomorphismCandidateIrV1 {
    /// Decoded schema version.
    pub schema_version: u32,
    /// Exact admitted source geometry.
    pub source_geometry: DerivedGeometryIdV1,
    /// Exact admitted target geometry.
    pub target_geometry: DerivedGeometryIdV1,
    /// Exact source local presentation.
    pub source_local_model: DerivedLocalModelIdV1,
    /// Exact target local presentation.
    pub target_local_model: DerivedLocalModelIdV1,
    /// Tangent, cotangent, or deformation-obstruction complex under comparison.
    pub complex_role: DerivedComplexRoleV1,
    /// Exact complex selected by the source local model for `complex_role`.
    pub source_complex: DerivedComplexIdV1,
    /// Exact complex selected by the target local model for `complex_role`.
    pub target_complex: DerivedComplexIdV1,
    /// Exact source finite resolution.
    pub source_resolution: DerivedResolutionIdV1,
    /// Exact target finite resolution.
    pub target_resolution: DerivedResolutionIdV1,
    /// Exact sealed homogeneous refinement path used as the directed candidate.
    pub refinement_path: DerivedMorphismIdV1,
    /// Nominal theorem-card identity for the claimed cohomology isomorphism.
    pub nominal_theorem: DerivedTheoremIdV1,
    /// Nominal independent-checker identity; not executed by RD.1b.
    pub nominal_checker: DerivedCheckerIdV1,
    /// Nominal external check receipt; not authenticated by RD.1b.
    pub nominal_check_receipt: DerivedWitnessIdV1,
    /// Explicit artifact denying quasi-isomorphism, presentation-equivalence,
    /// and physical-equivalence authority to this structural token.
    pub no_authority: DerivedNoClaimIdV1,
}

/// Structured refusal from RD.1b structural morphism admission/composition.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DerivedMorphismErrorV1 {
    /// Unsupported decoded schema version.
    UnsupportedSchemaVersion {
        /// Supplied version.
        found: u32,
        /// Sole supported version.
        supported: u32,
    },
    /// The raw request did not name the exact supplied admitted endpoint.
    EndpointMismatch {
        /// Stable endpoint field.
        field: &'static str,
    },
    /// A required nominal identity used the all-zero sentinel.
    MissingIdentity {
        /// Stable identity field.
        field: &'static str,
    },
    /// An identity arrow changed its object, transport, or boundary.
    InvalidIdentity,
    /// Nonidentity source and target describe different physical subjects.
    SubjectMismatch,
    /// Nonidentity endpoints name different immutable model versions.
    ModelVersionMismatch,
    /// Nonidentity endpoints use different mathematical categories.
    CategoryMismatch,
    /// Nonidentity endpoints use different coefficient semantics.
    CoefficientMismatch,
    /// Nonidentity endpoints use different coordinate frames.
    FrameMismatch,
    /// Nonidentity endpoints use different unit systems.
    UnitSystemMismatch,
    /// A declared chart ID is not owned by its exact endpoint geometry.
    MissingChart {
        /// Stable source/target chart field.
        field: &'static str,
    },
    /// Declared chart-map coordinate or ambient dimensions differ.
    ChartDimensionMismatch,
    /// Declared chart-map frames or coordinate-unit bindings differ.
    ChartConventionMismatch,
    /// A declared finite-complex ID is not owned by its exact endpoint geometry.
    MissingComplex {
        /// Stable source/target complex field.
        field: &'static str,
    },
    /// A finite-complex refinement violates one structural rank-envelope rule.
    ComplexRefinementMismatch {
        /// Stable failed refinement field or relation.
        field: &'static str,
    },
    /// Evidence endpoints contradict the declared variance.
    EvidenceOrientationMismatch,
    /// Declared evidence metadata attempted to strengthen authority.
    EvidenceStrengthening {
        /// Consumed authority.
        input: ColorRank,
        /// Attempted published authority.
        output: ColorRank,
    },
    /// A non-identity map attempted to claim identity/equivalence authority.
    EquivalenceLaundering,
    /// Two arrows do not share the exact middle object.
    CompositionEndpointMismatch,
    /// Nonidentity paths use incompatible evidence variance.
    CompositionVarianceMismatch,
    /// A sealed primitive path is missing or internally inconsistent.
    CompositionClassMismatch,
    /// Adjacent chart-map primitives do not share the exact middle chart.
    CompositionChartMismatch,
    /// Adjacent refinement primitives do not share an exact complex/resolution seam.
    CompositionRefinementMismatch,
    /// Evidence artifact identity or declared rank at the seam is inconsistent.
    CompositionEvidenceMismatch,
    /// Flattened lineage/no-claim retention exceeded the hard ceiling.
    ResourceLimit {
        /// Stable retained collection.
        field: &'static str,
        /// Requested entries.
        requested: usize,
        /// Hard limit.
        limit: usize,
    },
    /// A fallible lineage/no-claim allocation was refused.
    AllocationRefused {
        /// Stable retained collection.
        field: &'static str,
    },
    /// Cooperative cancellation was observed before publication.
    Cancelled {
        /// Stable admission/composition stage.
        stage: &'static str,
    },
    /// Canonical identity construction failed.
    Identity(CanonicalError),
}

impl fmt::Display for DerivedMorphismErrorV1 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "derived geometry morphism refused: {self:?}")
    }
}

impl core::error::Error for DerivedMorphismErrorV1 {}

/// Structured refusal from standalone stratum-scoped admission/composition.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DerivedStratumMorphismErrorV1 {
    /// Unsupported decoded schema version.
    UnsupportedSchemaVersion {
        /// Supplied version.
        found: u32,
        /// Sole supported version.
        supported: u32,
    },
    /// A raw stratum object did not name the exact supplied admitted geometry.
    EndpointMismatch {
        /// Stable source/target geometry field.
        field: &'static str,
    },
    /// A required nominal identity used the all-zero sentinel.
    MissingIdentity {
        /// Stable identity field.
        field: &'static str,
    },
    /// A stratification selector is not the exact one owned by its geometry.
    MissingStratification {
        /// Stable source/target stratification field.
        field: &'static str,
    },
    /// A stratum selector is not owned by its exact finite stratification.
    MissingStratum {
        /// Stable source/target stratum field.
        field: &'static str,
    },
    /// An identity request changed its exact object or authority boundary.
    InvalidIdentity,
    /// Nonidentity endpoints describe different physical subjects.
    SubjectMismatch,
    /// Nonidentity endpoints name different immutable model versions.
    ModelVersionMismatch,
    /// Nonidentity endpoints use different mathematical categories.
    CategoryMismatch,
    /// Nonidentity endpoints use different coefficient semantics.
    CoefficientMismatch,
    /// Nonidentity endpoints use different coordinate frames.
    FrameMismatch,
    /// Nonidentity endpoints use different unit systems.
    UnitSystemMismatch,
    /// A component declaration attempted to claim identity authority.
    AuthorityLaundering,
    /// Two arrows do not share the exact middle stratum object.
    CompositionEndpointMismatch,
    /// A sealed token has an internally inconsistent identity/path class.
    CompositionClassMismatch,
    /// Flattened lineage/no-claim retention exceeded the hard ceiling.
    ResourceLimit {
        /// Stable retained collection.
        field: &'static str,
        /// Requested entries.
        requested: usize,
        /// Hard limit.
        limit: usize,
    },
    /// A fallible lineage/no-claim allocation was refused.
    AllocationRefused {
        /// Stable retained collection.
        field: &'static str,
    },
    /// Cooperative cancellation was observed before publication.
    Cancelled {
        /// Stable admission/composition stage.
        stage: &'static str,
    },
    /// Canonical identity construction failed.
    Identity(CanonicalError),
}

impl fmt::Display for DerivedStratumMorphismErrorV1 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "derived stratum morphism refused: {self:?}")
    }
}

impl core::error::Error for DerivedStratumMorphismErrorV1 {}

/// Structured refusal from exhaustive finite stratified-map assembly admission.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DerivedExhaustiveStratifiedMapCandidateErrorV1 {
    /// Unsupported decoded schema version.
    UnsupportedSchemaVersion {
        /// Supplied version.
        found: u32,
        /// Sole supported version.
        supported: u32,
    },
    /// A raw geometry selector did not name the exact supplied admitted endpoint.
    EndpointMismatch {
        /// Stable source/target geometry field.
        field: &'static str,
    },
    /// A required nominal identity used the all-zero sentinel.
    MissingIdentity {
        /// Stable identity field.
        field: &'static str,
    },
    /// A raw stratification selector is not owned by its exact geometry.
    StratificationMismatch {
        /// Stable source/target stratification field.
        field: &'static str,
    },
    /// Raw or supplied component count differs from exact source-stratum coverage.
    ComponentCountMismatch {
        /// Stable raw/supplied component collection.
        field: &'static str,
        /// Required entries.
        expected: usize,
        /// Supplied entries.
        found: usize,
    },
    /// A raw component ID does not name the supplied sealed component.
    ComponentIdentityMismatch {
        /// Canonical source-stratum index.
        index: usize,
    },
    /// A raw binding or sealed component contradicts an exact endpoint selector.
    ComponentEndpointMismatch {
        /// Canonical source-stratum index.
        index: usize,
        /// Stable failed endpoint relation.
        field: &'static str,
    },
    /// A supplied component is a composite path rather than one direct declaration.
    CompositeComponent {
        /// Canonical source-stratum index.
        index: usize,
    },
    /// Component retention exceeded the hard ceiling.
    ResourceLimit {
        /// Stable retained collection.
        field: &'static str,
        /// Requested entries.
        requested: usize,
        /// Hard limit.
        limit: usize,
    },
    /// A fallible component allocation was refused.
    AllocationRefused {
        /// Stable retained collection.
        field: &'static str,
    },
    /// Cooperative cancellation was observed before publication.
    Cancelled {
        /// Stable admission stage.
        stage: &'static str,
    },
    /// Canonical identity construction failed.
    Identity(CanonicalError),
}

impl fmt::Display for DerivedExhaustiveStratifiedMapCandidateErrorV1 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "exhaustive stratified-map candidate refused: {self:?}")
    }
}

impl core::error::Error for DerivedExhaustiveStratifiedMapCandidateErrorV1 {}

/// Structured refusal from standalone declared-span admission.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DerivedSpanCorrespondenceErrorV1 {
    /// Unsupported decoded schema version.
    UnsupportedSchemaVersion {
        /// Supplied version.
        found: u32,
        /// Sole supported version.
        supported: u32,
    },
    /// The explicit no-authority artifact used the all-zero sentinel.
    MissingIdentity {
        /// Stable identity field.
        field: &'static str,
    },
    /// A raw leg ID does not name the supplied sealed admitted leg.
    LegIdentityMismatch {
        /// Stable left/right leg field.
        field: &'static str,
    },
    /// One sealed leg does not have the required apex/endpoint orientation.
    LegOrientationMismatch {
        /// Stable failed orientation relation.
        field: &'static str,
    },
    /// Cooperative cancellation was observed before publication.
    Cancelled {
        /// Stable admission stage.
        stage: &'static str,
    },
    /// Canonical identity construction failed.
    Identity(CanonicalError),
}

impl fmt::Display for DerivedSpanCorrespondenceErrorV1 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "derived span correspondence refused: {self:?}")
    }
}

impl core::error::Error for DerivedSpanCorrespondenceErrorV1 {}

/// Structured refusal from fixed-resolution quasi-isomorphism candidate admission.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1 {
    /// Unsupported decoded schema version.
    UnsupportedSchemaVersion {
        /// Supplied version.
        found: u32,
        /// Sole supported version.
        supported: u32,
    },
    /// A raw endpoint does not name the supplied admitted geometry.
    EndpointMismatch {
        /// Stable source/target geometry field.
        field: &'static str,
    },
    /// A required opaque identity used the all-zero sentinel.
    MissingIdentity {
        /// Stable identity field.
        field: &'static str,
    },
    /// The raw refinement-path ID does not name the supplied sealed path.
    PathIdentityMismatch,
    /// The supplied sealed path is not a homogeneous complex-refinement path.
    PathClassMismatch {
        /// Actual sealed path family.
        found: AdmittedDerivedMorphismClassV1,
    },
    /// A retained refinement primitive contradicts the sealed path or candidate scope.
    PathShapeMismatch {
        /// Stable failed path relation.
        field: &'static str,
        /// Zero-based primitive index, or zero for an empty path.
        index: usize,
    },
    /// A selected local-model ID is not owned by its exact endpoint geometry.
    MissingLocalModel {
        /// Stable source/target local-model field.
        field: &'static str,
    },
    /// A selected complex ID is not owned by its exact endpoint geometry.
    MissingCandidateComplex {
        /// Stable source/target complex field.
        field: &'static str,
    },
    /// A selected local model does not bind the complex under the declared role.
    ComplexRoleMismatch {
        /// Stable failed source/target role relation.
        field: &'static str,
    },
    /// A selected local model and complex disagree on their exact chart.
    LocalPresentationMismatch {
        /// Stable failed source/target presentation relation.
        field: &'static str,
    },
    /// A complex, local-model scope, and raw resolution selector disagree.
    ResolutionScopeMismatch {
        /// Stable failed source/target resolution relation.
        field: &'static str,
    },
    /// A local model is not explicitly scoped to one fixed finite resolution.
    PresentationScopeMismatch {
        /// Stable source/target presentation field.
        field: &'static str,
    },
    /// The endpoint local presentations do not name one exact locality.
    LocalityMismatch,
    /// Cooperative cancellation was observed before publication.
    Cancelled {
        /// Stable admission stage.
        stage: &'static str,
    },
    /// Canonical identity construction failed.
    Identity(CanonicalError),
}

impl fmt::Display for DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1 {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "fixed-resolution quasi-isomorphism candidate refused: {self:?}"
        )
    }
}

impl core::error::Error for DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1 {}

/// Sealed, content-addressed structural morphism.
#[derive(Debug, PartialEq, Eq)]
pub struct AdmittedDerivedMorphismV1 {
    source: DerivedGeometryIdV1,
    target: DerivedGeometryIdV1,
    class: AdmittedDerivedMorphismClassV1,
    chart_path: Option<DerivedChartPathV1>,
    primitive_path: Vec<AdmittedDerivedPrimitiveV1>,
    declared_chart_maps: Vec<DeclaredChartMapPrimitiveV1>,
    evidence: DerivedEvidenceTransportV1,
    no_equivalence_claims: Vec<DerivedNoClaimIdV1>,
    primitive_factors: Vec<DerivedMorphismIdV1>,
    receipt: IdentityReceipt<DerivedMorphismIdV1>,
}

impl AdmittedDerivedMorphismV1 {
    /// Exact source geometry.
    #[must_use]
    pub const fn source(&self) -> DerivedGeometryIdV1 {
        self.source
    }

    /// Exact target geometry.
    #[must_use]
    pub const fn target(&self) -> DerivedGeometryIdV1 {
        self.target
    }

    /// Identity versus admitted nonidentity path family.
    #[must_use]
    pub const fn class(&self) -> AdmittedDerivedMorphismClassV1 {
        self.class
    }

    /// Chart endpoints for a homogeneous declared chart-map path.
    #[must_use]
    pub const fn chart_path(&self) -> Option<DerivedChartPathV1> {
        self.chart_path
    }

    /// Ordered typed nonidentity primitives, including their exact endpoints.
    #[must_use]
    pub fn primitive_path(&self) -> &[AdmittedDerivedPrimitiveV1] {
        &self.primitive_path
    }

    /// Ordered chart-map primitives retained from the full typed path.
    #[must_use]
    pub fn declared_chart_maps(&self) -> &[DeclaredChartMapPrimitiveV1] {
        &self.declared_chart_maps
    }

    /// Checked structural direction, nominal artifact identities, and declared ranks.
    #[must_use]
    pub const fn evidence(&self) -> DerivedEvidenceTransportV1 {
        self.evidence
    }

    /// Ordered primitive factors after associative flattening.
    #[must_use]
    pub fn primitive_factors(&self) -> &[DerivedMorphismIdV1] {
        &self.primitive_factors
    }

    /// Ordered retained no-equivalence artifacts.
    #[must_use]
    pub fn no_equivalence_claims(&self) -> &[DerivedNoClaimIdV1] {
        &self.no_equivalence_claims
    }

    /// Typed structural morphism identity.
    #[must_use]
    pub const fn id(&self) -> DerivedMorphismIdV1 {
        self.receipt.id()
    }

    /// Canonical receipt and construction limits.
    #[must_use]
    pub const fn identity_receipt(&self) -> IdentityReceipt<DerivedMorphismIdV1> {
        self.receipt
    }
}

/// Sealed morphism in the standalone stratum-scoped category.
///
/// Its endpoints are exact `(geometry, stratification, stratum)` objects. The
/// token is intentionally not an `AdmittedDerivedMorphismV1`, exposes no
/// geometry-wide evidence transport, and cannot establish a map on any
/// unlisted stratum.
#[derive(Debug, PartialEq, Eq)]
pub struct AdmittedDerivedStratumMorphismV1 {
    source: DerivedStratumObjectV1,
    target: DerivedStratumObjectV1,
    class: AdmittedDerivedStratumMorphismClassV1,
    primitive_path: Vec<DeclaredStratumMapPrimitiveV1>,
    no_authority_claims: Vec<DerivedNoClaimIdV1>,
    primitive_factors: Vec<DerivedStratumMorphismIdV1>,
    receipt: IdentityReceipt<DerivedStratumMorphismIdV1>,
}

impl AdmittedDerivedStratumMorphismV1 {
    /// Exact source stratum object.
    #[must_use]
    pub const fn source(&self) -> DerivedStratumObjectV1 {
        self.source
    }

    /// Exact target stratum object.
    #[must_use]
    pub const fn target(&self) -> DerivedStratumObjectV1 {
        self.target
    }

    /// Identity versus nonempty declared component path.
    #[must_use]
    pub const fn class(&self) -> AdmittedDerivedStratumMorphismClassV1 {
        self.class
    }

    /// Ordered component declarations with exact scoped endpoints.
    #[must_use]
    pub fn primitive_path(&self) -> &[DeclaredStratumMapPrimitiveV1] {
        &self.primitive_path
    }

    /// Ordered primitive receipt identities after associative flattening.
    #[must_use]
    pub fn primitive_factors(&self) -> &[DerivedStratumMorphismIdV1] {
        &self.primitive_factors
    }

    /// Ordered artifacts denying whole-map and theorem authority.
    #[must_use]
    pub fn no_authority_claims(&self) -> &[DerivedNoClaimIdV1] {
        &self.no_authority_claims
    }

    /// Typed stratum-scoped morphism identity.
    #[must_use]
    pub const fn id(&self) -> DerivedStratumMorphismIdV1 {
        self.receipt.id()
    }

    /// Canonical receipt and construction limits.
    #[must_use]
    pub const fn identity_receipt(&self) -> IdentityReceipt<DerivedStratumMorphismIdV1> {
        self.receipt
    }
}

/// Sealed exhaustive finite source-stratum coverage candidate.
///
/// This token binds one direct sealed component to every source stratum in
/// canonical order. It does not admit a whole-geometry map, authenticate the
/// nominal assembly or constructibility declarations, or expose composition,
/// inversion, evidence transport, continuity, or incidence authority.
#[derive(Debug, PartialEq, Eq)]
pub struct AdmittedDerivedExhaustiveStratifiedMapCandidateV1 {
    source_geometry: DerivedGeometryIdV1,
    source_stratification: StratificationIdV1,
    target_geometry: DerivedGeometryIdV1,
    target_stratification: StratificationIdV1,
    components: Vec<DerivedStratifiedMapComponentBindingV1>,
    nominal_assembly: DerivedStratifiedMapAssemblyIdV1,
    nominal_constructibility: DerivedGlobalConstructibilityDeclarationIdV1,
    no_authority: DerivedNoClaimIdV1,
    receipt: IdentityReceipt<DerivedExhaustiveStratifiedMapCandidateIdV1>,
}

impl AdmittedDerivedExhaustiveStratifiedMapCandidateV1 {
    /// Exact source geometry.
    #[must_use]
    pub const fn source_geometry(&self) -> DerivedGeometryIdV1 {
        self.source_geometry
    }

    /// Exact source stratification.
    #[must_use]
    pub const fn source_stratification(&self) -> StratificationIdV1 {
        self.source_stratification
    }

    /// Exact target geometry.
    #[must_use]
    pub const fn target_geometry(&self) -> DerivedGeometryIdV1 {
        self.target_geometry
    }

    /// Exact target stratification.
    #[must_use]
    pub const fn target_stratification(&self) -> StratificationIdV1 {
        self.target_stratification
    }

    /// Canonically source-ordered direct component bindings.
    #[must_use]
    pub fn components(&self) -> &[DerivedStratifiedMapComponentBindingV1] {
        &self.components
    }

    /// Nominal whole-family artifact; not authenticated by this token.
    #[must_use]
    pub const fn nominal_assembly(&self) -> DerivedStratifiedMapAssemblyIdV1 {
        self.nominal_assembly
    }

    /// Nominal global constructibility declaration; not authenticated here.
    #[must_use]
    pub const fn nominal_constructibility(&self) -> DerivedGlobalConstructibilityDeclarationIdV1 {
        self.nominal_constructibility
    }

    /// Explicit artifact denying global map and theorem authority.
    #[must_use]
    pub const fn no_authority(&self) -> DerivedNoClaimIdV1 {
        self.no_authority
    }

    /// Typed exhaustive assembly-candidate identity.
    #[must_use]
    pub const fn id(&self) -> DerivedExhaustiveStratifiedMapCandidateIdV1 {
        self.receipt.id()
    }

    /// Canonical receipt and construction limits.
    #[must_use]
    pub const fn identity_receipt(
        &self,
    ) -> IdentityReceipt<DerivedExhaustiveStratifiedMapCandidateIdV1> {
        self.receipt
    }
}

/// Sealed, content-addressed standalone span `source <- apex -> target`.
///
/// This token proves only exact leg identity and orientation binding. It does
/// not expose a direct source-to-target evidence transport or certify that the
/// span is functional, invertible, a pullback, or physically meaningful.
#[derive(Debug, PartialEq, Eq)]
pub struct AdmittedDerivedSpanCorrespondenceV1 {
    source: DerivedGeometryIdV1,
    apex: DerivedGeometryIdV1,
    target: DerivedGeometryIdV1,
    left_leg: DerivedMorphismIdV1,
    right_leg: DerivedMorphismIdV1,
    no_claim: DerivedNoClaimIdV1,
    receipt: IdentityReceipt<DerivedSpanCorrespondenceIdV1>,
}

impl AdmittedDerivedSpanCorrespondenceV1 {
    /// Exact left endpoint.
    #[must_use]
    pub const fn source(&self) -> DerivedGeometryIdV1 {
        self.source
    }

    /// Exact common apex.
    #[must_use]
    pub const fn apex(&self) -> DerivedGeometryIdV1 {
        self.apex
    }

    /// Exact right endpoint.
    #[must_use]
    pub const fn target(&self) -> DerivedGeometryIdV1 {
        self.target
    }

    /// Exact admitted `apex -> source` morphism identity.
    #[must_use]
    pub const fn left_leg(&self) -> DerivedMorphismIdV1 {
        self.left_leg
    }

    /// Exact admitted `apex -> target` morphism identity.
    #[must_use]
    pub const fn right_leg(&self) -> DerivedMorphismIdV1 {
        self.right_leg
    }

    /// Explicit no-correspondence-authority artifact.
    #[must_use]
    pub const fn no_claim(&self) -> DerivedNoClaimIdV1 {
        self.no_claim
    }

    /// Typed standalone span identity.
    #[must_use]
    pub const fn id(&self) -> DerivedSpanCorrespondenceIdV1 {
        self.receipt.id()
    }

    /// Canonical receipt and construction limits.
    #[must_use]
    pub const fn identity_receipt(&self) -> IdentityReceipt<DerivedSpanCorrespondenceIdV1> {
        self.receipt
    }
}

/// Sealed structural declaration of one fixed-resolution quasi-isomorphism candidate.
///
/// This token binds exact local presentations and nominal external-check metadata
/// to an already sealed refinement path. It deliberately exposes no equivalence,
/// inverse, homotopy, composition, evidence-transport, or authority capability.
#[derive(Debug, PartialEq, Eq)]
pub struct AdmittedDerivedFixedResolutionQuasiIsomorphismCandidateV1 {
    refinement_path: DerivedMorphismIdV1,
    source_geometry: DerivedGeometryIdV1,
    target_geometry: DerivedGeometryIdV1,
    source_local_model: DerivedLocalModelIdV1,
    target_local_model: DerivedLocalModelIdV1,
    complex_role: DerivedComplexRoleV1,
    source_complex: DerivedComplexIdV1,
    target_complex: DerivedComplexIdV1,
    source_resolution: DerivedResolutionIdV1,
    target_resolution: DerivedResolutionIdV1,
    source_scope_witness: DerivedWitnessIdV1,
    target_scope_witness: DerivedWitnessIdV1,
    nominal_theorem: DerivedTheoremIdV1,
    nominal_checker: DerivedCheckerIdV1,
    nominal_check_receipt: DerivedWitnessIdV1,
    no_authority: DerivedNoClaimIdV1,
    receipt: IdentityReceipt<DerivedFixedResolutionQuasiIsomorphismCandidateIdV1>,
}

impl AdmittedDerivedFixedResolutionQuasiIsomorphismCandidateV1 {
    /// Exact sealed refinement path used as the directed candidate.
    #[must_use]
    pub const fn refinement_path(&self) -> DerivedMorphismIdV1 {
        self.refinement_path
    }

    /// Exact source geometry.
    #[must_use]
    pub const fn source_geometry(&self) -> DerivedGeometryIdV1 {
        self.source_geometry
    }

    /// Exact target geometry.
    #[must_use]
    pub const fn target_geometry(&self) -> DerivedGeometryIdV1 {
        self.target_geometry
    }

    /// Exact source local presentation.
    #[must_use]
    pub const fn source_local_model(&self) -> DerivedLocalModelIdV1 {
        self.source_local_model
    }

    /// Exact target local presentation.
    #[must_use]
    pub const fn target_local_model(&self) -> DerivedLocalModelIdV1 {
        self.target_local_model
    }

    /// Complex role under comparison.
    #[must_use]
    pub const fn complex_role(&self) -> DerivedComplexRoleV1 {
        self.complex_role
    }

    /// Exact source complex.
    #[must_use]
    pub const fn source_complex(&self) -> DerivedComplexIdV1 {
        self.source_complex
    }

    /// Exact target complex.
    #[must_use]
    pub const fn target_complex(&self) -> DerivedComplexIdV1 {
        self.target_complex
    }

    /// Exact source finite resolution.
    #[must_use]
    pub const fn source_resolution(&self) -> DerivedResolutionIdV1 {
        self.source_resolution
    }

    /// Exact target finite resolution.
    #[must_use]
    pub const fn target_resolution(&self) -> DerivedResolutionIdV1 {
        self.target_resolution
    }

    /// Fixed-resolution scope witness retained by the source local model.
    #[must_use]
    pub const fn source_scope_witness(&self) -> DerivedWitnessIdV1 {
        self.source_scope_witness
    }

    /// Fixed-resolution scope witness retained by the target local model.
    #[must_use]
    pub const fn target_scope_witness(&self) -> DerivedWitnessIdV1 {
        self.target_scope_witness
    }

    /// Nominal theorem card; not authenticated by this token.
    #[must_use]
    pub const fn nominal_theorem(&self) -> DerivedTheoremIdV1 {
        self.nominal_theorem
    }

    /// Nominal checker; not executed by this token.
    #[must_use]
    pub const fn nominal_checker(&self) -> DerivedCheckerIdV1 {
        self.nominal_checker
    }

    /// Nominal external check receipt; not authenticated by this token.
    #[must_use]
    pub const fn nominal_check_receipt(&self) -> DerivedWitnessIdV1 {
        self.nominal_check_receipt
    }

    /// Explicit artifact denying authority to this structural candidate.
    #[must_use]
    pub const fn no_authority(&self) -> DerivedNoClaimIdV1 {
        self.no_authority
    }

    /// Typed candidate identity.
    #[must_use]
    pub const fn id(&self) -> DerivedFixedResolutionQuasiIsomorphismCandidateIdV1 {
        self.receipt.id()
    }

    /// Canonical receipt and construction limits.
    #[must_use]
    pub const fn identity_receipt(
        &self,
    ) -> IdentityReceipt<DerivedFixedResolutionQuasiIsomorphismCandidateIdV1> {
        self.receipt
    }
}

#[derive(Debug, Clone, Copy)]
struct GeometryEndpointV1<'a> {
    id: DerivedGeometryIdV1,
    subject: DerivedSubjectIdV1,
    model_version: DerivedModelVersionIdV1,
    category: GeometricCategoryV1,
    coefficients: CoefficientSystemV1,
    frame: DerivedFrameIdV1,
    unit_system: DerivedUnitSystemIdV1,
    charts: &'a [ConfigurationChartV1],
    complexes: &'a [FiniteDerivedComplexV1],
}

impl<'a> GeometryEndpointV1<'a> {
    fn from_admitted(value: &'a AdmittedDerivedGeometryV1) -> Self {
        let ir = value.ir();
        Self {
            id: value.id(),
            subject: ir.subject,
            model_version: ir.model_version,
            category: ir.category,
            coefficients: ir.coefficients,
            frame: ir.frame,
            unit_system: ir.unit_system,
            charts: &ir.charts,
            complexes: &ir.complexes,
        }
    }
}

#[derive(Debug, Clone, Copy)]
enum ReceiptClassV1 {
    Identity,
    PrimitiveStrict(DerivedWitnessIdV1),
    CompositeStrict,
    PrimitiveDeclaredChartMap {
        source_chart: ConfigurationChartIdV1,
        target_chart: ConfigurationChartIdV1,
        overlap: DerivedChartOverlapIdV1,
        map: DerivedChartMapIdV1,
    },
    CompositeDeclaredChartMap,
    CompositeHeterogeneous,
    PrimitiveDeclaredComplexRefinement {
        source_complex: DerivedComplexIdV1,
        target_complex: DerivedComplexIdV1,
        source_resolution: DerivedResolutionIdV1,
        target_resolution: DerivedResolutionIdV1,
        prolongation: DerivedComplexRefinementMapIdV1,
        commutation: DerivedWitnessIdV1,
    },
    CompositeDeclaredComplexRefinement,
    PrimitiveDeclaredInclusion {
        map: DerivedInclusionMapIdV1,
        containment: DerivedWitnessIdV1,
    },
    CompositeDeclaredInclusion,
}

#[derive(Debug, Clone, Copy)]
enum StratumReceiptClassV1 {
    Identity,
    PrimitiveDeclaredComponent {
        map: DerivedStratumMapIdV1,
        constructibility: DerivedConstructibilityDeclarationIdV1,
    },
    CompositeDeclaredPath,
}

fn is_zero(bytes: &[u8; 32]) -> bool {
    bytes.iter().all(|byte| *byte == 0)
}

fn rank_tag(rank: ColorRank) -> u8 {
    match rank {
        ColorRank::Estimated => 0,
        ColorRank::Validated => 1,
        ColorRank::Verified => 2,
    }
}

enum EvidenceBytesV1 {
    Identity([u8; 1]),
    Strict([u8; 131]),
}

impl EvidenceBytesV1 {
    fn as_slice(&self) -> &[u8] {
        match self {
            Self::Identity(bytes) => bytes,
            Self::Strict(bytes) => bytes,
        }
    }
}

fn strict_evidence_bytes(
    tag: u8,
    input_geometry: DerivedGeometryIdV1,
    output_geometry: DerivedGeometryIdV1,
    input_evidence: DerivedEvidenceArtifactIdV1,
    output_evidence: DerivedEvidenceArtifactIdV1,
    input_rank: ColorRank,
    output_rank: ColorRank,
) -> EvidenceBytesV1 {
    let mut bytes = [0_u8; 131];
    bytes[0] = tag;
    bytes[1..33].copy_from_slice(input_geometry.as_bytes());
    bytes[33..65].copy_from_slice(output_geometry.as_bytes());
    bytes[65..97].copy_from_slice(input_evidence.as_bytes());
    bytes[97..129].copy_from_slice(output_evidence.as_bytes());
    bytes[129] = rank_tag(input_rank);
    bytes[130] = rank_tag(output_rank);
    EvidenceBytesV1::Strict(bytes)
}

fn evidence_bytes(evidence: DerivedEvidenceTransportV1) -> EvidenceBytesV1 {
    match evidence {
        DerivedEvidenceTransportV1::Identity => EvidenceBytesV1::Identity([0]),
        DerivedEvidenceTransportV1::RestrictionContravariant {
            input_geometry,
            output_geometry,
            input_evidence,
            output_evidence,
            input_rank,
            output_rank,
        } => strict_evidence_bytes(
            1,
            input_geometry,
            output_geometry,
            input_evidence,
            output_evidence,
            input_rank,
            output_rank,
        ),
        DerivedEvidenceTransportV1::BalanceCorestrictionCovariant {
            input_geometry,
            output_geometry,
            input_evidence,
            output_evidence,
            input_rank,
            output_rank,
        } => strict_evidence_bytes(
            2,
            input_geometry,
            output_geometry,
            input_evidence,
            output_evidence,
            input_rank,
            output_rank,
        ),
    }
}

enum ClassBytesV1 {
    Tag([u8; 1]),
    Primitive([u8; 33]),
    DeclaredInclusion([u8; 65]),
    DeclaredChartMap([u8; 129]),
    DeclaredComplexRefinement([u8; 193]),
}

impl ClassBytesV1 {
    fn as_slice(&self) -> &[u8] {
        match self {
            Self::Tag(bytes) => bytes,
            Self::Primitive(bytes) => bytes,
            Self::DeclaredInclusion(bytes) => bytes,
            Self::DeclaredChartMap(bytes) => bytes,
            Self::DeclaredComplexRefinement(bytes) => bytes,
        }
    }
}

fn class_bytes(class: ReceiptClassV1) -> ClassBytesV1 {
    match class {
        ReceiptClassV1::Identity => ClassBytesV1::Tag([0]),
        ReceiptClassV1::PrimitiveStrict(witness) => {
            let mut bytes = [0_u8; 33];
            bytes[0] = 1;
            bytes[1..].copy_from_slice(witness.as_bytes());
            ClassBytesV1::Primitive(bytes)
        }
        ReceiptClassV1::CompositeStrict => ClassBytesV1::Tag([2]),
        ReceiptClassV1::PrimitiveDeclaredChartMap {
            source_chart,
            target_chart,
            overlap,
            map,
        } => {
            let mut bytes = [0_u8; 129];
            bytes[0] = 3;
            bytes[1..33].copy_from_slice(source_chart.as_bytes());
            bytes[33..65].copy_from_slice(target_chart.as_bytes());
            bytes[65..97].copy_from_slice(overlap.as_bytes());
            bytes[97..129].copy_from_slice(map.as_bytes());
            ClassBytesV1::DeclaredChartMap(bytes)
        }
        ReceiptClassV1::CompositeDeclaredChartMap => ClassBytesV1::Tag([4]),
        ReceiptClassV1::CompositeHeterogeneous => ClassBytesV1::Tag([5]),
        ReceiptClassV1::PrimitiveDeclaredComplexRefinement {
            source_complex,
            target_complex,
            source_resolution,
            target_resolution,
            prolongation,
            commutation,
        } => {
            let mut bytes = [0_u8; 193];
            bytes[0] = 6;
            bytes[1..33].copy_from_slice(source_complex.as_bytes());
            bytes[33..65].copy_from_slice(target_complex.as_bytes());
            bytes[65..97].copy_from_slice(source_resolution.as_bytes());
            bytes[97..129].copy_from_slice(target_resolution.as_bytes());
            bytes[129..161].copy_from_slice(prolongation.as_bytes());
            bytes[161..193].copy_from_slice(commutation.as_bytes());
            ClassBytesV1::DeclaredComplexRefinement(bytes)
        }
        ReceiptClassV1::CompositeDeclaredComplexRefinement => ClassBytesV1::Tag([7]),
        ReceiptClassV1::PrimitiveDeclaredInclusion { map, containment } => {
            let mut bytes = [0_u8; 65];
            bytes[0] = 8;
            bytes[1..33].copy_from_slice(map.as_bytes());
            bytes[33..65].copy_from_slice(containment.as_bytes());
            ClassBytesV1::DeclaredInclusion(bytes)
        }
        ReceiptClassV1::CompositeDeclaredInclusion => ClassBytesV1::Tag([9]),
    }
}

enum StratumClassBytesV1 {
    Tag([u8; 1]),
    Primitive([u8; 65]),
}

impl StratumClassBytesV1 {
    fn as_slice(&self) -> &[u8] {
        match self {
            Self::Tag(bytes) => bytes,
            Self::Primitive(bytes) => bytes,
        }
    }
}

fn stratum_class_bytes(class: StratumReceiptClassV1) -> StratumClassBytesV1 {
    match class {
        StratumReceiptClassV1::Identity => StratumClassBytesV1::Tag([0]),
        StratumReceiptClassV1::PrimitiveDeclaredComponent {
            map,
            constructibility,
        } => {
            let mut bytes = [0_u8; 65];
            bytes[0] = 1;
            bytes[1..33].copy_from_slice(map.as_bytes());
            bytes[33..65].copy_from_slice(constructibility.as_bytes());
            StratumClassBytesV1::Primitive(bytes)
        }
        StratumReceiptClassV1::CompositeDeclaredPath => StratumClassBytesV1::Tag([2]),
    }
}

fn morphism_receipt(
    source: DerivedGeometryIdV1,
    target: DerivedGeometryIdV1,
    class: ReceiptClassV1,
    evidence: DerivedEvidenceTransportV1,
    no_claims: &[DerivedNoClaimIdV1],
    factors: &[DerivedMorphismIdV1],
    cx: &Cx<'_>,
) -> Result<IdentityReceipt<DerivedMorphismIdV1>, DerivedMorphismErrorV1> {
    if cx.checkpoint().is_err() {
        return Err(DerivedMorphismErrorV1::Cancelled {
            stage: "identity-entry",
        });
    }
    let class = class_bytes(class);
    let evidence = evidence_bytes(evidence);
    let map_identity_error = |error| match error {
        CanonicalError::Cancelled { .. } => DerivedMorphismErrorV1::Cancelled { stage: "identity" },
        other => DerivedMorphismErrorV1::Identity(other),
    };
    CanonicalEncoder::<DerivedMorphismIdV1, _>::new(DERIVED_MORPHISM_IDENTITY_LIMITS_V1, || {
        cx.checkpoint().is_err()
    })
    .map_err(map_identity_error)?
    .bytes(Field::new(0, "source"), source.as_bytes())
    .and_then(|encoder| encoder.bytes(Field::new(1, "target"), target.as_bytes()))
    .and_then(|encoder| encoder.bytes(Field::new(2, "class"), class.as_slice()))
    .and_then(|encoder| encoder.bytes(Field::new(3, "evidence-transport"), evidence.as_slice()))
    .and_then(|encoder| {
        encoder.ordered_bytes(
            Field::new(4, "no-equivalence-claims"),
            no_claims.len() as u64,
            no_claims.iter().map(|claim| &claim.as_bytes()[..]),
        )
    })
    .and_then(|encoder| {
        encoder.ordered_bytes(
            Field::new(5, "primitive-lineage"),
            factors.len() as u64,
            factors.iter().map(|factor| &factor.as_bytes()[..]),
        )
    })
    .and_then(|encoder| encoder.finish())
    .map_err(map_identity_error)
}

fn stratum_morphism_receipt(
    source: DerivedStratumObjectV1,
    target: DerivedStratumObjectV1,
    class: StratumReceiptClassV1,
    no_claims: &[DerivedNoClaimIdV1],
    factors: &[DerivedStratumMorphismIdV1],
    cx: &Cx<'_>,
) -> Result<IdentityReceipt<DerivedStratumMorphismIdV1>, DerivedStratumMorphismErrorV1> {
    if cx.checkpoint().is_err() {
        return Err(DerivedStratumMorphismErrorV1::Cancelled {
            stage: "stratum-identity-entry",
        });
    }
    let class = stratum_class_bytes(class);
    let map_identity_error = |error| match error {
        CanonicalError::Cancelled { .. } => DerivedStratumMorphismErrorV1::Cancelled {
            stage: "stratum-identity",
        },
        other => DerivedStratumMorphismErrorV1::Identity(other),
    };
    CanonicalEncoder::<DerivedStratumMorphismIdV1, _>::new(
        DERIVED_STRATUM_MORPHISM_IDENTITY_LIMITS_V1,
        || cx.checkpoint().is_err(),
    )
    .map_err(map_identity_error)?
    .bytes(Field::new(0, "source-geometry"), source.geometry.as_bytes())
    .and_then(|encoder| {
        encoder.bytes(
            Field::new(1, "source-stratification"),
            source.stratification.as_bytes(),
        )
    })
    .and_then(|encoder| encoder.bytes(Field::new(2, "source-stratum"), source.stratum.as_bytes()))
    .and_then(|encoder| encoder.bytes(Field::new(3, "target-geometry"), target.geometry.as_bytes()))
    .and_then(|encoder| {
        encoder.bytes(
            Field::new(4, "target-stratification"),
            target.stratification.as_bytes(),
        )
    })
    .and_then(|encoder| encoder.bytes(Field::new(5, "target-stratum"), target.stratum.as_bytes()))
    .and_then(|encoder| encoder.bytes(Field::new(6, "class"), class.as_slice()))
    .and_then(|encoder| {
        encoder.ordered_bytes(
            Field::new(7, "no-authority-claims"),
            no_claims.len() as u64,
            no_claims.iter().map(|claim| &claim.as_bytes()[..]),
        )
    })
    .and_then(|encoder| {
        encoder.ordered_bytes(
            Field::new(8, "primitive-lineage"),
            factors.len() as u64,
            factors.iter().map(|factor| &factor.as_bytes()[..]),
        )
    })
    .and_then(|encoder| encoder.finish())
    .map_err(map_identity_error)
}

fn exhaustive_stratified_map_candidate_receipt(
    ir: &DerivedExhaustiveStratifiedMapCandidateIrV1,
    components: &[DerivedStratifiedMapComponentBindingV1],
    cx: &Cx<'_>,
) -> Result<
    IdentityReceipt<DerivedExhaustiveStratifiedMapCandidateIdV1>,
    DerivedExhaustiveStratifiedMapCandidateErrorV1,
> {
    if cx.checkpoint().is_err() {
        return Err(DerivedExhaustiveStratifiedMapCandidateErrorV1::Cancelled {
            stage: "stratified-assembly-identity-entry",
        });
    }
    let map_identity_error = |error| match error {
        CanonicalError::Cancelled { .. } => {
            DerivedExhaustiveStratifiedMapCandidateErrorV1::Cancelled {
                stage: "stratified-assembly-identity",
            }
        }
        other => DerivedExhaustiveStratifiedMapCandidateErrorV1::Identity(other),
    };
    CanonicalEncoder::<DerivedExhaustiveStratifiedMapCandidateIdV1, _>::new(
        DERIVED_EXHAUSTIVE_STRATIFIED_MAP_CANDIDATE_IDENTITY_LIMITS_V1,
        || cx.checkpoint().is_err(),
    )
    .map_err(map_identity_error)?
    .bytes(
        Field::new(0, "source-geometry"),
        ir.source_geometry.as_bytes(),
    )
    .and_then(|encoder| {
        encoder.bytes(
            Field::new(1, "source-stratification"),
            ir.source_stratification.as_bytes(),
        )
    })
    .and_then(|encoder| {
        encoder.bytes(
            Field::new(2, "target-geometry"),
            ir.target_geometry.as_bytes(),
        )
    })
    .and_then(|encoder| {
        encoder.bytes(
            Field::new(3, "target-stratification"),
            ir.target_stratification.as_bytes(),
        )
    })
    .and_then(|encoder| {
        encoder.ordered_bytes(
            Field::new(4, "components"),
            components.len() as u64,
            components
                .iter()
                .map(|binding| &binding.component.as_bytes()[..]),
        )
    })
    .and_then(|encoder| {
        encoder.bytes(
            Field::new(5, "nominal-assembly"),
            ir.nominal_assembly.as_bytes(),
        )
    })
    .and_then(|encoder| {
        encoder.bytes(
            Field::new(6, "nominal-constructibility"),
            ir.nominal_constructibility.as_bytes(),
        )
    })
    .and_then(|encoder| encoder.bytes(Field::new(7, "no-authority"), ir.no_authority.as_bytes()))
    .and_then(|encoder| encoder.finish())
    .map_err(map_identity_error)
}

fn validate_evidence(
    source: DerivedGeometryIdV1,
    target: DerivedGeometryIdV1,
    evidence: DerivedEvidenceTransportV1,
) -> Result<(), DerivedMorphismErrorV1> {
    let (
        input_geometry,
        output_geometry,
        input_evidence,
        output_evidence,
        input_rank,
        output_rank,
        expected_input,
        expected_output,
    ) = match evidence {
        DerivedEvidenceTransportV1::Identity => {
            return if source == target {
                Ok(())
            } else {
                Err(DerivedMorphismErrorV1::EvidenceOrientationMismatch)
            };
        }
        DerivedEvidenceTransportV1::RestrictionContravariant {
            input_geometry,
            output_geometry,
            input_evidence,
            output_evidence,
            input_rank,
            output_rank,
        } => (
            input_geometry,
            output_geometry,
            input_evidence,
            output_evidence,
            input_rank,
            output_rank,
            target,
            source,
        ),
        DerivedEvidenceTransportV1::BalanceCorestrictionCovariant {
            input_geometry,
            output_geometry,
            input_evidence,
            output_evidence,
            input_rank,
            output_rank,
        } => (
            input_geometry,
            output_geometry,
            input_evidence,
            output_evidence,
            input_rank,
            output_rank,
            source,
            target,
        ),
    };
    if input_geometry != expected_input || output_geometry != expected_output {
        return Err(DerivedMorphismErrorV1::EvidenceOrientationMismatch);
    }
    if is_zero(input_evidence.as_bytes()) {
        return Err(DerivedMorphismErrorV1::MissingIdentity {
            field: "input-evidence",
        });
    }
    if is_zero(output_evidence.as_bytes()) {
        return Err(DerivedMorphismErrorV1::MissingIdentity {
            field: "output-evidence",
        });
    }
    if output_rank > input_rank {
        return Err(DerivedMorphismErrorV1::EvidenceStrengthening {
            input: input_rank,
            output: output_rank,
        });
    }
    Ok(())
}

fn shared_nonidentity_compatibility(
    source: GeometryEndpointV1<'_>,
    target: GeometryEndpointV1<'_>,
) -> Result<(), DerivedMorphismErrorV1> {
    if source.subject != target.subject {
        return Err(DerivedMorphismErrorV1::SubjectMismatch);
    }
    if source.model_version != target.model_version {
        return Err(DerivedMorphismErrorV1::ModelVersionMismatch);
    }
    if source.category != target.category {
        return Err(DerivedMorphismErrorV1::CategoryMismatch);
    }
    if source.coefficients != target.coefficients {
        return Err(DerivedMorphismErrorV1::CoefficientMismatch);
    }
    if source.frame != target.frame {
        return Err(DerivedMorphismErrorV1::FrameMismatch);
    }
    if source.unit_system != target.unit_system {
        return Err(DerivedMorphismErrorV1::UnitSystemMismatch);
    }
    Ok(())
}

fn declared_chart_path(
    source: GeometryEndpointV1<'_>,
    target: GeometryEndpointV1<'_>,
    source_chart: ConfigurationChartIdV1,
    target_chart: ConfigurationChartIdV1,
) -> Result<DerivedChartPathV1, DerivedMorphismErrorV1> {
    if is_zero(source_chart.as_bytes()) {
        return Err(DerivedMorphismErrorV1::MissingIdentity {
            field: "source-chart",
        });
    }
    if is_zero(target_chart.as_bytes()) {
        return Err(DerivedMorphismErrorV1::MissingIdentity {
            field: "target-chart",
        });
    }
    let source_spec = source
        .charts
        .iter()
        .find(|chart| chart.id == source_chart)
        .ok_or(DerivedMorphismErrorV1::MissingChart {
            field: "source-chart",
        })?;
    let target_spec = target
        .charts
        .iter()
        .find(|chart| chart.id == target_chart)
        .ok_or(DerivedMorphismErrorV1::MissingChart {
            field: "target-chart",
        })?;
    if source_spec.coordinate_dimension != target_spec.coordinate_dimension
        || source_spec.ambient_dimension != target_spec.ambient_dimension
    {
        return Err(DerivedMorphismErrorV1::ChartDimensionMismatch);
    }
    if source_spec.frame != source.frame
        || target_spec.frame != target.frame
        || source_spec.coordinates.system != source.unit_system
        || target_spec.coordinates.system != target.unit_system
        || source_spec.coordinates.quantity != target_spec.coordinates.quantity
        || source_spec.coordinates.scale_to_canonical.to_bits()
            != target_spec.coordinates.scale_to_canonical.to_bits()
    {
        return Err(DerivedMorphismErrorV1::ChartConventionMismatch);
    }
    Ok(DerivedChartPathV1 {
        source_chart,
        target_chart,
    })
}

const fn truncation_refinement_progress(source: u32, target: u32) -> Option<bool> {
    match (source, target) {
        (0, 0) => Some(false),
        (0, _) => None,
        (_, 0) => Some(true),
        (source, target) if target >= source => Some(target > source),
        _ => None,
    }
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)] // One bounded rank-envelope admission.
fn declared_complex_refinement(
    source: GeometryEndpointV1<'_>,
    target: GeometryEndpointV1<'_>,
    source_complex: DerivedComplexIdV1,
    target_complex: DerivedComplexIdV1,
    source_resolution: DerivedResolutionIdV1,
    target_resolution: DerivedResolutionIdV1,
    prolongation: DerivedComplexRefinementMapIdV1,
    commutation: DerivedWitnessIdV1,
    cx: &Cx<'_>,
) -> Result<DeclaredComplexRefinementPrimitiveV1, DerivedMorphismErrorV1> {
    let source_spec = source
        .complexes
        .iter()
        .find(|complex| complex.id == source_complex)
        .ok_or(DerivedMorphismErrorV1::MissingComplex {
            field: "source-complex",
        })?;
    let target_spec = target
        .complexes
        .iter()
        .find(|complex| complex.id == target_complex)
        .ok_or(DerivedMorphismErrorV1::MissingComplex {
            field: "target-complex",
        })?;
    for (matches, field) in [
        (
            source_spec.resolution.id == source_resolution,
            "source-resolution",
        ),
        (
            target_spec.resolution.id == target_resolution,
            "target-resolution",
        ),
        (source_spec.role == target_spec.role, "complex-role"),
        (source_spec.chart == target_spec.chart, "complex-chart"),
    ] {
        if !matches {
            return Err(DerivedMorphismErrorV1::ComplexRefinementMismatch { field });
        }
    }
    let source_chart = source
        .charts
        .iter()
        .find(|chart| chart.id == source_spec.chart)
        .ok_or(DerivedMorphismErrorV1::MissingChart {
            field: "source-complex-chart",
        })?;
    let target_chart = target
        .charts
        .iter()
        .find(|chart| chart.id == target_spec.chart)
        .ok_or(DerivedMorphismErrorV1::MissingChart {
            field: "target-complex-chart",
        })?;
    if source_chart != target_chart {
        return Err(DerivedMorphismErrorV1::ComplexRefinementMismatch {
            field: "chart-semantics",
        });
    }

    let mut target_index = 0usize;
    let mut strict_progress = truncation_refinement_progress(
        source_spec.resolution.truncation_order,
        target_spec.resolution.truncation_order,
    )
    .ok_or(DerivedMorphismErrorV1::ComplexRefinementMismatch {
        field: "truncation-policy",
    })?;
    for (completed, source_space) in source_spec.spaces.iter().enumerate() {
        if completed.is_multiple_of(DERIVED_MORPHISM_CANCELLATION_STRIDE_V1)
            && cx.checkpoint().is_err()
        {
            return Err(DerivedMorphismErrorV1::Cancelled {
                stage: "complex-refinement-shape",
            });
        }
        while let Some(extra_space) = target_spec
            .spaces
            .get(target_index)
            .filter(|space| space.degree < source_space.degree)
        {
            if target_index.is_multiple_of(DERIVED_MORPHISM_CANCELLATION_STRIDE_V1)
                && cx.checkpoint().is_err()
            {
                return Err(DerivedMorphismErrorV1::Cancelled {
                    stage: "complex-refinement-shape",
                });
            }
            strict_progress |= extra_space.dimension > 0;
            target_index += 1;
        }
        let target_space = target_spec.spaces.get(target_index).ok_or(
            DerivedMorphismErrorV1::ComplexRefinementMismatch {
                field: "degree-coverage",
            },
        )?;
        if target_space.degree != source_space.degree {
            return Err(DerivedMorphismErrorV1::ComplexRefinementMismatch {
                field: "degree-coverage",
            });
        }
        if target_space.quantity != source_space.quantity {
            return Err(DerivedMorphismErrorV1::ComplexRefinementMismatch {
                field: "degree-quantity",
            });
        }
        if target_space.dimension < source_space.dimension {
            return Err(DerivedMorphismErrorV1::ComplexRefinementMismatch {
                field: "degree-rank",
            });
        }
        strict_progress |= target_space.dimension > source_space.dimension;
        target_index += 1;
    }
    for (completed, extra_space) in target_spec.spaces[target_index..].iter().enumerate() {
        if completed.is_multiple_of(DERIVED_MORPHISM_CANCELLATION_STRIDE_V1)
            && cx.checkpoint().is_err()
        {
            return Err(DerivedMorphismErrorV1::Cancelled {
                stage: "complex-refinement-shape",
            });
        }
        strict_progress |= extra_space.dimension > 0;
    }
    if !strict_progress {
        return Err(DerivedMorphismErrorV1::ComplexRefinementMismatch {
            field: "strict-progress",
        });
    }
    if cx.checkpoint().is_err() {
        return Err(DerivedMorphismErrorV1::Cancelled {
            stage: "complex-refinement-shape",
        });
    }
    Ok(DeclaredComplexRefinementPrimitiveV1 {
        source_geometry: source.id,
        target_geometry: target.id,
        source_complex,
        target_complex,
        source_resolution,
        target_resolution,
        prolongation,
        commutation,
    })
}

#[derive(Debug, Clone, Copy)]
struct ValidatedMorphismClassV1 {
    admitted: AdmittedDerivedMorphismClassV1,
    receipt: ReceiptClassV1,
    no_claim: Option<DerivedNoClaimIdV1>,
    chart_path: Option<DerivedChartPathV1>,
    primitive: Option<AdmittedDerivedPrimitiveV1>,
    chart_primitive: Option<DeclaredChartMapPrimitiveV1>,
}

#[allow(clippy::too_many_lines)] // One exhaustive primitive-family admission dispatch.
fn validate_morphism_class(
    ir: DerivedMorphismIrV1,
    source: GeometryEndpointV1<'_>,
    target: GeometryEndpointV1<'_>,
    cx: &Cx<'_>,
) -> Result<ValidatedMorphismClassV1, DerivedMorphismErrorV1> {
    match (ir.kind, ir.equivalence) {
        (DerivedMorphismKindV1::Identity, DerivedEquivalenceBoundaryV1::IdentityOnly) => {
            if source.id != target.id || ir.evidence != DerivedEvidenceTransportV1::Identity {
                return Err(DerivedMorphismErrorV1::InvalidIdentity);
            }
            Ok(ValidatedMorphismClassV1 {
                admitted: AdmittedDerivedMorphismClassV1::Identity,
                receipt: ReceiptClassV1::Identity,
                no_claim: None,
                chart_path: None,
                primitive: None,
                chart_primitive: None,
            })
        }
        (
            DerivedMorphismKindV1::Strict { witness },
            DerivedEquivalenceBoundaryV1::NoClaim { artifact },
        ) => {
            if is_zero(witness.as_bytes()) {
                return Err(DerivedMorphismErrorV1::MissingIdentity {
                    field: "strict-witness",
                });
            }
            if is_zero(artifact.as_bytes()) {
                return Err(DerivedMorphismErrorV1::MissingIdentity {
                    field: "no-equivalence-claim",
                });
            }
            if ir.evidence == DerivedEvidenceTransportV1::Identity {
                return Err(DerivedMorphismErrorV1::EvidenceOrientationMismatch);
            }
            shared_nonidentity_compatibility(source, target)?;
            Ok(ValidatedMorphismClassV1 {
                admitted: AdmittedDerivedMorphismClassV1::Strict,
                receipt: ReceiptClassV1::PrimitiveStrict(witness),
                no_claim: Some(artifact),
                chart_path: None,
                primitive: Some(AdmittedDerivedPrimitiveV1::Strict {
                    source_geometry: source.id,
                    target_geometry: target.id,
                    witness,
                }),
                chart_primitive: None,
            })
        }
        (
            DerivedMorphismKindV1::DeclaredChartMap {
                source_chart,
                target_chart,
                overlap,
                map,
            },
            DerivedEquivalenceBoundaryV1::NoClaim { artifact },
        ) => {
            for (bytes, field) in [
                (overlap.as_bytes(), "chart-overlap"),
                (map.as_bytes(), "chart-map"),
                (artifact.as_bytes(), "no-equivalence-claim"),
            ] {
                if is_zero(bytes) {
                    return Err(DerivedMorphismErrorV1::MissingIdentity { field });
                }
            }
            if ir.evidence == DerivedEvidenceTransportV1::Identity {
                return Err(DerivedMorphismErrorV1::EvidenceOrientationMismatch);
            }
            shared_nonidentity_compatibility(source, target)?;
            let chart_path = declared_chart_path(source, target, source_chart, target_chart)?;
            let chart_primitive = DeclaredChartMapPrimitiveV1 {
                source_geometry: source.id,
                target_geometry: target.id,
                source_chart,
                target_chart,
                overlap,
                map,
            };
            Ok(ValidatedMorphismClassV1 {
                admitted: AdmittedDerivedMorphismClassV1::DeclaredChartMapPath,
                receipt: ReceiptClassV1::PrimitiveDeclaredChartMap {
                    source_chart,
                    target_chart,
                    overlap,
                    map,
                },
                no_claim: Some(artifact),
                chart_path: Some(chart_path),
                primitive: Some(AdmittedDerivedPrimitiveV1::DeclaredChartMap(
                    chart_primitive,
                )),
                chart_primitive: Some(chart_primitive),
            })
        }
        (
            DerivedMorphismKindV1::DeclaredComplexRefinement {
                source_complex,
                target_complex,
                source_resolution,
                target_resolution,
                prolongation,
                commutation,
            },
            DerivedEquivalenceBoundaryV1::NoClaim { artifact },
        ) => {
            for (bytes, field) in [
                (source_complex.as_bytes(), "source-complex"),
                (target_complex.as_bytes(), "target-complex"),
                (source_resolution.as_bytes(), "source-resolution"),
                (target_resolution.as_bytes(), "target-resolution"),
                (prolongation.as_bytes(), "complex-refinement-prolongation"),
                (commutation.as_bytes(), "complex-refinement-commutation"),
                (artifact.as_bytes(), "no-equivalence-claim"),
            ] {
                if is_zero(bytes) {
                    return Err(DerivedMorphismErrorV1::MissingIdentity { field });
                }
            }
            if ir.evidence == DerivedEvidenceTransportV1::Identity {
                return Err(DerivedMorphismErrorV1::EvidenceOrientationMismatch);
            }
            shared_nonidentity_compatibility(source, target)?;
            let primitive = declared_complex_refinement(
                source,
                target,
                source_complex,
                target_complex,
                source_resolution,
                target_resolution,
                prolongation,
                commutation,
                cx,
            )?;
            Ok(ValidatedMorphismClassV1 {
                admitted: AdmittedDerivedMorphismClassV1::DeclaredComplexRefinementPath,
                receipt: ReceiptClassV1::PrimitiveDeclaredComplexRefinement {
                    source_complex,
                    target_complex,
                    source_resolution,
                    target_resolution,
                    prolongation,
                    commutation,
                },
                no_claim: Some(artifact),
                chart_path: None,
                primitive: Some(AdmittedDerivedPrimitiveV1::DeclaredComplexRefinement(
                    primitive,
                )),
                chart_primitive: None,
            })
        }
        (
            DerivedMorphismKindV1::DeclaredInclusion { map, containment },
            DerivedEquivalenceBoundaryV1::NoClaim { artifact },
        ) => {
            for (bytes, field) in [
                (map.as_bytes(), "inclusion-map"),
                (containment.as_bytes(), "inclusion-containment"),
                (artifact.as_bytes(), "no-equivalence-claim"),
            ] {
                if is_zero(bytes) {
                    return Err(DerivedMorphismErrorV1::MissingIdentity { field });
                }
            }
            if ir.evidence == DerivedEvidenceTransportV1::Identity {
                return Err(DerivedMorphismErrorV1::EvidenceOrientationMismatch);
            }
            shared_nonidentity_compatibility(source, target)?;
            let primitive = DeclaredInclusionPrimitiveV1 {
                source_geometry: source.id,
                target_geometry: target.id,
                map,
                containment,
            };
            Ok(ValidatedMorphismClassV1 {
                admitted: AdmittedDerivedMorphismClassV1::DeclaredInclusionPath,
                receipt: ReceiptClassV1::PrimitiveDeclaredInclusion { map, containment },
                no_claim: Some(artifact),
                chart_path: None,
                primitive: Some(AdmittedDerivedPrimitiveV1::DeclaredInclusion(primitive)),
                chart_primitive: None,
            })
        }
        (DerivedMorphismKindV1::Identity, _) => Err(DerivedMorphismErrorV1::InvalidIdentity),
        (
            DerivedMorphismKindV1::Strict { .. }
            | DerivedMorphismKindV1::DeclaredChartMap { .. }
            | DerivedMorphismKindV1::DeclaredComplexRefinement { .. }
            | DerivedMorphismKindV1::DeclaredInclusion { .. },
            _,
        ) => Err(DerivedMorphismErrorV1::EquivalenceLaundering),
    }
}

fn retain_no_claim(
    no_claim: Option<DerivedNoClaimIdV1>,
) -> Result<Vec<DerivedNoClaimIdV1>, DerivedMorphismErrorV1> {
    let mut retained = Vec::new();
    if let Some(claim) = no_claim {
        retained
            .try_reserve_exact(1)
            .map_err(|_| DerivedMorphismErrorV1::AllocationRefused {
                field: "no-equivalence-claims",
            })?;
        retained.push(claim);
    }
    Ok(retained)
}

fn retain_chart_primitive(
    primitive: Option<DeclaredChartMapPrimitiveV1>,
) -> Result<Vec<DeclaredChartMapPrimitiveV1>, DerivedMorphismErrorV1> {
    let mut retained = Vec::new();
    if let Some(primitive) = primitive {
        retained
            .try_reserve_exact(1)
            .map_err(|_| DerivedMorphismErrorV1::AllocationRefused {
                field: "declared-chart-map-lineage",
            })?;
        retained.push(primitive);
    }
    Ok(retained)
}

fn retain_typed_primitive(
    primitive: Option<AdmittedDerivedPrimitiveV1>,
) -> Result<Vec<AdmittedDerivedPrimitiveV1>, DerivedMorphismErrorV1> {
    let mut retained = Vec::new();
    if let Some(primitive) = primitive {
        retained
            .try_reserve_exact(1)
            .map_err(|_| DerivedMorphismErrorV1::AllocationRefused {
                field: "typed-primitive-lineage",
            })?;
        retained.push(primitive);
    }
    Ok(retained)
}

fn admit_between_endpoints(
    ir: DerivedMorphismIrV1,
    source: GeometryEndpointV1<'_>,
    target: GeometryEndpointV1<'_>,
    cx: &Cx<'_>,
) -> Result<AdmittedDerivedMorphismV1, DerivedMorphismErrorV1> {
    if cx.checkpoint().is_err() {
        return Err(DerivedMorphismErrorV1::Cancelled {
            stage: "admission-entry",
        });
    }
    if ir.schema_version != DERIVED_MORPHISM_SCHEMA_VERSION_V1 {
        return Err(DerivedMorphismErrorV1::UnsupportedSchemaVersion {
            found: ir.schema_version,
            supported: DERIVED_MORPHISM_SCHEMA_VERSION_V1,
        });
    }
    if ir.source != source.id {
        return Err(DerivedMorphismErrorV1::EndpointMismatch { field: "source" });
    }
    if ir.target != target.id {
        return Err(DerivedMorphismErrorV1::EndpointMismatch { field: "target" });
    }
    let validated = validate_morphism_class(ir, source, target, cx)?;
    validate_evidence(source.id, target.id, ir.evidence)?;

    if cx.checkpoint().is_err() {
        return Err(DerivedMorphismErrorV1::Cancelled { stage: "admission" });
    }
    let no_equivalence_claims = retain_no_claim(validated.no_claim)?;
    let primitive_path = retain_typed_primitive(validated.primitive)?;
    let declared_chart_maps = retain_chart_primitive(validated.chart_primitive)?;
    let receipt = morphism_receipt(
        source.id,
        target.id,
        validated.receipt,
        ir.evidence,
        &no_equivalence_claims,
        &[],
        cx,
    )?;
    let mut primitive_factors = Vec::new();
    if validated.admitted != AdmittedDerivedMorphismClassV1::Identity {
        primitive_factors.try_reserve_exact(1).map_err(|_| {
            DerivedMorphismErrorV1::AllocationRefused {
                field: "primitive-lineage",
            }
        })?;
        primitive_factors.push(receipt.id());
    }
    if cx.checkpoint().is_err() {
        return Err(DerivedMorphismErrorV1::Cancelled {
            stage: "publication",
        });
    }
    Ok(AdmittedDerivedMorphismV1 {
        source: source.id,
        target: target.id,
        class: validated.admitted,
        chart_path: validated.chart_path,
        primitive_path,
        declared_chart_maps,
        evidence: ir.evidence,
        no_equivalence_claims,
        primitive_factors,
        receipt,
    })
}

/// Admit one primitive identity, strict map, chart map, complex refinement, or inclusion.
///
/// This validates structural endpoint compatibility, caller-declared evidence
/// rank monotonicity, and family-specific chart or finite graded-rank envelopes.
/// Nominal evidence/map/witness identities are retained but not authenticated.
/// No primitive can mint equivalence, inverse, quasi-isomorphism, chart-map
/// invertibility, overlap coverage, chain commutation, subset containment,
/// injectivity, numerical error reduction, physical correspondence, or theorem
/// authority.
///
/// # Errors
/// Returns a typed refusal for endpoint/model/category/coefficient/frame/unit,
/// chart/complex ownership and shape, evidence-direction/rank,
/// equivalence-boundary, cancellation, allocation, or canonical-identity
/// defects.
#[must_use = "a raw morphism request has no structural authority"]
pub fn admit_derived_morphism_v1(
    ir: DerivedMorphismIrV1,
    source: &AdmittedDerivedGeometryV1,
    target: &AdmittedDerivedGeometryV1,
    cx: &Cx<'_>,
) -> Result<AdmittedDerivedMorphismV1, DerivedMorphismErrorV1> {
    admit_between_endpoints(
        ir,
        GeometryEndpointV1::from_admitted(source),
        GeometryEndpointV1::from_admitted(target),
        cx,
    )
}

/// Mint the exact identity arrow on one admitted geometry.
///
/// # Errors
/// Returns a cancellation, allocation, or canonical-identity refusal.
#[must_use = "identity construction must complete before composition"]
pub fn identity_derived_morphism_v1(
    object: &AdmittedDerivedGeometryV1,
    cx: &Cx<'_>,
) -> Result<AdmittedDerivedMorphismV1, DerivedMorphismErrorV1> {
    let endpoint = GeometryEndpointV1::from_admitted(object);
    admit_between_endpoints(
        DerivedMorphismIrV1 {
            schema_version: DERIVED_MORPHISM_SCHEMA_VERSION_V1,
            source: endpoint.id,
            target: endpoint.id,
            kind: DerivedMorphismKindV1::Identity,
            evidence: DerivedEvidenceTransportV1::Identity,
            equivalence: DerivedEquivalenceBoundaryV1::IdentityOnly,
        },
        endpoint,
        endpoint,
        cx,
    )
}

fn checked_combined_len(
    field: &'static str,
    left: usize,
    right: usize,
) -> Result<usize, DerivedMorphismErrorV1> {
    let requested = left
        .checked_add(right)
        .ok_or(DerivedMorphismErrorV1::ResourceLimit {
            field,
            requested: usize::MAX,
            limit: DERIVED_MORPHISM_MAX_FACTORS_V1,
        })?;
    if requested > DERIVED_MORPHISM_MAX_FACTORS_V1 {
        return Err(DerivedMorphismErrorV1::ResourceLimit {
            field,
            requested,
            limit: DERIVED_MORPHISM_MAX_FACTORS_V1,
        });
    }
    Ok(requested)
}

fn combine_slices<T: Copy>(
    field: &'static str,
    left: &[T],
    right: &[T],
    cx: &Cx<'_>,
) -> Result<Vec<T>, DerivedMorphismErrorV1> {
    let len = checked_combined_len(field, left.len(), right.len())?;
    if cx.checkpoint().is_err() {
        return Err(DerivedMorphismErrorV1::Cancelled { stage: field });
    }
    let mut out = Vec::new();
    out.try_reserve_exact(len)
        .map_err(|_| DerivedMorphismErrorV1::AllocationRefused { field })?;
    for (index, value) in left.iter().chain(right).copied().enumerate() {
        out.push(value);
        if (index + 1) % DERIVED_MORPHISM_CANCELLATION_STRIDE_V1 == 0 && cx.checkpoint().is_err() {
            return Err(DerivedMorphismErrorV1::Cancelled { stage: field });
        }
    }
    if cx.checkpoint().is_err() {
        return Err(DerivedMorphismErrorV1::Cancelled { stage: field });
    }
    Ok(out)
}

fn copy_admitted_morphism(
    value: &AdmittedDerivedMorphismV1,
    cx: &Cx<'_>,
) -> Result<AdmittedDerivedMorphismV1, DerivedMorphismErrorV1> {
    let primitive_path = combine_slices(
        "typed-primitive-lineage-copy",
        &value.primitive_path,
        &[],
        cx,
    )?;
    let primitive_factors =
        combine_slices("primitive-lineage-copy", &value.primitive_factors, &[], cx)?;
    let no_equivalence_claims = combine_slices(
        "no-equivalence-claims-copy",
        &value.no_equivalence_claims,
        &[],
        cx,
    )?;
    let declared_chart_maps = combine_slices(
        "declared-chart-map-lineage-copy",
        &value.declared_chart_maps,
        &[],
        cx,
    )?;
    if cx.checkpoint().is_err() {
        return Err(DerivedMorphismErrorV1::Cancelled {
            stage: "composition-publication",
        });
    }
    Ok(AdmittedDerivedMorphismV1 {
        source: value.source,
        target: value.target,
        class: value.class,
        chart_path: value.chart_path,
        primitive_path,
        declared_chart_maps,
        evidence: value.evidence,
        no_equivalence_claims,
        primitive_factors,
        receipt: value.receipt,
    })
}

fn validate_typed_primitive_seam(
    first: &AdmittedDerivedMorphismV1,
    second: &AdmittedDerivedMorphismV1,
) -> Result<(), DerivedMorphismErrorV1> {
    let first_start = first
        .primitive_path
        .first()
        .ok_or(DerivedMorphismErrorV1::CompositionClassMismatch)?;
    let first_end = first
        .primitive_path
        .last()
        .ok_or(DerivedMorphismErrorV1::CompositionClassMismatch)?;
    let second_start = second
        .primitive_path
        .first()
        .ok_or(DerivedMorphismErrorV1::CompositionClassMismatch)?;
    let second_end = second
        .primitive_path
        .last()
        .ok_or(DerivedMorphismErrorV1::CompositionClassMismatch)?;
    if first_start.source_geometry() != first.source
        || first_end.target_geometry() != first.target
        || second_start.source_geometry() != second.source
        || second_end.target_geometry() != second.target
        || first_end.target_geometry() != second_start.source_geometry()
    {
        return Err(DerivedMorphismErrorV1::CompositionClassMismatch);
    }
    match (first_end, second_start) {
        (
            AdmittedDerivedPrimitiveV1::DeclaredChartMap(first_chart),
            AdmittedDerivedPrimitiveV1::DeclaredChartMap(second_chart),
        ) if first_chart.target_chart != second_chart.source_chart => {
            return Err(DerivedMorphismErrorV1::CompositionChartMismatch);
        }
        (
            AdmittedDerivedPrimitiveV1::DeclaredComplexRefinement(first_refinement),
            AdmittedDerivedPrimitiveV1::DeclaredComplexRefinement(second_refinement),
        ) if first_refinement.target_complex != second_refinement.source_complex
            || first_refinement.target_resolution != second_refinement.source_resolution =>
        {
            return Err(DerivedMorphismErrorV1::CompositionRefinementMismatch);
        }
        _ => {}
    }
    Ok(())
}

fn compose_class(
    first: &AdmittedDerivedMorphismV1,
    second: &AdmittedDerivedMorphismV1,
) -> Result<
    (
        AdmittedDerivedMorphismClassV1,
        ReceiptClassV1,
        Option<DerivedChartPathV1>,
    ),
    DerivedMorphismErrorV1,
> {
    validate_typed_primitive_seam(first, second)?;
    match (first.class, second.class) {
        (AdmittedDerivedMorphismClassV1::Strict, AdmittedDerivedMorphismClassV1::Strict) => Ok((
            AdmittedDerivedMorphismClassV1::Strict,
            ReceiptClassV1::CompositeStrict,
            None,
        )),
        (
            AdmittedDerivedMorphismClassV1::DeclaredChartMapPath,
            AdmittedDerivedMorphismClassV1::DeclaredChartMapPath,
        ) => {
            let first_path = first
                .chart_path
                .ok_or(DerivedMorphismErrorV1::CompositionClassMismatch)?;
            let second_path = second
                .chart_path
                .ok_or(DerivedMorphismErrorV1::CompositionClassMismatch)?;
            let first_primitive = first
                .declared_chart_maps
                .last()
                .ok_or(DerivedMorphismErrorV1::CompositionClassMismatch)?;
            let second_primitive = second
                .declared_chart_maps
                .first()
                .ok_or(DerivedMorphismErrorV1::CompositionClassMismatch)?;
            if first_path.target_chart != second_path.source_chart
                || first_primitive.target_geometry != second_primitive.source_geometry
                || first_primitive.target_geometry != first.target
                || second_primitive.source_geometry != second.source
                || first_primitive.target_chart != second_primitive.source_chart
                || first_path.target_chart != first_primitive.target_chart
                || second_path.source_chart != second_primitive.source_chart
            {
                return Err(DerivedMorphismErrorV1::CompositionChartMismatch);
            }
            Ok((
                AdmittedDerivedMorphismClassV1::DeclaredChartMapPath,
                ReceiptClassV1::CompositeDeclaredChartMap,
                Some(DerivedChartPathV1 {
                    source_chart: first_path.source_chart,
                    target_chart: second_path.target_chart,
                }),
            ))
        }
        (
            AdmittedDerivedMorphismClassV1::DeclaredComplexRefinementPath,
            AdmittedDerivedMorphismClassV1::DeclaredComplexRefinementPath,
        ) => Ok((
            AdmittedDerivedMorphismClassV1::DeclaredComplexRefinementPath,
            ReceiptClassV1::CompositeDeclaredComplexRefinement,
            None,
        )),
        (
            AdmittedDerivedMorphismClassV1::DeclaredInclusionPath,
            AdmittedDerivedMorphismClassV1::DeclaredInclusionPath,
        ) => Ok((
            AdmittedDerivedMorphismClassV1::DeclaredInclusionPath,
            ReceiptClassV1::CompositeDeclaredInclusion,
            None,
        )),
        (AdmittedDerivedMorphismClassV1::Identity, _)
        | (_, AdmittedDerivedMorphismClassV1::Identity) => {
            Err(DerivedMorphismErrorV1::CompositionClassMismatch)
        }
        _ => Ok((
            AdmittedDerivedMorphismClassV1::HeterogeneousPath,
            ReceiptClassV1::CompositeHeterogeneous,
            None,
        )),
    }
}

fn compose_evidence(
    first: DerivedEvidenceTransportV1,
    second: DerivedEvidenceTransportV1,
) -> Result<DerivedEvidenceTransportV1, DerivedMorphismErrorV1> {
    match (first, second) {
        (
            DerivedEvidenceTransportV1::BalanceCorestrictionCovariant {
                input_geometry,
                output_geometry: seam_from_first,
                input_evidence,
                output_evidence: seam_evidence_from_first,
                input_rank,
                output_rank: seam_rank_from_first,
            },
            DerivedEvidenceTransportV1::BalanceCorestrictionCovariant {
                input_geometry: seam_from_second,
                output_geometry,
                input_evidence: seam_evidence_from_second,
                output_evidence,
                input_rank: seam_rank_from_second,
                output_rank,
            },
        ) => {
            if seam_from_first != seam_from_second
                || seam_evidence_from_first != seam_evidence_from_second
                || seam_rank_from_first != seam_rank_from_second
            {
                return Err(DerivedMorphismErrorV1::CompositionEvidenceMismatch);
            }
            Ok(DerivedEvidenceTransportV1::BalanceCorestrictionCovariant {
                input_geometry,
                output_geometry,
                input_evidence,
                output_evidence,
                input_rank,
                output_rank,
            })
        }
        (
            DerivedEvidenceTransportV1::RestrictionContravariant {
                input_geometry: seam_from_first,
                output_geometry,
                input_evidence: seam_evidence_from_first,
                output_evidence,
                input_rank: seam_rank_from_first,
                output_rank,
            },
            DerivedEvidenceTransportV1::RestrictionContravariant {
                input_geometry,
                output_geometry: seam_from_second,
                input_evidence,
                output_evidence: seam_evidence_from_second,
                input_rank,
                output_rank: seam_rank_from_second,
            },
        ) => {
            if seam_from_second != seam_from_first
                || seam_evidence_from_second != seam_evidence_from_first
                || seam_rank_from_second != seam_rank_from_first
            {
                return Err(DerivedMorphismErrorV1::CompositionEvidenceMismatch);
            }
            Ok(DerivedEvidenceTransportV1::RestrictionContravariant {
                input_geometry,
                output_geometry,
                input_evidence,
                output_evidence,
                input_rank,
                output_rank,
            })
        }
        _ => Err(DerivedMorphismErrorV1::CompositionVarianceMismatch),
    }
}

/// Compose `first: X -> Y` followed by `second: Y -> Z`.
///
/// Typed primitive factors and no-equivalence artifacts are flattened in
/// semantic order, so parenthesization does not change the receipt. Adjacent
/// declared chart maps require an exact middle chart, and adjacent complex
/// refinements require exact complex and resolution seams, even inside
/// heterogeneous paths. Inclusion declarations compose only through the exact
/// middle geometry and retain factor-local map/containment artifacts without
/// asserting transitive containment. Identity arrows are unique per geometry
/// and rank-neutral, so they are exact composition units.
///
/// # Errors
/// Returns a typed refusal for endpoint, chart/refinement/variance/evidence
/// seams, lineage caps, allocation, cancellation, or identity defects.
#[must_use = "composition refusal must not be treated as a morphism"]
pub fn compose_derived_morphisms_v1(
    first: &AdmittedDerivedMorphismV1,
    second: &AdmittedDerivedMorphismV1,
    cx: &Cx<'_>,
) -> Result<AdmittedDerivedMorphismV1, DerivedMorphismErrorV1> {
    if cx.checkpoint().is_err() {
        return Err(DerivedMorphismErrorV1::Cancelled {
            stage: "composition-entry",
        });
    }
    if first.target != second.source {
        return Err(DerivedMorphismErrorV1::CompositionEndpointMismatch);
    }

    if first.class == AdmittedDerivedMorphismClassV1::Identity {
        return copy_admitted_morphism(second, cx);
    }
    if second.class == AdmittedDerivedMorphismClassV1::Identity {
        return copy_admitted_morphism(first, cx);
    }

    let (class, receipt_class, chart_path) = compose_class(first, second)?;
    let evidence = compose_evidence(first.evidence, second.evidence)?;
    let primitive_path = combine_slices(
        "typed-primitive-lineage",
        &first.primitive_path,
        &second.primitive_path,
        cx,
    )?;
    let primitive_factors = combine_slices(
        "primitive-lineage",
        &first.primitive_factors,
        &second.primitive_factors,
        cx,
    )?;
    let no_equivalence_claims = combine_slices(
        "no-equivalence-claims",
        &first.no_equivalence_claims,
        &second.no_equivalence_claims,
        cx,
    )?;
    let declared_chart_maps = combine_slices(
        "declared-chart-map-lineage",
        &first.declared_chart_maps,
        &second.declared_chart_maps,
        cx,
    )?;
    let receipt = morphism_receipt(
        first.source,
        second.target,
        receipt_class,
        evidence,
        &no_equivalence_claims,
        &primitive_factors,
        cx,
    )?;
    if cx.checkpoint().is_err() {
        return Err(DerivedMorphismErrorV1::Cancelled {
            stage: "composition-publication",
        });
    }
    Ok(AdmittedDerivedMorphismV1 {
        source: first.source,
        target: second.target,
        class,
        chart_path,
        primitive_path,
        declared_chart_maps,
        evidence,
        no_equivalence_claims,
        primitive_factors,
        receipt,
    })
}

fn validate_stratum_object(
    object: DerivedStratumObjectV1,
    geometry: &AdmittedDerivedGeometryV1,
    geometry_field: &'static str,
    stratification_field: &'static str,
    stratum_field: &'static str,
    cx: &Cx<'_>,
) -> Result<(), DerivedStratumMorphismErrorV1> {
    if object.geometry != geometry.id() {
        return Err(DerivedStratumMorphismErrorV1::EndpointMismatch {
            field: geometry_field,
        });
    }
    if is_zero(object.stratification.as_bytes()) {
        return Err(DerivedStratumMorphismErrorV1::MissingIdentity {
            field: stratification_field,
        });
    }
    if object.stratification != geometry.ir().stratification.id {
        return Err(DerivedStratumMorphismErrorV1::MissingStratification {
            field: stratification_field,
        });
    }
    if is_zero(object.stratum.as_bytes()) {
        return Err(DerivedStratumMorphismErrorV1::MissingIdentity {
            field: stratum_field,
        });
    }
    let mut owns_stratum = false;
    for (index, stratum) in geometry.ir().stratification.strata.iter().enumerate() {
        if index.is_multiple_of(DERIVED_MORPHISM_CANCELLATION_STRIDE_V1) && cx.checkpoint().is_err()
        {
            return Err(DerivedStratumMorphismErrorV1::Cancelled {
                stage: "stratum-selector-resolution",
            });
        }
        if stratum.id == object.stratum {
            owns_stratum = true;
            break;
        }
    }
    if !owns_stratum {
        return Err(DerivedStratumMorphismErrorV1::MissingStratum {
            field: stratum_field,
        });
    }
    Ok(())
}

fn validate_stratum_nonidentity_compatibility(
    source: &AdmittedDerivedGeometryV1,
    target: &AdmittedDerivedGeometryV1,
) -> Result<(), DerivedStratumMorphismErrorV1> {
    if source.ir().subject != target.ir().subject {
        return Err(DerivedStratumMorphismErrorV1::SubjectMismatch);
    }
    if source.ir().model_version != target.ir().model_version {
        return Err(DerivedStratumMorphismErrorV1::ModelVersionMismatch);
    }
    if source.ir().category != target.ir().category {
        return Err(DerivedStratumMorphismErrorV1::CategoryMismatch);
    }
    if source.ir().coefficients != target.ir().coefficients {
        return Err(DerivedStratumMorphismErrorV1::CoefficientMismatch);
    }
    if source.ir().frame != target.ir().frame {
        return Err(DerivedStratumMorphismErrorV1::FrameMismatch);
    }
    if source.ir().unit_system != target.ir().unit_system {
        return Err(DerivedStratumMorphismErrorV1::UnitSystemMismatch);
    }
    Ok(())
}

fn retain_stratum_value<T>(
    field: &'static str,
    value: Option<T>,
) -> Result<Vec<T>, DerivedStratumMorphismErrorV1> {
    let mut retained = Vec::new();
    if let Some(value) = value {
        retained
            .try_reserve_exact(1)
            .map_err(|_| DerivedStratumMorphismErrorV1::AllocationRefused { field })?;
        retained.push(value);
    }
    Ok(retained)
}

/// Admit one identity or nominal component in the stratum-scoped category.
///
/// Both selectors are resolved against their exact supplied admitted geometry.
/// A component may connect strata with different charts or dimensions, but its
/// two geometries must retain the same subject, immutable model version,
/// category, coefficient, frame, and unit semantics. The result is not a
/// whole-geometry morphism and carries no evidence-transport capability.
///
/// # Errors
/// Returns a typed refusal for schema, endpoint ownership, compatibility,
/// authority, allocation, cancellation, or canonical-identity defects.
#[must_use = "a raw stratum-map declaration has no structural authority"]
#[allow(clippy::too_many_lines)] // One exhaustive identity/component admission dispatch.
pub fn admit_derived_stratum_morphism_v1(
    ir: &DerivedStratumMorphismIrV1,
    source_geometry: &AdmittedDerivedGeometryV1,
    target_geometry: &AdmittedDerivedGeometryV1,
    cx: &Cx<'_>,
) -> Result<AdmittedDerivedStratumMorphismV1, DerivedStratumMorphismErrorV1> {
    if cx.checkpoint().is_err() {
        return Err(DerivedStratumMorphismErrorV1::Cancelled {
            stage: "stratum-admission-entry",
        });
    }
    if ir.schema_version != DERIVED_STRATUM_MORPHISM_SCHEMA_VERSION_V1 {
        return Err(DerivedStratumMorphismErrorV1::UnsupportedSchemaVersion {
            found: ir.schema_version,
            supported: DERIVED_STRATUM_MORPHISM_SCHEMA_VERSION_V1,
        });
    }
    validate_stratum_object(
        ir.source,
        source_geometry,
        "source-geometry",
        "source-stratification",
        "source-stratum",
        cx,
    )?;
    validate_stratum_object(
        ir.target,
        target_geometry,
        "target-geometry",
        "target-stratification",
        "target-stratum",
        cx,
    )?;

    let (class, receipt_class, no_claim, primitive) = match (ir.kind, ir.authority) {
        (
            DerivedStratumMorphismKindV1::Identity,
            DerivedStratumAuthorityBoundaryV1::IdentityOnly,
        ) => {
            if ir.source != ir.target {
                return Err(DerivedStratumMorphismErrorV1::InvalidIdentity);
            }
            (
                AdmittedDerivedStratumMorphismClassV1::Identity,
                StratumReceiptClassV1::Identity,
                None,
                None,
            )
        }
        (
            DerivedStratumMorphismKindV1::DeclaredComponent {
                map,
                constructibility,
            },
            DerivedStratumAuthorityBoundaryV1::NoClaim { artifact },
        ) => {
            for (bytes, field) in [
                (map.as_bytes(), "stratum-map"),
                (constructibility.as_bytes(), "constructibility-declaration"),
                (artifact.as_bytes(), "no-stratum-map-authority"),
            ] {
                if is_zero(bytes) {
                    return Err(DerivedStratumMorphismErrorV1::MissingIdentity { field });
                }
            }
            validate_stratum_nonidentity_compatibility(source_geometry, target_geometry)?;
            (
                AdmittedDerivedStratumMorphismClassV1::DeclaredPath,
                StratumReceiptClassV1::PrimitiveDeclaredComponent {
                    map,
                    constructibility,
                },
                Some(artifact),
                Some(DeclaredStratumMapPrimitiveV1 {
                    source: ir.source,
                    target: ir.target,
                    map,
                    constructibility,
                }),
            )
        }
        (DerivedStratumMorphismKindV1::Identity, _) => {
            return Err(DerivedStratumMorphismErrorV1::InvalidIdentity);
        }
        (DerivedStratumMorphismKindV1::DeclaredComponent { .. }, _) => {
            return Err(DerivedStratumMorphismErrorV1::AuthorityLaundering);
        }
    };

    if cx.checkpoint().is_err() {
        return Err(DerivedStratumMorphismErrorV1::Cancelled {
            stage: "stratum-admission",
        });
    }
    let no_authority_claims = retain_stratum_value("stratum-no-authority-claims", no_claim)?;
    let primitive_path = retain_stratum_value("stratum-typed-primitive-lineage", primitive)?;
    let receipt = stratum_morphism_receipt(
        ir.source,
        ir.target,
        receipt_class,
        &no_authority_claims,
        &[],
        cx,
    )?;
    let mut primitive_factors = Vec::new();
    if class == AdmittedDerivedStratumMorphismClassV1::DeclaredPath {
        primitive_factors.try_reserve_exact(1).map_err(|_| {
            DerivedStratumMorphismErrorV1::AllocationRefused {
                field: "stratum-primitive-lineage",
            }
        })?;
        primitive_factors.push(receipt.id());
    }
    if cx.checkpoint().is_err() {
        return Err(DerivedStratumMorphismErrorV1::Cancelled {
            stage: "stratum-publication",
        });
    }
    Ok(AdmittedDerivedStratumMorphismV1 {
        source: ir.source,
        target: ir.target,
        class,
        primitive_path,
        no_authority_claims,
        primitive_factors,
        receipt,
    })
}

/// Mint the exact identity on one stratum of an admitted geometry.
///
/// # Errors
/// Returns a typed refusal if the stratum is not owned by the geometry, or if
/// cancellation, allocation, or canonical identity construction fails.
#[must_use = "identity construction must complete before composition"]
pub fn identity_derived_stratum_morphism_v1(
    geometry: &AdmittedDerivedGeometryV1,
    stratum: StratumIdV1,
    cx: &Cx<'_>,
) -> Result<AdmittedDerivedStratumMorphismV1, DerivedStratumMorphismErrorV1> {
    let object = DerivedStratumObjectV1 {
        geometry: geometry.id(),
        stratification: geometry.ir().stratification.id,
        stratum,
    };
    admit_derived_stratum_morphism_v1(
        &DerivedStratumMorphismIrV1 {
            schema_version: DERIVED_STRATUM_MORPHISM_SCHEMA_VERSION_V1,
            source: object,
            target: object,
            kind: DerivedStratumMorphismKindV1::Identity,
            authority: DerivedStratumAuthorityBoundaryV1::IdentityOnly,
        },
        geometry,
        geometry,
        cx,
    )
}

fn checked_stratum_combined_len(
    field: &'static str,
    left: usize,
    right: usize,
) -> Result<usize, DerivedStratumMorphismErrorV1> {
    let requested =
        left.checked_add(right)
            .ok_or(DerivedStratumMorphismErrorV1::ResourceLimit {
                field,
                requested: usize::MAX,
                limit: DERIVED_MORPHISM_MAX_FACTORS_V1,
            })?;
    if requested > DERIVED_MORPHISM_MAX_FACTORS_V1 {
        return Err(DerivedStratumMorphismErrorV1::ResourceLimit {
            field,
            requested,
            limit: DERIVED_MORPHISM_MAX_FACTORS_V1,
        });
    }
    Ok(requested)
}

fn combine_stratum_slices<T: Copy>(
    field: &'static str,
    left: &[T],
    right: &[T],
    cx: &Cx<'_>,
) -> Result<Vec<T>, DerivedStratumMorphismErrorV1> {
    let len = checked_stratum_combined_len(field, left.len(), right.len())?;
    if cx.checkpoint().is_err() {
        return Err(DerivedStratumMorphismErrorV1::Cancelled { stage: field });
    }
    let mut out = Vec::new();
    out.try_reserve_exact(len)
        .map_err(|_| DerivedStratumMorphismErrorV1::AllocationRefused { field })?;
    for (index, value) in left.iter().chain(right).copied().enumerate() {
        out.push(value);
        if (index + 1).is_multiple_of(DERIVED_MORPHISM_CANCELLATION_STRIDE_V1)
            && cx.checkpoint().is_err()
        {
            return Err(DerivedStratumMorphismErrorV1::Cancelled { stage: field });
        }
    }
    if cx.checkpoint().is_err() {
        return Err(DerivedStratumMorphismErrorV1::Cancelled { stage: field });
    }
    Ok(out)
}

fn validate_admitted_stratum_path(
    value: &AdmittedDerivedStratumMorphismV1,
    cx: &Cx<'_>,
) -> Result<(), DerivedStratumMorphismErrorV1> {
    match value.class {
        AdmittedDerivedStratumMorphismClassV1::Identity => {
            if value.source != value.target
                || !value.primitive_path.is_empty()
                || !value.no_authority_claims.is_empty()
                || !value.primitive_factors.is_empty()
            {
                return Err(DerivedStratumMorphismErrorV1::CompositionClassMismatch);
            }
        }
        AdmittedDerivedStratumMorphismClassV1::DeclaredPath => {
            if value.primitive_path.is_empty()
                || value.primitive_path.len() != value.primitive_factors.len()
                || value.primitive_path.len() != value.no_authority_claims.len()
                || value
                    .primitive_path
                    .first()
                    .map(|primitive| primitive.source)
                    != Some(value.source)
                || value
                    .primitive_path
                    .last()
                    .map(|primitive| primitive.target)
                    != Some(value.target)
            {
                return Err(DerivedStratumMorphismErrorV1::CompositionClassMismatch);
            }
            for (index, pair) in value.primitive_path.windows(2).enumerate() {
                if index.is_multiple_of(DERIVED_MORPHISM_CANCELLATION_STRIDE_V1)
                    && cx.checkpoint().is_err()
                {
                    return Err(DerivedStratumMorphismErrorV1::Cancelled {
                        stage: "stratum-path-validation",
                    });
                }
                if pair[0].target != pair[1].source {
                    return Err(DerivedStratumMorphismErrorV1::CompositionClassMismatch);
                }
            }
        }
    }
    Ok(())
}

fn copy_admitted_stratum_morphism(
    value: &AdmittedDerivedStratumMorphismV1,
    cx: &Cx<'_>,
) -> Result<AdmittedDerivedStratumMorphismV1, DerivedStratumMorphismErrorV1> {
    validate_admitted_stratum_path(value, cx)?;
    let primitive_path = combine_stratum_slices(
        "stratum-typed-primitive-lineage-copy",
        &value.primitive_path,
        &[],
        cx,
    )?;
    let no_authority_claims = combine_stratum_slices(
        "stratum-no-authority-claims-copy",
        &value.no_authority_claims,
        &[],
        cx,
    )?;
    let primitive_factors = combine_stratum_slices(
        "stratum-primitive-lineage-copy",
        &value.primitive_factors,
        &[],
        cx,
    )?;
    Ok(AdmittedDerivedStratumMorphismV1 {
        source: value.source,
        target: value.target,
        class: value.class,
        primitive_path,
        no_authority_claims,
        primitive_factors,
        receipt: value.receipt,
    })
}

/// Compose `first: S -> T` followed by `second: T -> U` in the standalone
/// stratum-scoped category.
///
/// The complete middle `(geometry, stratification, stratum)` object must match.
/// Composition flattens component and no-authority lineage associatively. It
/// cannot compose with `AdmittedDerivedMorphismV1` and creates no evidence
/// transport or whole-stratification map.
///
/// # Errors
/// Returns a typed refusal for an inexact seam, inconsistent sealed class,
/// lineage cap, allocation, cancellation, or canonical identity defect.
#[must_use = "composition refusal must not be treated as a stratum morphism"]
pub fn compose_derived_stratum_morphisms_v1(
    first: &AdmittedDerivedStratumMorphismV1,
    second: &AdmittedDerivedStratumMorphismV1,
    cx: &Cx<'_>,
) -> Result<AdmittedDerivedStratumMorphismV1, DerivedStratumMorphismErrorV1> {
    if cx.checkpoint().is_err() {
        return Err(DerivedStratumMorphismErrorV1::Cancelled {
            stage: "stratum-composition-entry",
        });
    }
    validate_admitted_stratum_path(first, cx)?;
    validate_admitted_stratum_path(second, cx)?;
    if first.target != second.source {
        return Err(DerivedStratumMorphismErrorV1::CompositionEndpointMismatch);
    }
    if first.class == AdmittedDerivedStratumMorphismClassV1::Identity {
        return copy_admitted_stratum_morphism(second, cx);
    }
    if second.class == AdmittedDerivedStratumMorphismClassV1::Identity {
        return copy_admitted_stratum_morphism(first, cx);
    }

    let primitive_path = combine_stratum_slices(
        "stratum-typed-primitive-lineage",
        &first.primitive_path,
        &second.primitive_path,
        cx,
    )?;
    let no_authority_claims = combine_stratum_slices(
        "stratum-no-authority-claims",
        &first.no_authority_claims,
        &second.no_authority_claims,
        cx,
    )?;
    let primitive_factors = combine_stratum_slices(
        "stratum-primitive-lineage",
        &first.primitive_factors,
        &second.primitive_factors,
        cx,
    )?;
    let receipt = stratum_morphism_receipt(
        first.source,
        second.target,
        StratumReceiptClassV1::CompositeDeclaredPath,
        &no_authority_claims,
        &primitive_factors,
        cx,
    )?;
    if cx.checkpoint().is_err() {
        return Err(DerivedStratumMorphismErrorV1::Cancelled {
            stage: "stratum-composition-publication",
        });
    }
    Ok(AdmittedDerivedStratumMorphismV1 {
        source: first.source,
        target: second.target,
        class: AdmittedDerivedStratumMorphismClassV1::DeclaredPath,
        primitive_path,
        no_authority_claims,
        primitive_factors,
        receipt,
    })
}

/// Admit exhaustive finite source-stratum coverage by direct sealed components.
///
/// Raw bindings and supplied tokens must follow the admitted source
/// stratification's canonical stratum order. Exactly one direct component is
/// required for each source stratum; target strata may repeat and need not be
/// exhaustive. An exact stratum identity is direct; a declared path is direct
/// only when it retains exactly one primitive. The nominal assembly and global
/// constructibility identities are retained but not resolved or authenticated.
/// Success creates a standalone candidate with no whole-geometry map, evidence
/// transport, composition, continuity, incidence/frontier, gluing,
/// constructibility, or equivalence authority.
///
/// # Errors
/// Returns a typed refusal for schema, endpoint/stratification binding, missing
/// nominal IDs, incomplete or misordered source coverage, raw/sealed component
/// mismatch, composite components, resource/allocation limits, cancellation, or
/// canonical identity defects.
#[must_use = "an exhaustive component packet has no global map authority"]
#[allow(clippy::too_many_lines)] // One bounded exhaustive coverage admission scan.
pub fn admit_derived_exhaustive_stratified_map_candidate_v1(
    ir: &DerivedExhaustiveStratifiedMapCandidateIrV1,
    source: &AdmittedDerivedGeometryV1,
    target: &AdmittedDerivedGeometryV1,
    components: &[AdmittedDerivedStratumMorphismV1],
    cx: &Cx<'_>,
) -> Result<
    AdmittedDerivedExhaustiveStratifiedMapCandidateV1,
    DerivedExhaustiveStratifiedMapCandidateErrorV1,
> {
    if cx.checkpoint().is_err() {
        return Err(DerivedExhaustiveStratifiedMapCandidateErrorV1::Cancelled {
            stage: "stratified-assembly-admission-entry",
        });
    }
    if ir.schema_version != DERIVED_EXHAUSTIVE_STRATIFIED_MAP_CANDIDATE_SCHEMA_VERSION_V1 {
        return Err(
            DerivedExhaustiveStratifiedMapCandidateErrorV1::UnsupportedSchemaVersion {
                found: ir.schema_version,
                supported: DERIVED_EXHAUSTIVE_STRATIFIED_MAP_CANDIDATE_SCHEMA_VERSION_V1,
            },
        );
    }
    for (matches, field) in [
        (ir.source_geometry == source.id(), "source-geometry"),
        (ir.target_geometry == target.id(), "target-geometry"),
    ] {
        if !matches {
            return Err(DerivedExhaustiveStratifiedMapCandidateErrorV1::EndpointMismatch { field });
        }
    }
    for (bytes, field) in [
        (ir.source_stratification.as_bytes(), "source-stratification"),
        (ir.target_stratification.as_bytes(), "target-stratification"),
        (ir.nominal_assembly.as_bytes(), "nominal-assembly"),
        (
            ir.nominal_constructibility.as_bytes(),
            "nominal-global-constructibility",
        ),
        (ir.no_authority.as_bytes(), "no-global-map-authority"),
    ] {
        if is_zero(bytes) {
            return Err(DerivedExhaustiveStratifiedMapCandidateErrorV1::MissingIdentity { field });
        }
    }
    for (matches, field) in [
        (
            ir.source_stratification == source.ir().stratification.id,
            "source-stratification",
        ),
        (
            ir.target_stratification == target.ir().stratification.id,
            "target-stratification",
        ),
    ] {
        if !matches {
            return Err(
                DerivedExhaustiveStratifiedMapCandidateErrorV1::StratificationMismatch { field },
            );
        }
    }

    let expected = source.ir().stratification.strata.len();
    for (field, found) in [
        ("component-bindings", ir.components.len()),
        ("sealed-components", components.len()),
    ] {
        if found > DERIVED_MORPHISM_MAX_FACTORS_V1 {
            return Err(
                DerivedExhaustiveStratifiedMapCandidateErrorV1::ResourceLimit {
                    field,
                    requested: found,
                    limit: DERIVED_MORPHISM_MAX_FACTORS_V1,
                },
            );
        }
        if found != expected {
            return Err(
                DerivedExhaustiveStratifiedMapCandidateErrorV1::ComponentCountMismatch {
                    field,
                    expected,
                    found,
                },
            );
        }
    }
    let mut retained = Vec::new();
    retained.try_reserve_exact(expected).map_err(|_| {
        DerivedExhaustiveStratifiedMapCandidateErrorV1::AllocationRefused {
            field: "component-bindings",
        }
    })?;
    for (index, ((source_stratum, binding), component)) in source
        .ir()
        .stratification
        .strata
        .iter()
        .zip(&ir.components)
        .zip(components)
        .enumerate()
    {
        if index.is_multiple_of(DERIVED_MORPHISM_CANCELLATION_STRIDE_V1) && cx.checkpoint().is_err()
        {
            return Err(DerivedExhaustiveStratifiedMapCandidateErrorV1::Cancelled {
                stage: "stratified-assembly-components",
            });
        }
        if binding.component != component.id() {
            return Err(
                DerivedExhaustiveStratifiedMapCandidateErrorV1::ComponentIdentityMismatch { index },
            );
        }
        for (matches, field) in [
            (
                binding.source_stratum == source_stratum.id,
                "binding-source-stratum-order",
            ),
            (
                component.source().geometry == ir.source_geometry,
                "component-source-geometry",
            ),
            (
                component.source().stratification == ir.source_stratification,
                "component-source-stratification",
            ),
            (
                component.source().stratum == binding.source_stratum,
                "component-source-stratum",
            ),
            (
                component.target().geometry == ir.target_geometry,
                "component-target-geometry",
            ),
            (
                component.target().stratification == ir.target_stratification,
                "component-target-stratification",
            ),
            (
                component.target().stratum == binding.target_stratum,
                "component-target-stratum",
            ),
        ] {
            if !matches {
                return Err(
                    DerivedExhaustiveStratifiedMapCandidateErrorV1::ComponentEndpointMismatch {
                        index,
                        field,
                    },
                );
            }
        }
        match component.class() {
            AdmittedDerivedStratumMorphismClassV1::Identity => {}
            AdmittedDerivedStratumMorphismClassV1::DeclaredPath
                if component.primitive_path().len() == 1 => {}
            AdmittedDerivedStratumMorphismClassV1::DeclaredPath => {
                return Err(
                    DerivedExhaustiveStratifiedMapCandidateErrorV1::CompositeComponent { index },
                );
            }
        }
        retained.push(*binding);
    }
    if cx.checkpoint().is_err() {
        return Err(DerivedExhaustiveStratifiedMapCandidateErrorV1::Cancelled {
            stage: "stratified-assembly-admission",
        });
    }
    let receipt = exhaustive_stratified_map_candidate_receipt(ir, &retained, cx)?;
    if cx.checkpoint().is_err() {
        return Err(DerivedExhaustiveStratifiedMapCandidateErrorV1::Cancelled {
            stage: "stratified-assembly-publication",
        });
    }
    Ok(AdmittedDerivedExhaustiveStratifiedMapCandidateV1 {
        source_geometry: ir.source_geometry,
        source_stratification: ir.source_stratification,
        target_geometry: ir.target_geometry,
        target_stratification: ir.target_stratification,
        components: retained,
        nominal_assembly: ir.nominal_assembly,
        nominal_constructibility: ir.nominal_constructibility,
        no_authority: ir.no_authority,
        receipt,
    })
}

fn span_correspondence_receipt(
    ir: DerivedSpanCorrespondenceIrV1,
    cx: &Cx<'_>,
) -> Result<IdentityReceipt<DerivedSpanCorrespondenceIdV1>, DerivedSpanCorrespondenceErrorV1> {
    if cx.checkpoint().is_err() {
        return Err(DerivedSpanCorrespondenceErrorV1::Cancelled {
            stage: "span-identity-entry",
        });
    }
    let map_identity_error = |error| match error {
        CanonicalError::Cancelled { .. } => DerivedSpanCorrespondenceErrorV1::Cancelled {
            stage: "span-identity",
        },
        other => DerivedSpanCorrespondenceErrorV1::Identity(other),
    };
    CanonicalEncoder::<DerivedSpanCorrespondenceIdV1, _>::new(
        DERIVED_MORPHISM_IDENTITY_LIMITS_V1,
        || cx.checkpoint().is_err(),
    )
    .map_err(map_identity_error)?
    .bytes(Field::new(0, "source"), ir.source.as_bytes())
    .and_then(|encoder| encoder.bytes(Field::new(1, "apex"), ir.apex.as_bytes()))
    .and_then(|encoder| encoder.bytes(Field::new(2, "target"), ir.target.as_bytes()))
    .and_then(|encoder| encoder.bytes(Field::new(3, "left-leg"), ir.left_leg.as_bytes()))
    .and_then(|encoder| encoder.bytes(Field::new(4, "right-leg"), ir.right_leg.as_bytes()))
    .and_then(|encoder| encoder.bytes(Field::new(5, "no-claim"), ir.no_claim.as_bytes()))
    .and_then(|encoder| encoder.finish())
    .map_err(map_identity_error)
}

/// Admit one standalone declared span `source <- apex -> target`.
///
/// Both supplied legs must already be sealed RD.1b morphisms. Admission binds
/// their exact identities and checks orientations `apex -> source` and
/// `apex -> target`. It does not create a direct source-to-target map or
/// evidence transport, authenticate leg payloads beyond their sealed structural
/// receipts, or prove functionality, pullback, equivalence, or physical
/// correspondence. Span composition requires separately admitted pullback data
/// and is intentionally absent from v1.
///
/// # Errors
/// Returns a typed refusal for schema, no-claim identity, raw-leg binding,
/// sealed-leg orientation, cancellation, or canonical-identity defects.
#[must_use = "a raw span request has no correspondence authority"]
pub fn admit_derived_span_correspondence_v1(
    ir: DerivedSpanCorrespondenceIrV1,
    left_leg: &AdmittedDerivedMorphismV1,
    right_leg: &AdmittedDerivedMorphismV1,
    cx: &Cx<'_>,
) -> Result<AdmittedDerivedSpanCorrespondenceV1, DerivedSpanCorrespondenceErrorV1> {
    if cx.checkpoint().is_err() {
        return Err(DerivedSpanCorrespondenceErrorV1::Cancelled {
            stage: "span-admission-entry",
        });
    }
    if ir.schema_version != DERIVED_SPAN_CORRESPONDENCE_SCHEMA_VERSION_V1 {
        return Err(DerivedSpanCorrespondenceErrorV1::UnsupportedSchemaVersion {
            found: ir.schema_version,
            supported: DERIVED_SPAN_CORRESPONDENCE_SCHEMA_VERSION_V1,
        });
    }
    if is_zero(ir.no_claim.as_bytes()) {
        return Err(DerivedSpanCorrespondenceErrorV1::MissingIdentity {
            field: "no-correspondence-claim",
        });
    }
    for (matches, field) in [
        (ir.left_leg == left_leg.id(), "left-leg"),
        (ir.right_leg == right_leg.id(), "right-leg"),
    ] {
        if !matches {
            return Err(DerivedSpanCorrespondenceErrorV1::LegIdentityMismatch { field });
        }
    }
    for (matches, field) in [
        (left_leg.source() == ir.apex, "left-source-apex"),
        (left_leg.target() == ir.source, "left-target-source"),
        (right_leg.source() == ir.apex, "right-source-apex"),
        (right_leg.target() == ir.target, "right-target-target"),
    ] {
        if !matches {
            return Err(DerivedSpanCorrespondenceErrorV1::LegOrientationMismatch { field });
        }
    }
    let receipt = span_correspondence_receipt(ir, cx)?;
    if cx.checkpoint().is_err() {
        return Err(DerivedSpanCorrespondenceErrorV1::Cancelled {
            stage: "span-publication",
        });
    }
    Ok(AdmittedDerivedSpanCorrespondenceV1 {
        source: ir.source,
        apex: ir.apex,
        target: ir.target,
        left_leg: ir.left_leg,
        right_leg: ir.right_leg,
        no_claim: ir.no_claim,
        receipt,
    })
}

#[derive(Debug, Clone, Copy)]
struct FixedResolutionQuasiIsomorphismCandidateBindingV1 {
    refinement_path: DerivedMorphismIdV1,
    source_geometry: DerivedGeometryIdV1,
    target_geometry: DerivedGeometryIdV1,
    source_local_model: DerivedLocalModelIdV1,
    target_local_model: DerivedLocalModelIdV1,
    complex_role: DerivedComplexRoleV1,
    source_complex: DerivedComplexIdV1,
    target_complex: DerivedComplexIdV1,
    source_resolution: DerivedResolutionIdV1,
    target_resolution: DerivedResolutionIdV1,
    source_scope_witness: DerivedWitnessIdV1,
    target_scope_witness: DerivedWitnessIdV1,
    nominal_theorem: DerivedTheoremIdV1,
    nominal_checker: DerivedCheckerIdV1,
    nominal_check_receipt: DerivedWitnessIdV1,
    no_authority: DerivedNoClaimIdV1,
}

const fn complex_role_tag(role: DerivedComplexRoleV1) -> [u8; 1] {
    match role {
        DerivedComplexRoleV1::Tangent => [0],
        DerivedComplexRoleV1::Cotangent => [1],
        DerivedComplexRoleV1::DeformationObstruction => [2],
    }
}

const fn local_model_complex_for_role(
    model: &DerivedLocalModelV1,
    role: DerivedComplexRoleV1,
) -> DerivedComplexIdV1 {
    match role {
        DerivedComplexRoleV1::Tangent => model.tangent_complex,
        DerivedComplexRoleV1::Cotangent => model.cotangent_complex,
        DerivedComplexRoleV1::DeformationObstruction => model.deformation_complex,
    }
}

fn fixed_resolution_scope_witness(
    model: &DerivedLocalModelV1,
    expected_resolution: DerivedResolutionIdV1,
    presentation_field: &'static str,
    resolution_field: &'static str,
    witness_field: &'static str,
) -> Result<DerivedWitnessIdV1, DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1> {
    match model.presentation {
        PresentationScopeV1::FixedResolution {
            resolution,
            witness,
        } => {
            if resolution != expected_resolution {
                return Err(
                    DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1::ResolutionScopeMismatch {
                        field: resolution_field,
                    },
                );
            }
            if is_zero(witness.as_bytes()) {
                return Err(
                    DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1::MissingIdentity {
                        field: witness_field,
                    },
                );
            }
            Ok(witness)
        }
        PresentationScopeV1::Literal { .. } | PresentationScopeV1::ExternallyChecked { .. } => Err(
            DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1::PresentationScopeMismatch {
                field: presentation_field,
            },
        ),
    }
}

fn validate_fixed_resolution_candidate_path(
    ir: &DerivedFixedResolutionQuasiIsomorphismCandidateIrV1,
    refinement_path: &AdmittedDerivedMorphismV1,
    cx: &Cx<'_>,
) -> Result<(), DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1> {
    if refinement_path.class() != AdmittedDerivedMorphismClassV1::DeclaredComplexRefinementPath {
        return Err(
            DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1::PathClassMismatch {
                found: refinement_path.class(),
            },
        );
    }
    let primitives = refinement_path.primitive_path();
    if primitives.is_empty() {
        return Err(
            DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1::PathShapeMismatch {
                field: "empty-refinement-path",
                index: 0,
            },
        );
    }

    let mut previous_target = None;
    for (index, primitive) in primitives.iter().enumerate() {
        if index.is_multiple_of(DERIVED_MORPHISM_CANCELLATION_STRIDE_V1) && cx.checkpoint().is_err()
        {
            return Err(
                DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1::Cancelled {
                    stage: "candidate-path-scan",
                },
            );
        }
        let AdmittedDerivedPrimitiveV1::DeclaredComplexRefinement(refinement) = primitive else {
            return Err(
                DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1::PathShapeMismatch {
                    field: "non-refinement-primitive",
                    index,
                },
            );
        };
        if let Some((geometry, complex, resolution)) = previous_target
            && (
                refinement.source_geometry,
                refinement.source_complex,
                refinement.source_resolution,
            ) != (geometry, complex, resolution)
        {
            return Err(
                DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1::PathShapeMismatch {
                    field: "refinement-seam",
                    index,
                },
            );
        }
        if index == 0
            && (
                refinement.source_geometry,
                refinement.source_complex,
                refinement.source_resolution,
            ) != (ir.source_geometry, ir.source_complex, ir.source_resolution)
        {
            return Err(
                DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1::PathShapeMismatch {
                    field: "source-selector",
                    index,
                },
            );
        }
        if index + 1 == primitives.len()
            && (
                refinement.target_geometry,
                refinement.target_complex,
                refinement.target_resolution,
            ) != (ir.target_geometry, ir.target_complex, ir.target_resolution)
        {
            return Err(
                DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1::PathShapeMismatch {
                    field: "target-selector",
                    index,
                },
            );
        }
        previous_target = Some((
            refinement.target_geometry,
            refinement.target_complex,
            refinement.target_resolution,
        ));
    }
    Ok(())
}

fn fixed_resolution_quasi_isomorphism_candidate_receipt(
    binding: &FixedResolutionQuasiIsomorphismCandidateBindingV1,
    cx: &Cx<'_>,
) -> Result<
    IdentityReceipt<DerivedFixedResolutionQuasiIsomorphismCandidateIdV1>,
    DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1,
> {
    if cx.checkpoint().is_err() {
        return Err(
            DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1::Cancelled {
                stage: "candidate-identity-entry",
            },
        );
    }
    let role = complex_role_tag(binding.complex_role);
    let map_identity_error = |error| match error {
        CanonicalError::Cancelled { .. } => {
            DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1::Cancelled {
                stage: "candidate-identity",
            }
        }
        other => DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1::Identity(other),
    };
    CanonicalEncoder::<DerivedFixedResolutionQuasiIsomorphismCandidateIdV1, _>::new(
        DERIVED_FIXED_RESOLUTION_QUASI_ISOMORPHISM_CANDIDATE_IDENTITY_LIMITS_V1,
        || cx.checkpoint().is_err(),
    )
    .map_err(map_identity_error)?
    .bytes(
        Field::new(0, "refinement-path"),
        binding.refinement_path.as_bytes(),
    )
    .and_then(|encoder| {
        encoder.bytes(
            Field::new(1, "source-geometry"),
            binding.source_geometry.as_bytes(),
        )
    })
    .and_then(|encoder| {
        encoder.bytes(
            Field::new(2, "target-geometry"),
            binding.target_geometry.as_bytes(),
        )
    })
    .and_then(|encoder| {
        encoder.bytes(
            Field::new(3, "source-local-model"),
            binding.source_local_model.as_bytes(),
        )
    })
    .and_then(|encoder| {
        encoder.bytes(
            Field::new(4, "target-local-model"),
            binding.target_local_model.as_bytes(),
        )
    })
    .and_then(|encoder| encoder.bytes(Field::new(5, "complex-role"), &role))
    .and_then(|encoder| {
        encoder.bytes(
            Field::new(6, "source-complex"),
            binding.source_complex.as_bytes(),
        )
    })
    .and_then(|encoder| {
        encoder.bytes(
            Field::new(7, "target-complex"),
            binding.target_complex.as_bytes(),
        )
    })
    .and_then(|encoder| {
        encoder.bytes(
            Field::new(8, "source-resolution"),
            binding.source_resolution.as_bytes(),
        )
    })
    .and_then(|encoder| {
        encoder.bytes(
            Field::new(9, "target-resolution"),
            binding.target_resolution.as_bytes(),
        )
    })
    .and_then(|encoder| {
        encoder.bytes(
            Field::new(10, "source-scope-witness"),
            binding.source_scope_witness.as_bytes(),
        )
    })
    .and_then(|encoder| {
        encoder.bytes(
            Field::new(11, "target-scope-witness"),
            binding.target_scope_witness.as_bytes(),
        )
    })
    .and_then(|encoder| {
        encoder.bytes(
            Field::new(12, "nominal-theorem"),
            binding.nominal_theorem.as_bytes(),
        )
    })
    .and_then(|encoder| {
        encoder.bytes(
            Field::new(13, "nominal-checker"),
            binding.nominal_checker.as_bytes(),
        )
    })
    .and_then(|encoder| {
        encoder.bytes(
            Field::new(14, "nominal-check-receipt"),
            binding.nominal_check_receipt.as_bytes(),
        )
    })
    .and_then(|encoder| {
        encoder.bytes(
            Field::new(15, "no-authority"),
            binding.no_authority.as_bytes(),
        )
    })
    .and_then(|encoder| encoder.finish())
    .map_err(map_identity_error)
}

#[allow(clippy::too_many_lines)] // One fail-closed structural candidate admission.
fn fixed_resolution_quasi_isomorphism_candidate_binding(
    ir: &DerivedFixedResolutionQuasiIsomorphismCandidateIrV1,
    source: &AdmittedDerivedGeometryV1,
    target: &AdmittedDerivedGeometryV1,
    refinement_path: &AdmittedDerivedMorphismV1,
    cx: &Cx<'_>,
) -> Result<
    FixedResolutionQuasiIsomorphismCandidateBindingV1,
    DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1,
> {
    if ir.source_geometry != source.id() {
        return Err(
            DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1::EndpointMismatch {
                field: "source-geometry",
            },
        );
    }
    if ir.target_geometry != target.id() {
        return Err(
            DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1::EndpointMismatch {
                field: "target-geometry",
            },
        );
    }
    for (bytes, field) in [
        (ir.source_local_model.as_bytes(), "source-local-model"),
        (ir.target_local_model.as_bytes(), "target-local-model"),
        (ir.source_complex.as_bytes(), "source-complex"),
        (ir.target_complex.as_bytes(), "target-complex"),
        (ir.source_resolution.as_bytes(), "source-resolution"),
        (ir.target_resolution.as_bytes(), "target-resolution"),
        (ir.refinement_path.as_bytes(), "refinement-path"),
        (ir.nominal_theorem.as_bytes(), "nominal-theorem"),
        (ir.nominal_checker.as_bytes(), "nominal-checker"),
        (ir.nominal_check_receipt.as_bytes(), "nominal-check-receipt"),
        (ir.no_authority.as_bytes(), "no-authority"),
    ] {
        if is_zero(bytes) {
            return Err(
                DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1::MissingIdentity { field },
            );
        }
    }
    if ir.refinement_path != refinement_path.id() {
        return Err(DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1::PathIdentityMismatch);
    }
    for (matches, field) in [
        (
            refinement_path.source() == ir.source_geometry,
            "refinement-path-source",
        ),
        (
            refinement_path.target() == ir.target_geometry,
            "refinement-path-target",
        ),
    ] {
        if !matches {
            return Err(
                DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1::EndpointMismatch { field },
            );
        }
    }
    validate_fixed_resolution_candidate_path(ir, refinement_path, cx)?;

    let source_model = source
        .ir()
        .local_models
        .iter()
        .find(|model| model.id == ir.source_local_model)
        .ok_or(
            DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1::MissingLocalModel {
                field: "source-local-model",
            },
        )?;
    let target_model = target
        .ir()
        .local_models
        .iter()
        .find(|model| model.id == ir.target_local_model)
        .ok_or(
            DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1::MissingLocalModel {
                field: "target-local-model",
            },
        )?;
    let source_complex = source
        .ir()
        .complexes
        .iter()
        .find(|complex| complex.id == ir.source_complex)
        .ok_or(
            DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1::MissingCandidateComplex {
                field: "source-complex",
            },
        )?;
    let target_complex = target
        .ir()
        .complexes
        .iter()
        .find(|complex| complex.id == ir.target_complex)
        .ok_or(
            DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1::MissingCandidateComplex {
                field: "target-complex",
            },
        )?;

    for (matches, field) in [
        (
            source_complex.role == ir.complex_role,
            "source-complex-role",
        ),
        (
            target_complex.role == ir.complex_role,
            "target-complex-role",
        ),
        (
            local_model_complex_for_role(source_model, ir.complex_role) == ir.source_complex,
            "source-local-model-complex",
        ),
        (
            local_model_complex_for_role(target_model, ir.complex_role) == ir.target_complex,
            "target-local-model-complex",
        ),
    ] {
        if !matches {
            return Err(
                DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1::ComplexRoleMismatch {
                    field,
                },
            );
        }
    }
    for (matches, field) in [
        (
            source_model.chart == source_complex.chart,
            "source-model-complex-chart",
        ),
        (
            target_model.chart == target_complex.chart,
            "target-model-complex-chart",
        ),
        (
            source_model.chart == target_model.chart,
            "endpoint-local-model-chart",
        ),
    ] {
        if !matches {
            return Err(
                DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1::LocalPresentationMismatch {
                    field,
                },
            );
        }
    }
    for (matches, field) in [
        (
            source_complex.resolution.id == ir.source_resolution,
            "source-complex-resolution",
        ),
        (
            target_complex.resolution.id == ir.target_resolution,
            "target-complex-resolution",
        ),
    ] {
        if !matches {
            return Err(
                DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1::ResolutionScopeMismatch {
                    field,
                },
            );
        }
    }
    if source_model.locality != target_model.locality {
        return Err(DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1::LocalityMismatch);
    }
    let source_scope_witness = fixed_resolution_scope_witness(
        source_model,
        ir.source_resolution,
        "source-presentation",
        "source-presentation-resolution",
        "source-scope-witness",
    )?;
    let target_scope_witness = fixed_resolution_scope_witness(
        target_model,
        ir.target_resolution,
        "target-presentation",
        "target-presentation-resolution",
        "target-scope-witness",
    )?;

    Ok(FixedResolutionQuasiIsomorphismCandidateBindingV1 {
        refinement_path: ir.refinement_path,
        source_geometry: ir.source_geometry,
        target_geometry: ir.target_geometry,
        source_local_model: ir.source_local_model,
        target_local_model: ir.target_local_model,
        complex_role: ir.complex_role,
        source_complex: ir.source_complex,
        target_complex: ir.target_complex,
        source_resolution: ir.source_resolution,
        target_resolution: ir.target_resolution,
        source_scope_witness,
        target_scope_witness,
        nominal_theorem: ir.nominal_theorem,
        nominal_checker: ir.nominal_checker,
        nominal_check_receipt: ir.nominal_check_receipt,
        no_authority: ir.no_authority,
    })
}

/// Admit a structural fixed-resolution quasi-isomorphism candidate.
///
/// The supplied path must be an exact sealed homogeneous refinement path. Each
/// endpoint selector must resolve to a local model that owns the selected role
/// complex and explicitly declares its own matching `FixedResolution` scope.
/// The endpoint localities must be identical. The resulting token retains the
/// exact scope witnesses and nominal theorem/checker/receipt IDs, but does not
/// authenticate them or prove a chain map, commutation, cohomology isomorphism,
/// inverse, homotopy, presentation equivalence, refinement invariance, evidence
/// preservation, or physical equivalence. RD.1c owns independent promotion.
///
/// # Errors
/// Returns a typed refusal for schema, endpoint/path, local-model/complex/role,
/// fixed-resolution/locality, nominal-identity, cancellation, or canonical-
/// identity defects.
#[must_use = "a raw quasi-isomorphism candidate has no theorem authority"]
pub fn admit_derived_fixed_resolution_quasi_isomorphism_candidate_v1(
    ir: &DerivedFixedResolutionQuasiIsomorphismCandidateIrV1,
    source: &AdmittedDerivedGeometryV1,
    target: &AdmittedDerivedGeometryV1,
    refinement_path: &AdmittedDerivedMorphismV1,
    cx: &Cx<'_>,
) -> Result<
    AdmittedDerivedFixedResolutionQuasiIsomorphismCandidateV1,
    DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1,
> {
    if cx.checkpoint().is_err() {
        return Err(
            DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1::Cancelled {
                stage: "candidate-admission-entry",
            },
        );
    }
    if ir.schema_version != DERIVED_FIXED_RESOLUTION_QUASI_ISOMORPHISM_CANDIDATE_SCHEMA_VERSION_V1 {
        return Err(
            DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1::UnsupportedSchemaVersion {
                found: ir.schema_version,
                supported: DERIVED_FIXED_RESOLUTION_QUASI_ISOMORPHISM_CANDIDATE_SCHEMA_VERSION_V1,
            },
        );
    }
    let binding = fixed_resolution_quasi_isomorphism_candidate_binding(
        ir,
        source,
        target,
        refinement_path,
        cx,
    )?;
    let receipt = fixed_resolution_quasi_isomorphism_candidate_receipt(&binding, cx)?;
    if cx.checkpoint().is_err() {
        return Err(
            DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1::Cancelled {
                stage: "candidate-publication",
            },
        );
    }
    Ok(AdmittedDerivedFixedResolutionQuasiIsomorphismCandidateV1 {
        refinement_path: binding.refinement_path,
        source_geometry: binding.source_geometry,
        target_geometry: binding.target_geometry,
        source_local_model: binding.source_local_model,
        target_local_model: binding.target_local_model,
        complex_role: binding.complex_role,
        source_complex: binding.source_complex,
        target_complex: binding.target_complex,
        source_resolution: binding.source_resolution,
        target_resolution: binding.target_resolution,
        source_scope_witness: binding.source_scope_witness,
        target_scope_witness: binding.target_scope_witness,
        nominal_theorem: binding.nominal_theorem,
        nominal_checker: binding.nominal_checker,
        nominal_check_receipt: binding.nominal_check_receipt,
        no_authority: binding.no_authority,
        receipt,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use asupersync::types::Budget;
    use fs_exec::{CancelGate, ExecMode, StreamKey};

    use crate::derived::{
        CompactnessV1, ComplexDifferentialV1, ConfigurationChartClassV1, DerivedAdmissionBudgetV1,
        DerivedComplexRoleV1, DerivedGeometryIrV1, DerivedLinearMapIdV1, DerivedLocalModelClassV1,
        DerivedProofStateV1, DerivedQuantityKindIdV1, FiniteComputabilityV1, FiniteResolutionV1,
        GradedSpaceV1, LocalityScopeV1, RegularityClassV1, StratificationClassV1,
        StratificationIdV1, StratificationV1, StratumIdV1, StratumSpecV1, UnitBindingV1,
        admit_derived_geometry_v1,
    };

    fn with_cx<R>(cancelled: bool, f: impl FnOnce(&Cx<'_>) -> R) -> R {
        let gate = CancelGate::new_clock_free();
        if cancelled {
            gate.request();
        }
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 0x5244_3162,
                    kernel_id: 1,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            f(&cx)
        })
    }

    fn geometry_id(seed: u8) -> DerivedGeometryIdV1 {
        DerivedGeometryIdV1::parse_slice(&[seed; 32]).expect("nonzero geometry identity")
    }

    fn endpoint(seed: u8) -> GeometryEndpointV1<'static> {
        GeometryEndpointV1 {
            id: geometry_id(seed),
            subject: DerivedSubjectIdV1::from_bytes([1; 32]),
            model_version: DerivedModelVersionIdV1::from_bytes([4; 32]),
            category: GeometricCategoryV1::Semialgebraic,
            coefficients: CoefficientSystemV1::RationalReal,
            frame: DerivedFrameIdV1::from_bytes([2; 32]),
            unit_system: DerivedUnitSystemIdV1::from_bytes([3; 32]),
            charts: &[],
            complexes: &[],
        }
    }

    fn stratification_id(seed: u8) -> StratificationIdV1 {
        StratificationIdV1::from_bytes([seed; 32])
    }

    fn stratum_id(seed: u8) -> StratumIdV1 {
        StratumIdV1::from_bytes([seed; 32])
    }

    fn chart_id(seed: u8) -> ConfigurationChartIdV1 {
        ConfigurationChartIdV1::from_bytes([seed; 32])
    }

    fn chart(
        seed: u8,
        coordinate_dimension: u32,
        ambient_dimension: u32,
        quantity_seed: u8,
        scale_to_canonical: f64,
    ) -> ConfigurationChartV1 {
        let id = chart_id(seed);
        ConfigurationChartV1 {
            id,
            class: ConfigurationChartClassV1::Semialgebraic,
            coordinate_dimension,
            ambient_dimension,
            frame: DerivedFrameIdV1::from_bytes([2; 32]),
            coordinates: UnitBindingV1 {
                system: DerivedUnitSystemIdV1::from_bytes([3; 32]),
                quantity: DerivedQuantityKindIdV1::from_bytes([quantity_seed; 32]),
                scale_to_canonical,
            },
            locality: LocalityScopeV1::GermAt {
                chart: id,
                point: DerivedWitnessIdV1::from_bytes([seed.wrapping_add(1); 32]),
            },
            compactness: CompactnessV1::RelativelyCompact {
                witness: DerivedWitnessIdV1::from_bytes([seed.wrapping_add(2); 32]),
            },
            regularity: RegularityClassV1::Polynomial,
            computability: FiniteComputabilityV1::ExactFinite {
                kernel: DerivedWitnessIdV1::from_bytes([seed.wrapping_add(3); 32]),
            },
        }
    }

    fn endpoint_with_charts<'a>(
        seed: u8,
        charts: &'a [ConfigurationChartV1],
    ) -> GeometryEndpointV1<'a> {
        GeometryEndpointV1 {
            charts,
            ..endpoint(seed)
        }
    }

    fn complex_id(seed: u8) -> DerivedComplexIdV1 {
        DerivedComplexIdV1::from_bytes([seed; 32])
    }

    fn resolution_id(seed: u8) -> DerivedResolutionIdV1 {
        DerivedResolutionIdV1::from_bytes([seed; 32])
    }

    fn complex(
        seed: u8,
        resolution_seed: u8,
        chart: ConfigurationChartIdV1,
        role: DerivedComplexRoleV1,
        spaces: &[(i16, u32, u8)],
        truncation_order: u32,
    ) -> FiniteDerivedComplexV1 {
        let spaces = spaces
            .iter()
            .map(|(degree, dimension, quantity)| GradedSpaceV1 {
                degree: *degree,
                dimension: *dimension,
                quantity: DerivedQuantityKindIdV1::from_bytes([*quantity; 32]),
            })
            .collect::<Vec<_>>();
        let differentials = spaces
            .windows(2)
            .enumerate()
            .map(|(index, degrees)| ComplexDifferentialV1 {
                from_degree: degrees[0].degree,
                to_degree: degrees[1].degree,
                map: DerivedLinearMapIdV1::from_bytes([seed.wrapping_add(16 + index as u8); 32]),
                square_zero_witness: DerivedWitnessIdV1::from_bytes(
                    [seed.wrapping_add(32 + index as u8); 32],
                ),
            })
            .collect::<Vec<_>>();
        let min_degree = spaces.first().expect("nonempty complex fixture").degree;
        let max_degree = spaces.last().expect("nonempty complex fixture").degree;
        let max_basis_dimension = spaces
            .iter()
            .map(|space| space.dimension)
            .max()
            .unwrap_or(1)
            .max(1);
        FiniteDerivedComplexV1 {
            id: complex_id(seed),
            chart,
            role,
            spaces,
            differentials,
            resolution: FiniteResolutionV1 {
                id: resolution_id(resolution_seed),
                min_degree,
                max_degree,
                max_basis_dimension,
                truncation_order,
                remainder: (truncation_order > 0).then_some(DerivedWitnessIdV1::from_bytes(
                    [resolution_seed.wrapping_add(1); 32],
                )),
            },
            computability: FiniteComputabilityV1::ExactFinite {
                kernel: DerivedWitnessIdV1::from_bytes([seed.wrapping_add(48); 32]),
            },
        }
    }

    fn endpoint_with_parts<'a>(
        seed: u8,
        charts: &'a [ConfigurationChartV1],
        complexes: &'a [FiniteDerivedComplexV1],
    ) -> GeometryEndpointV1<'a> {
        GeometryEndpointV1 {
            charts,
            complexes,
            ..endpoint(seed)
        }
    }

    fn evidence_id(geometry: DerivedGeometryIdV1) -> DerivedEvidenceArtifactIdV1 {
        DerivedEvidenceArtifactIdV1::from_bytes(*geometry.as_bytes())
    }

    fn strict_ir(
        source: GeometryEndpointV1<'_>,
        target: GeometryEndpointV1<'_>,
        witness: u8,
        input_rank: ColorRank,
        output_rank: ColorRank,
    ) -> DerivedMorphismIrV1 {
        DerivedMorphismIrV1 {
            schema_version: DERIVED_MORPHISM_SCHEMA_VERSION_V1,
            source: source.id,
            target: target.id,
            kind: DerivedMorphismKindV1::Strict {
                witness: DerivedWitnessIdV1::from_bytes([witness; 32]),
            },
            evidence: DerivedEvidenceTransportV1::BalanceCorestrictionCovariant {
                input_geometry: source.id,
                output_geometry: target.id,
                input_evidence: evidence_id(source.id),
                output_evidence: evidence_id(target.id),
                input_rank,
                output_rank,
            },
            equivalence: DerivedEquivalenceBoundaryV1::NoClaim {
                artifact: DerivedNoClaimIdV1::from_bytes([witness.wrapping_add(64); 32]),
            },
        }
    }

    fn chart_map_ir(
        source: GeometryEndpointV1<'_>,
        target: GeometryEndpointV1<'_>,
        source_chart: ConfigurationChartIdV1,
        target_chart: ConfigurationChartIdV1,
        artifact_seed: u8,
        input_rank: ColorRank,
        output_rank: ColorRank,
    ) -> DerivedMorphismIrV1 {
        let mut ir = strict_ir(source, target, artifact_seed, input_rank, output_rank);
        ir.kind = DerivedMorphismKindV1::DeclaredChartMap {
            source_chart,
            target_chart,
            overlap: DerivedChartOverlapIdV1::from_bytes([artifact_seed.wrapping_add(1); 32]),
            map: DerivedChartMapIdV1::from_bytes([artifact_seed.wrapping_add(2); 32]),
        };
        ir
    }

    #[allow(clippy::too_many_arguments)]
    fn complex_refinement_ir(
        source: GeometryEndpointV1<'_>,
        target: GeometryEndpointV1<'_>,
        source_complex: DerivedComplexIdV1,
        target_complex: DerivedComplexIdV1,
        source_resolution: DerivedResolutionIdV1,
        target_resolution: DerivedResolutionIdV1,
        artifact_seed: u8,
        input_rank: ColorRank,
        output_rank: ColorRank,
    ) -> DerivedMorphismIrV1 {
        let mut ir = strict_ir(source, target, artifact_seed, input_rank, output_rank);
        ir.kind = DerivedMorphismKindV1::DeclaredComplexRefinement {
            source_complex,
            target_complex,
            source_resolution,
            target_resolution,
            prolongation: DerivedComplexRefinementMapIdV1::from_bytes(
                [artifact_seed.wrapping_add(1); 32],
            ),
            commutation: DerivedWitnessIdV1::from_bytes([artifact_seed.wrapping_add(2); 32]),
        };
        ir
    }

    fn inclusion_ir(
        source: GeometryEndpointV1<'_>,
        target: GeometryEndpointV1<'_>,
        artifact_seed: u8,
        input_rank: ColorRank,
        output_rank: ColorRank,
    ) -> DerivedMorphismIrV1 {
        let mut ir = strict_ir(source, target, artifact_seed, input_rank, output_rank);
        ir.kind = DerivedMorphismKindV1::DeclaredInclusion {
            map: DerivedInclusionMapIdV1::from_bytes([artifact_seed.wrapping_add(1); 32]),
            containment: DerivedWitnessIdV1::from_bytes([artifact_seed.wrapping_add(2); 32]),
        };
        ir
    }

    fn span_ir(
        source: GeometryEndpointV1<'_>,
        apex: GeometryEndpointV1<'_>,
        target: GeometryEndpointV1<'_>,
        left_leg: &AdmittedDerivedMorphismV1,
        right_leg: &AdmittedDerivedMorphismV1,
        no_claim_seed: u8,
    ) -> DerivedSpanCorrespondenceIrV1 {
        DerivedSpanCorrespondenceIrV1 {
            schema_version: DERIVED_SPAN_CORRESPONDENCE_SCHEMA_VERSION_V1,
            source: source.id,
            apex: apex.id,
            target: target.id,
            left_leg: left_leg.id(),
            right_leg: right_leg.id(),
            no_claim: DerivedNoClaimIdV1::from_bytes([no_claim_seed; 32]),
        }
    }

    fn restriction_ir(
        source: GeometryEndpointV1<'_>,
        target: GeometryEndpointV1<'_>,
        witness: u8,
        input_rank: ColorRank,
        output_rank: ColorRank,
    ) -> DerivedMorphismIrV1 {
        let mut ir = strict_ir(source, target, witness, input_rank, output_rank);
        ir.evidence = DerivedEvidenceTransportV1::RestrictionContravariant {
            input_geometry: target.id,
            output_geometry: source.id,
            input_evidence: evidence_id(target.id),
            output_evidence: evidence_id(source.id),
            input_rank,
            output_rank,
        };
        ir
    }

    fn admit_strict(
        source: GeometryEndpointV1<'_>,
        target: GeometryEndpointV1<'_>,
        witness: u8,
        input_rank: ColorRank,
        output_rank: ColorRank,
        cx: &Cx<'_>,
    ) -> AdmittedDerivedMorphismV1 {
        admit_between_endpoints(
            strict_ir(source, target, witness, input_rank, output_rank),
            source,
            target,
            cx,
        )
        .expect("valid strict morphism")
    }

    #[allow(clippy::too_many_arguments)]
    fn admit_chart_map(
        source: GeometryEndpointV1<'_>,
        target: GeometryEndpointV1<'_>,
        source_chart: ConfigurationChartIdV1,
        target_chart: ConfigurationChartIdV1,
        artifact_seed: u8,
        input_rank: ColorRank,
        output_rank: ColorRank,
        cx: &Cx<'_>,
    ) -> AdmittedDerivedMorphismV1 {
        admit_between_endpoints(
            chart_map_ir(
                source,
                target,
                source_chart,
                target_chart,
                artifact_seed,
                input_rank,
                output_rank,
            ),
            source,
            target,
            cx,
        )
        .expect("valid declared chart map")
    }

    #[allow(clippy::too_many_arguments)]
    fn admit_complex_refinement(
        source: GeometryEndpointV1<'_>,
        target: GeometryEndpointV1<'_>,
        source_complex: DerivedComplexIdV1,
        target_complex: DerivedComplexIdV1,
        source_resolution: DerivedResolutionIdV1,
        target_resolution: DerivedResolutionIdV1,
        artifact_seed: u8,
        input_rank: ColorRank,
        output_rank: ColorRank,
        cx: &Cx<'_>,
    ) -> AdmittedDerivedMorphismV1 {
        admit_between_endpoints(
            complex_refinement_ir(
                source,
                target,
                source_complex,
                target_complex,
                source_resolution,
                target_resolution,
                artifact_seed,
                input_rank,
                output_rank,
            ),
            source,
            target,
            cx,
        )
        .expect("valid declared finite-complex refinement")
    }

    fn admit_inclusion(
        source: GeometryEndpointV1<'_>,
        target: GeometryEndpointV1<'_>,
        artifact_seed: u8,
        input_rank: ColorRank,
        output_rank: ColorRank,
        cx: &Cx<'_>,
    ) -> AdmittedDerivedMorphismV1 {
        admit_between_endpoints(
            inclusion_ir(source, target, artifact_seed, input_rank, output_rank),
            source,
            target,
            cx,
        )
        .expect("valid declared inclusion")
    }

    fn try_complex_refinement(
        source: GeometryEndpointV1<'_>,
        target: GeometryEndpointV1<'_>,
        source_complex: &FiniteDerivedComplexV1,
        target_complex: &FiniteDerivedComplexV1,
        cx: &Cx<'_>,
    ) -> Result<AdmittedDerivedMorphismV1, DerivedMorphismErrorV1> {
        admit_between_endpoints(
            complex_refinement_ir(
                source,
                target,
                source_complex.id,
                target_complex.id,
                source_complex.resolution.id,
                target_complex.resolution.id,
                30,
                ColorRank::Validated,
                ColorRank::Estimated,
            ),
            source,
            target,
            cx,
        )
    }

    fn admit_identity(object: GeometryEndpointV1<'_>, cx: &Cx<'_>) -> AdmittedDerivedMorphismV1 {
        admit_between_endpoints(
            DerivedMorphismIrV1 {
                schema_version: DERIVED_MORPHISM_SCHEMA_VERSION_V1,
                source: object.id,
                target: object.id,
                kind: DerivedMorphismKindV1::Identity,
                evidence: DerivedEvidenceTransportV1::Identity,
                equivalence: DerivedEquivalenceBoundaryV1::IdentityOnly,
            },
            object,
            object,
            cx,
        )
        .expect("valid identity morphism")
    }

    fn fixed_resolution_geometry_ir(
        tangent_complex_seed: u8,
        tangent_resolution_seed: u8,
        local_model_seed: u8,
        tangent_rank: u32,
    ) -> DerivedGeometryIrV1 {
        let chart = chart(4, 2, 2, 8, 1.0);
        let chart_id = chart.id;
        let tangent = complex(
            tangent_complex_seed,
            tangent_resolution_seed,
            chart_id,
            DerivedComplexRoleV1::Tangent,
            &[(0, tangent_rank, 120), (1, 1, 121)],
            0,
        );
        let cotangent = complex(
            tangent_complex_seed.wrapping_add(1),
            tangent_resolution_seed.wrapping_add(1),
            chart_id,
            DerivedComplexRoleV1::Cotangent,
            &[(0, 1, 122), (1, 1, 123)],
            0,
        );
        let deformation = complex(
            tangent_complex_seed.wrapping_add(2),
            tangent_resolution_seed.wrapping_add(2),
            chart_id,
            DerivedComplexRoleV1::DeformationObstruction,
            &[(0, 1, 124), (1, 1, 125)],
            0,
        );
        let local_model = DerivedLocalModelIdV1::from_bytes([local_model_seed; 32]);
        let stratum = StratumIdV1::from_bytes([local_model_seed.wrapping_add(1); 32]);
        let locality = LocalityScopeV1::GermAt {
            chart: chart_id,
            point: DerivedWitnessIdV1::from_bytes([5; 32]),
        };
        DerivedGeometryIrV1 {
            schema_version: crate::derived::DERIVED_GEOMETRY_SCHEMA_VERSION_V1,
            subject: DerivedSubjectIdV1::from_bytes([1; 32]),
            model_version: DerivedModelVersionIdV1::from_bytes([4; 32]),
            category: GeometricCategoryV1::Semialgebraic,
            coefficients: CoefficientSystemV1::RationalReal,
            frame: DerivedFrameIdV1::from_bytes([2; 32]),
            unit_system: DerivedUnitSystemIdV1::from_bytes([3; 32]),
            locality,
            compactness: CompactnessV1::RelativelyCompact {
                witness: DerivedWitnessIdV1::from_bytes([6; 32]),
            },
            charts: vec![chart],
            equalities: Vec::new(),
            inequalities: Vec::new(),
            boundaries: Vec::new(),
            contacts: Vec::new(),
            constitutive_data: Vec::new(),
            complexes: vec![tangent, cotangent, deformation],
            local_models: vec![DerivedLocalModelV1 {
                id: local_model,
                chart: chart_id,
                class: DerivedLocalModelClassV1::GeneralFiniteDerived,
                equalities: Vec::new(),
                active_inequalities: Vec::new(),
                active_contacts: Vec::new(),
                constitutive_data: Vec::new(),
                tangent_complex: complex_id(tangent_complex_seed),
                cotangent_complex: complex_id(tangent_complex_seed.wrapping_add(1)),
                deformation_complex: complex_id(tangent_complex_seed.wrapping_add(2)),
                virtual_dimension: 1,
                locality,
                presentation: PresentationScopeV1::FixedResolution {
                    resolution: resolution_id(tangent_resolution_seed),
                    witness: DerivedWitnessIdV1::from_bytes(
                        [tangent_resolution_seed.wrapping_add(32); 32],
                    ),
                },
            }],
            stratification: StratificationV1 {
                id: StratificationIdV1::from_bytes([local_model_seed.wrapping_add(2); 32]),
                class: StratificationClassV1::FiniteIncidence,
                strata: vec![StratumSpecV1 {
                    id: stratum,
                    chart: chart_id,
                    local_model,
                    dimension: 1,
                    active_inequalities: Vec::new(),
                    active_contacts: Vec::new(),
                    relative_boundary: None,
                    regularity: RegularityClassV1::Polynomial,
                    compactness: CompactnessV1::RelativelyCompact {
                        witness: DerivedWitnessIdV1::from_bytes(
                            [local_model_seed.wrapping_add(3); 32],
                        ),
                    },
                }],
                incidences: Vec::new(),
                local_links: Vec::new(),
            },
            proof_state: DerivedProofStateV1::StructuralNoClaim {
                no_claim: DerivedNoClaimIdV1::from_bytes([local_model_seed.wrapping_add(4); 32]),
            },
        }
    }

    fn fixed_resolution_candidate_fixture(
        cx: &Cx<'_>,
    ) -> (
        AdmittedDerivedGeometryV1,
        AdmittedDerivedGeometryV1,
        AdmittedDerivedMorphismV1,
        DerivedFixedResolutionQuasiIsomorphismCandidateIrV1,
    ) {
        let source = admit_derived_geometry_v1(
            fixed_resolution_geometry_ir(70, 80, 90, 1),
            DerivedAdmissionBudgetV1::STANDARD,
            cx,
        )
        .expect("valid source fixed-resolution geometry");
        let target = admit_derived_geometry_v1(
            fixed_resolution_geometry_ir(73, 83, 93, 2),
            DerivedAdmissionBudgetV1::STANDARD,
            cx,
        )
        .expect("valid target fixed-resolution geometry");
        let path = admit_derived_morphism_v1(
            DerivedMorphismIrV1 {
                schema_version: DERIVED_MORPHISM_SCHEMA_VERSION_V1,
                source: source.id(),
                target: target.id(),
                kind: DerivedMorphismKindV1::DeclaredComplexRefinement {
                    source_complex: complex_id(70),
                    target_complex: complex_id(73),
                    source_resolution: resolution_id(80),
                    target_resolution: resolution_id(83),
                    prolongation: DerivedComplexRefinementMapIdV1::from_bytes([101; 32]),
                    commutation: DerivedWitnessIdV1::from_bytes([102; 32]),
                },
                evidence: DerivedEvidenceTransportV1::BalanceCorestrictionCovariant {
                    input_geometry: source.id(),
                    output_geometry: target.id(),
                    input_evidence: evidence_id(source.id()),
                    output_evidence: evidence_id(target.id()),
                    input_rank: ColorRank::Validated,
                    output_rank: ColorRank::Estimated,
                },
                equivalence: DerivedEquivalenceBoundaryV1::NoClaim {
                    artifact: DerivedNoClaimIdV1::from_bytes([103; 32]),
                },
            },
            &source,
            &target,
            cx,
        )
        .expect("valid fixed-resolution refinement path");
        let ir = DerivedFixedResolutionQuasiIsomorphismCandidateIrV1 {
            schema_version: DERIVED_FIXED_RESOLUTION_QUASI_ISOMORPHISM_CANDIDATE_SCHEMA_VERSION_V1,
            source_geometry: source.id(),
            target_geometry: target.id(),
            source_local_model: DerivedLocalModelIdV1::from_bytes([90; 32]),
            target_local_model: DerivedLocalModelIdV1::from_bytes([93; 32]),
            complex_role: DerivedComplexRoleV1::Tangent,
            source_complex: complex_id(70),
            target_complex: complex_id(73),
            source_resolution: resolution_id(80),
            target_resolution: resolution_id(83),
            refinement_path: path.id(),
            nominal_theorem: DerivedTheoremIdV1::from_bytes([104; 32]),
            nominal_checker: DerivedCheckerIdV1::from_bytes([105; 32]),
            nominal_check_receipt: DerivedWitnessIdV1::from_bytes([106; 32]),
            no_authority: DerivedNoClaimIdV1::from_bytes([107; 32]),
        };
        (source, target, path, ir)
    }

    fn admitted_fixed_resolution_geometry(
        tangent_complex_seed: u8,
        tangent_resolution_seed: u8,
        local_model_seed: u8,
        tangent_rank: u32,
        cx: &Cx<'_>,
    ) -> AdmittedDerivedGeometryV1 {
        admit_derived_geometry_v1(
            fixed_resolution_geometry_ir(
                tangent_complex_seed,
                tangent_resolution_seed,
                local_model_seed,
                tangent_rank,
            ),
            DerivedAdmissionBudgetV1::STANDARD,
            cx,
        )
        .expect("valid admitted stratum fixture")
    }

    fn sole_stratum_object(geometry: &AdmittedDerivedGeometryV1) -> DerivedStratumObjectV1 {
        DerivedStratumObjectV1 {
            geometry: geometry.id(),
            stratification: geometry.ir().stratification.id,
            stratum: geometry.ir().stratification.strata[0].id,
        }
    }

    fn stratum_component_ir(
        source: DerivedStratumObjectV1,
        target: DerivedStratumObjectV1,
        seed: u8,
    ) -> DerivedStratumMorphismIrV1 {
        DerivedStratumMorphismIrV1 {
            schema_version: DERIVED_STRATUM_MORPHISM_SCHEMA_VERSION_V1,
            source,
            target,
            kind: DerivedStratumMorphismKindV1::DeclaredComponent {
                map: DerivedStratumMapIdV1::from_bytes([seed; 32]),
                constructibility: DerivedConstructibilityDeclarationIdV1::from_bytes(
                    [seed.wrapping_add(1); 32],
                ),
            },
            authority: DerivedStratumAuthorityBoundaryV1::NoClaim {
                artifact: DerivedNoClaimIdV1::from_bytes([seed.wrapping_add(2); 32]),
            },
        }
    }

    fn admit_stratum_component(
        source: &AdmittedDerivedGeometryV1,
        target: &AdmittedDerivedGeometryV1,
        seed: u8,
        cx: &Cx<'_>,
    ) -> AdmittedDerivedStratumMorphismV1 {
        admit_stratum_component_between(
            source,
            target,
            sole_stratum_object(source),
            sole_stratum_object(target),
            seed,
            cx,
        )
    }

    fn admit_stratum_component_between(
        source: &AdmittedDerivedGeometryV1,
        target: &AdmittedDerivedGeometryV1,
        source_object: DerivedStratumObjectV1,
        target_object: DerivedStratumObjectV1,
        seed: u8,
        cx: &Cx<'_>,
    ) -> AdmittedDerivedStratumMorphismV1 {
        admit_derived_stratum_morphism_v1(
            &stratum_component_ir(source_object, target_object, seed),
            source,
            target,
            cx,
        )
        .expect("valid declared stratum component")
    }

    fn exhaustive_stratified_map_candidate_ir(
        source: &AdmittedDerivedGeometryV1,
        target: &AdmittedDerivedGeometryV1,
        components: &[AdmittedDerivedStratumMorphismV1],
        seed: u8,
    ) -> DerivedExhaustiveStratifiedMapCandidateIrV1 {
        DerivedExhaustiveStratifiedMapCandidateIrV1 {
            schema_version: DERIVED_EXHAUSTIVE_STRATIFIED_MAP_CANDIDATE_SCHEMA_VERSION_V1,
            source_geometry: source.id(),
            source_stratification: source.ir().stratification.id,
            target_geometry: target.id(),
            target_stratification: target.ir().stratification.id,
            components: components
                .iter()
                .map(|component| DerivedStratifiedMapComponentBindingV1 {
                    source_stratum: component.source().stratum,
                    target_stratum: component.target().stratum,
                    component: component.id(),
                })
                .collect(),
            nominal_assembly: DerivedStratifiedMapAssemblyIdV1::from_bytes([seed; 32]),
            nominal_constructibility: DerivedGlobalConstructibilityDeclarationIdV1::from_bytes(
                [seed.wrapping_add(1); 32],
            ),
            no_authority: DerivedNoClaimIdV1::from_bytes([seed.wrapping_add(2); 32]),
        }
    }

    #[test]
    fn v1_class_bytes_keep_old_tags_and_domain_separate_new_paths() {
        assert_eq!(
            <DerivedMorphismIdentitySchemaV1 as CanonicalSchema>::CONTEXT,
            "typed endpoints, strict map class, evidence variance, no-equivalence boundary, and ordered primitive lineage"
        );
        for (class, expected) in [
            (ReceiptClassV1::Identity, 0),
            (ReceiptClassV1::CompositeStrict, 2),
            (ReceiptClassV1::CompositeDeclaredChartMap, 4),
            (ReceiptClassV1::CompositeHeterogeneous, 5),
            (ReceiptClassV1::CompositeDeclaredComplexRefinement, 7),
            (ReceiptClassV1::CompositeDeclaredInclusion, 9),
        ] {
            assert_eq!(class_bytes(class).as_slice(), &[expected]);
        }
        let strict = class_bytes(ReceiptClassV1::PrimitiveStrict(
            DerivedWitnessIdV1::from_bytes([7; 32]),
        ));
        assert_eq!(strict.as_slice().len(), 33);
        assert_eq!(strict.as_slice()[0], 1);
        assert!(strict.as_slice()[1..].iter().all(|byte| *byte == 7));

        let chart = class_bytes(ReceiptClassV1::PrimitiveDeclaredChartMap {
            source_chart: chart_id(8),
            target_chart: chart_id(9),
            overlap: DerivedChartOverlapIdV1::from_bytes([10; 32]),
            map: DerivedChartMapIdV1::from_bytes([11; 32]),
        });
        assert_eq!(chart.as_slice().len(), 129);
        assert_eq!(chart.as_slice()[0], 3);
        assert_eq!(&chart.as_slice()[1..33], chart_id(8).as_bytes());
        assert_eq!(&chart.as_slice()[33..65], chart_id(9).as_bytes());
        assert_eq!(
            &chart.as_slice()[65..97],
            DerivedChartOverlapIdV1::from_bytes([10; 32]).as_bytes()
        );
        assert_eq!(
            &chart.as_slice()[97..129],
            DerivedChartMapIdV1::from_bytes([11; 32]).as_bytes()
        );

        let refinement = class_bytes(ReceiptClassV1::PrimitiveDeclaredComplexRefinement {
            source_complex: complex_id(12),
            target_complex: complex_id(13),
            source_resolution: resolution_id(14),
            target_resolution: resolution_id(15),
            prolongation: DerivedComplexRefinementMapIdV1::from_bytes([16; 32]),
            commutation: DerivedWitnessIdV1::from_bytes([17; 32]),
        });
        assert_eq!(refinement.as_slice().len(), 193);
        assert_eq!(refinement.as_slice()[0], 6);
        for (range, expected) in [
            (1..33, 12),
            (33..65, 13),
            (65..97, 14),
            (97..129, 15),
            (129..161, 16),
            (161..193, 17),
        ] {
            assert!(
                refinement.as_slice()[range]
                    .iter()
                    .all(|byte| *byte == expected)
            );
        }

        let inclusion = class_bytes(ReceiptClassV1::PrimitiveDeclaredInclusion {
            map: DerivedInclusionMapIdV1::from_bytes([18; 32]),
            containment: DerivedWitnessIdV1::from_bytes([19; 32]),
        });
        assert_eq!(inclusion.as_slice().len(), 65);
        assert_eq!(inclusion.as_slice()[0], 8);
        assert!(inclusion.as_slice()[1..33].iter().all(|byte| *byte == 18));
        assert!(inclusion.as_slice()[33..65].iter().all(|byte| *byte == 19));
    }

    #[test]
    fn stratum_morphisms_have_a_separate_domain_and_frozen_class_tags() {
        assert_ne!(
            <DerivedStratumMorphismIdentitySchemaV1 as CanonicalSchema>::DOMAIN,
            <DerivedMorphismIdentitySchemaV1 as CanonicalSchema>::DOMAIN
        );
        assert_ne!(
            <DerivedStratumMorphismIdentitySchemaV1 as CanonicalSchema>::DOMAIN,
            <DerivedSpanCorrespondenceIdentitySchemaV1 as CanonicalSchema>::DOMAIN
        );
        assert_eq!(
            <DerivedStratumMorphismIdentitySchemaV1 as CanonicalSchema>::FIELDS.len(),
            9
        );
        assert_eq!(DERIVED_STRATUM_MORPHISM_IDENTITY_LIMITS_V1.max_fields(), 9);
        assert_eq!(
            <DerivedFixedResolutionQuasiIsomorphismCandidateIdentitySchemaV1 as CanonicalSchema>::FIELDS
                .len(),
            16
        );
        assert_eq!(
            DERIVED_FIXED_RESOLUTION_QUASI_ISOMORPHISM_CANDIDATE_IDENTITY_LIMITS_V1.max_fields(),
            16
        );
        assert_eq!(
            stratum_class_bytes(StratumReceiptClassV1::Identity).as_slice(),
            &[0]
        );
        assert_eq!(
            stratum_class_bytes(StratumReceiptClassV1::CompositeDeclaredPath).as_slice(),
            &[2]
        );
        let primitive = stratum_class_bytes(StratumReceiptClassV1::PrimitiveDeclaredComponent {
            map: DerivedStratumMapIdV1::from_bytes([31; 32]),
            constructibility: DerivedConstructibilityDeclarationIdV1::from_bytes([32; 32]),
        });
        assert_eq!(primitive.as_slice().len(), 65);
        assert_eq!(primitive.as_slice()[0], 1);
        assert!(primitive.as_slice()[1..33].iter().all(|byte| *byte == 31));
        assert!(primitive.as_slice()[33..65].iter().all(|byte| *byte == 32));
    }

    #[test]
    fn stratum_component_admission_binds_exact_selectors_and_replays() {
        with_cx(false, |cx| {
            let source = admitted_fixed_resolution_geometry(70, 80, 90, 1, cx);
            let target = admitted_fixed_resolution_geometry(73, 83, 93, 2, cx);
            let source_object = sole_stratum_object(&source);
            let target_object = sole_stratum_object(&target);
            let ir = stratum_component_ir(source_object, target_object, 110);
            let first = admit_derived_stratum_morphism_v1(&ir, &source, &target, cx)
                .expect("valid stratum component");
            let replay = admit_derived_stratum_morphism_v1(&ir, &source, &target, cx)
                .expect("deterministic replay");

            assert_eq!(first, replay);
            assert_eq!(first.source(), source_object);
            assert_eq!(first.target(), target_object);
            assert_eq!(
                first.class(),
                AdmittedDerivedStratumMorphismClassV1::DeclaredPath
            );
            assert_eq!(first.primitive_factors(), &[first.id()]);
            assert_eq!(
                first.no_authority_claims(),
                &[DerivedNoClaimIdV1::from_bytes([112; 32])]
            );
            assert_eq!(
                first.primitive_path(),
                &[DeclaredStratumMapPrimitiveV1 {
                    source: source_object,
                    target: target_object,
                    map: DerivedStratumMapIdV1::from_bytes([110; 32]),
                    constructibility: DerivedConstructibilityDeclarationIdV1::from_bytes([111; 32]),
                }]
            );

            let changed = admit_stratum_component(&source, &target, 113, cx);
            assert_ne!(first.id(), changed.id());
        });
    }

    #[test]
    fn stratum_component_admission_rejects_unowned_selectors_and_authority() {
        with_cx(false, |cx| {
            let source = admitted_fixed_resolution_geometry(70, 80, 90, 1, cx);
            let target = admitted_fixed_resolution_geometry(73, 83, 93, 2, cx);
            let mut ir = stratum_component_ir(
                sole_stratum_object(&source),
                sole_stratum_object(&target),
                120,
            );

            ir.source.stratification = stratification_id(121);
            assert!(matches!(
                admit_derived_stratum_morphism_v1(&ir, &source, &target, cx),
                Err(DerivedStratumMorphismErrorV1::MissingStratification {
                    field: "source-stratification"
                })
            ));

            ir.source = sole_stratum_object(&source);
            ir.target.stratum = stratum_id(122);
            assert!(matches!(
                admit_derived_stratum_morphism_v1(&ir, &source, &target, cx),
                Err(DerivedStratumMorphismErrorV1::MissingStratum {
                    field: "target-stratum"
                })
            ));

            ir.target = sole_stratum_object(&target);
            ir.source.stratum = StratumIdV1::from_bytes([0; 32]);
            assert!(matches!(
                admit_derived_stratum_morphism_v1(&ir, &source, &target, cx),
                Err(DerivedStratumMorphismErrorV1::MissingIdentity {
                    field: "source-stratum"
                })
            ));

            ir.source = sole_stratum_object(&source);
            ir.kind = DerivedStratumMorphismKindV1::DeclaredComponent {
                map: DerivedStratumMapIdV1::from_bytes([0; 32]),
                constructibility: DerivedConstructibilityDeclarationIdV1::from_bytes([123; 32]),
            };
            assert!(matches!(
                admit_derived_stratum_morphism_v1(&ir, &source, &target, cx),
                Err(DerivedStratumMorphismErrorV1::MissingIdentity {
                    field: "stratum-map"
                })
            ));

            ir.kind = DerivedStratumMorphismKindV1::DeclaredComponent {
                map: DerivedStratumMapIdV1::from_bytes([124; 32]),
                constructibility: DerivedConstructibilityDeclarationIdV1::from_bytes([0; 32]),
            };
            assert!(matches!(
                admit_derived_stratum_morphism_v1(&ir, &source, &target, cx),
                Err(DerivedStratumMorphismErrorV1::MissingIdentity {
                    field: "constructibility-declaration"
                })
            ));

            ir.kind = DerivedStratumMorphismKindV1::DeclaredComponent {
                map: DerivedStratumMapIdV1::from_bytes([124; 32]),
                constructibility: DerivedConstructibilityDeclarationIdV1::from_bytes([125; 32]),
            };
            ir.authority = DerivedStratumAuthorityBoundaryV1::NoClaim {
                artifact: DerivedNoClaimIdV1::from_bytes([0; 32]),
            };
            assert!(matches!(
                admit_derived_stratum_morphism_v1(&ir, &source, &target, cx),
                Err(DerivedStratumMorphismErrorV1::MissingIdentity {
                    field: "no-stratum-map-authority"
                })
            ));

            ir.authority = DerivedStratumAuthorityBoundaryV1::IdentityOnly;
            assert_eq!(
                admit_derived_stratum_morphism_v1(&ir, &source, &target, cx),
                Err(DerivedStratumMorphismErrorV1::AuthorityLaundering)
            );
        });
    }

    #[test]
    fn stratum_composition_requires_the_exact_stratum_within_one_geometry() {
        with_cx(false, |cx| {
            let x = admitted_fixed_resolution_geometry(70, 80, 90, 1, cx);
            let mut middle_ir = fixed_resolution_geometry_ir(73, 83, 93, 2);
            let first_middle_stratum = middle_ir.stratification.strata[0].id;
            let second_middle_stratum = stratum_id(150);
            let mut second = middle_ir.stratification.strata[0].clone();
            second.id = second_middle_stratum;
            middle_ir.stratification.strata.push(second);
            let middle =
                admit_derived_geometry_v1(middle_ir, DerivedAdmissionBudgetV1::STANDARD, cx)
                    .expect("valid geometry with two strata");
            let z = admitted_fixed_resolution_geometry(76, 86, 96, 3, cx);
            let middle_object = |stratum| DerivedStratumObjectV1 {
                geometry: middle.id(),
                stratification: middle.ir().stratification.id,
                stratum,
            };

            let f = admit_derived_stratum_morphism_v1(
                &stratum_component_ir(
                    sole_stratum_object(&x),
                    middle_object(first_middle_stratum),
                    151,
                ),
                &x,
                &middle,
                cx,
            )
            .expect("component into first middle stratum");
            let exact = admit_derived_stratum_morphism_v1(
                &stratum_component_ir(
                    middle_object(first_middle_stratum),
                    sole_stratum_object(&z),
                    154,
                ),
                &middle,
                &z,
                cx,
            )
            .expect("component from exact middle stratum");
            let wrong = admit_derived_stratum_morphism_v1(
                &stratum_component_ir(
                    middle_object(second_middle_stratum),
                    sole_stratum_object(&z),
                    154,
                ),
                &middle,
                &z,
                cx,
            )
            .expect("component from other middle stratum");

            assert!(compose_derived_stratum_morphisms_v1(&f, &exact, cx).is_ok());
            assert_eq!(
                compose_derived_stratum_morphisms_v1(&f, &wrong, cx),
                Err(DerivedStratumMorphismErrorV1::CompositionEndpointMismatch)
            );
            assert_ne!(exact.id(), wrong.id());
        });
    }

    #[test]
    fn stratum_components_allow_different_finite_dimensions() {
        with_cx(false, |cx| {
            let source = admitted_fixed_resolution_geometry(70, 80, 90, 1, cx);
            let mut target_ir = fixed_resolution_geometry_ir(73, 83, 93, 2);
            target_ir.stratification.strata[0].dimension = 2;
            let target =
                admit_derived_geometry_v1(target_ir, DerivedAdmissionBudgetV1::STANDARD, cx)
                    .expect("valid dimension-two target stratum");

            assert_eq!(source.ir().stratification.strata[0].dimension, 1);
            assert_eq!(target.ir().stratification.strata[0].dimension, 2);
            assert!(
                admit_derived_stratum_morphism_v1(
                    &stratum_component_ir(
                        sole_stratum_object(&source),
                        sole_stratum_object(&target),
                        160,
                    ),
                    &source,
                    &target,
                    cx,
                )
                .is_ok()
            );

            let mut incompatible_ir = fixed_resolution_geometry_ir(76, 86, 96, 3);
            incompatible_ir.model_version = DerivedModelVersionIdV1::from_bytes([5; 32]);
            let incompatible =
                admit_derived_geometry_v1(incompatible_ir, DerivedAdmissionBudgetV1::STANDARD, cx)
                    .expect("valid independently versioned target");
            assert_eq!(
                admit_derived_stratum_morphism_v1(
                    &stratum_component_ir(
                        sole_stratum_object(&source),
                        sole_stratum_object(&incompatible),
                        163,
                    ),
                    &source,
                    &incompatible,
                    cx,
                ),
                Err(DerivedStratumMorphismErrorV1::ModelVersionMismatch)
            );
        });
    }

    #[test]
    fn stratum_identity_is_neutral_and_composition_is_associative() {
        with_cx(false, |cx| {
            let x = admitted_fixed_resolution_geometry(70, 80, 90, 1, cx);
            let y = admitted_fixed_resolution_geometry(73, 83, 93, 2, cx);
            let z = admitted_fixed_resolution_geometry(76, 86, 96, 3, cx);
            let w = admitted_fixed_resolution_geometry(79, 89, 99, 4, cx);
            let f = admit_stratum_component(&x, &y, 130, cx);
            let g = admit_stratum_component(&y, &z, 133, cx);
            let h = admit_stratum_component(&z, &w, 136, cx);
            let identity_x =
                identity_derived_stratum_morphism_v1(&x, sole_stratum_object(&x).stratum, cx)
                    .expect("source identity");
            let identity_w =
                identity_derived_stratum_morphism_v1(&w, sole_stratum_object(&w).stratum, cx)
                    .expect("target identity");

            assert_eq!(
                compose_derived_stratum_morphisms_v1(&identity_x, &f, cx).expect("left identity"),
                f
            );
            assert_eq!(
                compose_derived_stratum_morphisms_v1(&h, &identity_w, cx).expect("right identity"),
                h
            );

            let fg = compose_derived_stratum_morphisms_v1(&f, &g, cx).expect("f then g");
            let gh = compose_derived_stratum_morphisms_v1(&g, &h, cx).expect("g then h");
            let left = compose_derived_stratum_morphisms_v1(&fg, &h, cx).expect("(fg) then h");
            let right = compose_derived_stratum_morphisms_v1(&f, &gh, cx).expect("f then (gh)");
            assert_eq!(left, right);
            assert_eq!(left.primitive_path().len(), 3);
            assert_eq!(left.primitive_factors(), &[f.id(), g.id(), h.id()]);
            assert_eq!(left.no_authority_claims().len(), 3);

            assert_eq!(
                compose_derived_stratum_morphisms_v1(&f, &h, cx),
                Err(DerivedStratumMorphismErrorV1::CompositionEndpointMismatch)
            );

            let a = admit_stratum_component(&x, &x, 170, cx);
            let b = admit_stratum_component(&x, &x, 173, cx);
            let ab = compose_derived_stratum_morphisms_v1(&a, &b, cx).expect("a then b");
            let ba = compose_derived_stratum_morphisms_v1(&b, &a, cx).expect("b then a");
            assert_ne!(ab.id(), ba.id());
            assert_eq!(ab.primitive_factors(), &[a.id(), b.id()]);
        });
    }

    #[test]
    fn stratum_identity_rejects_changed_objects_and_cancelled_work_publishes_nothing() {
        let (source, target, admitted) = with_cx(false, |cx| {
            let source = admitted_fixed_resolution_geometry(70, 80, 90, 1, cx);
            let target = admitted_fixed_resolution_geometry(73, 83, 93, 2, cx);
            let admitted = admit_stratum_component(&source, &target, 140, cx);
            (source, target, admitted)
        });
        with_cx(false, |cx| {
            let object = sole_stratum_object(&source);
            let wrong_boundary = DerivedStratumMorphismIrV1 {
                schema_version: DERIVED_STRATUM_MORPHISM_SCHEMA_VERSION_V1,
                source: object,
                target: object,
                kind: DerivedStratumMorphismKindV1::Identity,
                authority: DerivedStratumAuthorityBoundaryV1::NoClaim {
                    artifact: DerivedNoClaimIdV1::from_bytes([141; 32]),
                },
            };
            assert_eq!(
                admit_derived_stratum_morphism_v1(&wrong_boundary, &source, &source, cx),
                Err(DerivedStratumMorphismErrorV1::InvalidIdentity)
            );

            let ir = DerivedStratumMorphismIrV1 {
                schema_version: DERIVED_STRATUM_MORPHISM_SCHEMA_VERSION_V1,
                source: sole_stratum_object(&source),
                target: sole_stratum_object(&target),
                kind: DerivedStratumMorphismKindV1::Identity,
                authority: DerivedStratumAuthorityBoundaryV1::IdentityOnly,
            };
            assert_eq!(
                admit_derived_stratum_morphism_v1(&ir, &source, &target, cx),
                Err(DerivedStratumMorphismErrorV1::InvalidIdentity)
            );
        });
        with_cx(true, |cx| {
            let ir = stratum_component_ir(
                sole_stratum_object(&source),
                sole_stratum_object(&target),
                143,
            );
            assert!(matches!(
                admit_derived_stratum_morphism_v1(&ir, &source, &target, cx),
                Err(DerivedStratumMorphismErrorV1::Cancelled {
                    stage: "stratum-admission-entry"
                })
            ));
            assert!(matches!(
                compose_derived_stratum_morphisms_v1(&admitted, &admitted, cx),
                Err(DerivedStratumMorphismErrorV1::Cancelled {
                    stage: "stratum-composition-entry"
                })
            ));
        });
    }

    #[test]
    #[allow(clippy::too_many_lines)] // Replay plus independent retained-field identity movement.
    fn exhaustive_stratified_map_candidate_is_domain_separate_and_replays() {
        with_cx(false, |cx| {
            assert_ne!(
                <DerivedExhaustiveStratifiedMapCandidateIdentitySchemaV1 as CanonicalSchema>::DOMAIN,
                <DerivedStratumMorphismIdentitySchemaV1 as CanonicalSchema>::DOMAIN
            );
            assert_eq!(
                <DerivedExhaustiveStratifiedMapCandidateIdentitySchemaV1 as CanonicalSchema>::FIELDS
                    .len(),
                8
            );
            assert_eq!(
                DERIVED_EXHAUSTIVE_STRATIFIED_MAP_CANDIDATE_IDENTITY_LIMITS_V1.max_fields(),
                8
            );

            let source = admitted_fixed_resolution_geometry(70, 80, 90, 1, cx);
            let target = admitted_fixed_resolution_geometry(73, 83, 93, 2, cx);
            let components = vec![admit_stratum_component(&source, &target, 180, cx)];
            let ir = exhaustive_stratified_map_candidate_ir(&source, &target, &components, 183);
            let first = admit_derived_exhaustive_stratified_map_candidate_v1(
                &ir,
                &source,
                &target,
                &components,
                cx,
            )
            .expect("valid exhaustive assembly candidate");
            let replay = admit_derived_exhaustive_stratified_map_candidate_v1(
                &ir,
                &source,
                &target,
                &components,
                cx,
            )
            .expect("deterministic assembly replay");

            assert_eq!(first, replay);
            assert_eq!(first.source_geometry(), source.id());
            assert_eq!(first.source_stratification(), source.ir().stratification.id);
            assert_eq!(first.target_geometry(), target.id());
            assert_eq!(first.target_stratification(), target.ir().stratification.id);
            assert_eq!(first.components(), ir.components.as_slice());
            assert_eq!(first.nominal_assembly(), ir.nominal_assembly);
            assert_eq!(
                first.nominal_constructibility(),
                ir.nominal_constructibility
            );
            assert_eq!(first.no_authority(), ir.no_authority);

            let alternate_components = vec![admit_stratum_component(&source, &target, 193, cx)];
            let alternate_ir = exhaustive_stratified_map_candidate_ir(
                &source,
                &target,
                &alternate_components,
                183,
            );
            let alternate = admit_derived_exhaustive_stratified_map_candidate_v1(
                &alternate_ir,
                &source,
                &target,
                &alternate_components,
                cx,
            )
            .expect("alternate component remains a structural candidate");
            assert_ne!(first.id(), alternate.id());

            for (field, changed_ir) in [
                (
                    "assembly",
                    DerivedExhaustiveStratifiedMapCandidateIrV1 {
                        nominal_assembly: DerivedStratifiedMapAssemblyIdV1::from_bytes([190; 32]),
                        ..ir.clone()
                    },
                ),
                (
                    "constructibility",
                    DerivedExhaustiveStratifiedMapCandidateIrV1 {
                        nominal_constructibility:
                            DerivedGlobalConstructibilityDeclarationIdV1::from_bytes([191; 32]),
                        ..ir.clone()
                    },
                ),
                (
                    "no-authority",
                    DerivedExhaustiveStratifiedMapCandidateIrV1 {
                        no_authority: DerivedNoClaimIdV1::from_bytes([192; 32]),
                        ..ir.clone()
                    },
                ),
            ] {
                let changed = admit_derived_exhaustive_stratified_map_candidate_v1(
                    &changed_ir,
                    &source,
                    &target,
                    &components,
                    cx,
                )
                .unwrap_or_else(|error| panic!("changed {field} remains structural: {error}"));
                assert_ne!(first.id(), changed.id(), "{field} must move identity");
            }
        });
    }

    #[test]
    #[allow(clippy::too_many_lines)] // Full canonical coverage and binding-refusal fixture.
    fn exhaustive_candidate_checks_full_canonical_source_coverage() {
        with_cx(false, |cx| {
            let mut source_ir = fixed_resolution_geometry_ir(70, 80, 90, 1);
            let first_source_stratum = source_ir.stratification.strata[0].id;
            let second_source_stratum = stratum_id(150);
            let mut second = source_ir.stratification.strata[0].clone();
            second.id = second_source_stratum;
            source_ir.stratification.strata.push(second);
            let source =
                admit_derived_geometry_v1(source_ir, DerivedAdmissionBudgetV1::STANDARD, cx)
                    .expect("valid two-stratum source");
            let target = admitted_fixed_resolution_geometry(73, 83, 93, 2, cx);
            let source_object = |stratum| DerivedStratumObjectV1 {
                geometry: source.id(),
                stratification: source.ir().stratification.id,
                stratum,
            };
            let target_object = sole_stratum_object(&target);
            let components = vec![
                admit_stratum_component_between(
                    &source,
                    &target,
                    source_object(first_source_stratum),
                    target_object,
                    193,
                    cx,
                ),
                admit_stratum_component_between(
                    &source,
                    &target,
                    source_object(second_source_stratum),
                    target_object,
                    196,
                    cx,
                ),
            ];
            let ir = exhaustive_stratified_map_candidate_ir(&source, &target, &components, 199);
            let admitted = admit_derived_exhaustive_stratified_map_candidate_v1(
                &ir,
                &source,
                &target,
                &components,
                cx,
            )
            .expect("complete canonical source coverage");
            assert_eq!(admitted.components().len(), 2);
            assert_eq!(
                admitted.components()[0].source_stratum,
                first_source_stratum
            );
            assert_eq!(
                admitted.components()[1].source_stratum,
                second_source_stratum
            );
            assert_eq!(
                admitted.components()[0].target_stratum,
                admitted.components()[1].target_stratum,
                "target coverage and uniqueness are intentionally not required"
            );

            let mut missing = ir.clone();
            missing.components.pop();
            assert!(matches!(
                admit_derived_exhaustive_stratified_map_candidate_v1(
                    &missing,
                    &source,
                    &target,
                    &components[..1],
                    cx
                ),
                Err(
                    DerivedExhaustiveStratifiedMapCandidateErrorV1::ComponentCountMismatch {
                        field: "component-bindings",
                        expected: 2,
                        found: 1
                    }
                )
            ));

            let mut wrong_order = ir.clone();
            wrong_order.components[0].source_stratum = second_source_stratum;
            assert!(matches!(
                admit_derived_exhaustive_stratified_map_candidate_v1(
                    &wrong_order,
                    &source,
                    &target,
                    &components,
                    cx
                ),
                Err(
                    DerivedExhaustiveStratifiedMapCandidateErrorV1::ComponentEndpointMismatch {
                        index: 0,
                        field: "binding-source-stratum-order"
                    }
                )
            ));

            let mut wrong_id = ir.clone();
            wrong_id.components[0].component = components[1].id();
            assert_eq!(
                admit_derived_exhaustive_stratified_map_candidate_v1(
                    &wrong_id,
                    &source,
                    &target,
                    &components,
                    cx
                ),
                Err(
                    DerivedExhaustiveStratifiedMapCandidateErrorV1::ComponentIdentityMismatch {
                        index: 0
                    }
                )
            );
        });
    }

    #[test]
    fn exhaustive_candidate_accepts_identity_components_but_not_composite_paths() {
        with_cx(false, |cx| {
            let geometry = admitted_fixed_resolution_geometry(70, 80, 90, 1, cx);
            let identity = identity_derived_stratum_morphism_v1(
                &geometry,
                sole_stratum_object(&geometry).stratum,
                cx,
            )
            .expect("valid stratum identity");
            let identities = vec![identity];
            let identity_ir =
                exhaustive_stratified_map_candidate_ir(&geometry, &geometry, &identities, 202);
            assert!(
                admit_derived_exhaustive_stratified_map_candidate_v1(
                    &identity_ir,
                    &geometry,
                    &geometry,
                    &identities,
                    cx,
                )
                .is_ok()
            );

            let source = admitted_fixed_resolution_geometry(73, 83, 93, 2, cx);
            let middle = admitted_fixed_resolution_geometry(76, 86, 96, 3, cx);
            let target = admitted_fixed_resolution_geometry(79, 89, 99, 4, cx);
            let first = admit_stratum_component(&source, &middle, 205, cx);
            let second = admit_stratum_component(&middle, &target, 208, cx);
            let composite =
                compose_derived_stratum_morphisms_v1(&first, &second, cx).expect("sealed path");
            let composites = vec![composite];
            let composite_ir =
                exhaustive_stratified_map_candidate_ir(&source, &target, &composites, 211);
            assert_eq!(
                admit_derived_exhaustive_stratified_map_candidate_v1(
                    &composite_ir,
                    &source,
                    &target,
                    &composites,
                    cx,
                ),
                Err(
                    DerivedExhaustiveStratifiedMapCandidateErrorV1::CompositeComponent { index: 0 }
                )
            );
        });
    }

    #[test]
    fn exhaustive_candidate_refuses_missing_authority_fields_caps_and_cancellation() {
        let (source, target, components, ir) = with_cx(false, |cx| {
            let source = admitted_fixed_resolution_geometry(70, 80, 90, 1, cx);
            let target = admitted_fixed_resolution_geometry(73, 83, 93, 2, cx);
            let components = vec![admit_stratum_component(&source, &target, 214, cx)];
            let ir = exhaustive_stratified_map_candidate_ir(&source, &target, &components, 217);
            (source, target, components, ir)
        });
        with_cx(false, |cx| {
            for (field, changed_ir) in [
                (
                    "nominal-assembly",
                    DerivedExhaustiveStratifiedMapCandidateIrV1 {
                        nominal_assembly: DerivedStratifiedMapAssemblyIdV1::from_bytes([0; 32]),
                        ..ir.clone()
                    },
                ),
                (
                    "nominal-global-constructibility",
                    DerivedExhaustiveStratifiedMapCandidateIrV1 {
                        nominal_constructibility:
                            DerivedGlobalConstructibilityDeclarationIdV1::from_bytes([0; 32]),
                        ..ir.clone()
                    },
                ),
                (
                    "no-global-map-authority",
                    DerivedExhaustiveStratifiedMapCandidateIrV1 {
                        no_authority: DerivedNoClaimIdV1::from_bytes([0; 32]),
                        ..ir.clone()
                    },
                ),
            ] {
                assert!(matches!(
                    admit_derived_exhaustive_stratified_map_candidate_v1(
                        &changed_ir,
                        &source,
                        &target,
                        &components,
                        cx
                    ),
                    Err(
                        DerivedExhaustiveStratifiedMapCandidateErrorV1::MissingIdentity {
                            field: found
                        }
                    ) if found == field
                ));
            }

            let mut oversized = ir.clone();
            oversized.components = vec![ir.components[0]; DERIVED_MORPHISM_MAX_FACTORS_V1 + 1];
            assert!(matches!(
                admit_derived_exhaustive_stratified_map_candidate_v1(
                    &oversized,
                    &source,
                    &target,
                    &components,
                    cx
                ),
                Err(
                    DerivedExhaustiveStratifiedMapCandidateErrorV1::ResourceLimit {
                        field: "component-bindings",
                        requested,
                        limit: DERIVED_MORPHISM_MAX_FACTORS_V1
                    }
                ) if requested == DERIVED_MORPHISM_MAX_FACTORS_V1 + 1
            ));
        });
        with_cx(true, |cx| {
            assert!(matches!(
                admit_derived_exhaustive_stratified_map_candidate_v1(
                    &ir,
                    &source,
                    &target,
                    &components,
                    cx
                ),
                Err(DerivedExhaustiveStratifiedMapCandidateErrorV1::Cancelled {
                    stage: "stratified-assembly-admission-entry"
                })
            ));
        });
    }

    #[test]
    fn refinement_truncation_order_is_fail_closed_at_untruncated_zero() {
        assert_eq!(truncation_refinement_progress(0, 0), Some(false));
        assert_eq!(truncation_refinement_progress(0, 1), None);
        assert_eq!(truncation_refinement_progress(2, 1), None);
        assert_eq!(truncation_refinement_progress(2, 2), Some(false));
        assert_eq!(truncation_refinement_progress(2, 3), Some(true));
        assert_eq!(truncation_refinement_progress(2, 0), Some(true));
    }

    #[test]
    fn identity_is_neutral_and_composition_is_associative_by_receipt() {
        with_cx(false, |cx| {
            let x = endpoint(10);
            let y = endpoint(11);
            let z = endpoint(12);
            let w = endpoint(13);
            let f = admit_strict(x, y, 20, ColorRank::Verified, ColorRank::Validated, cx);
            let g = admit_strict(y, z, 21, ColorRank::Validated, ColorRank::Estimated, cx);
            let h = admit_strict(z, w, 22, ColorRank::Estimated, ColorRank::Estimated, cx);
            let left_identity = admit_identity(x, cx);
            let right_identity = admit_identity(w, cx);
            assert_eq!(left_identity, admit_identity(x, cx));
            assert_eq!(
                compose_derived_morphisms_v1(&left_identity, &f, cx).expect("left identity"),
                f
            );
            assert_eq!(
                compose_derived_morphisms_v1(&h, &right_identity, cx).expect("right identity"),
                h
            );

            let fg = compose_derived_morphisms_v1(&f, &g, cx).expect("f then g");
            let gh = compose_derived_morphisms_v1(&g, &h, cx).expect("g then h");
            let left = compose_derived_morphisms_v1(&fg, &h, cx).expect("(fg)h");
            let right = compose_derived_morphisms_v1(&f, &gh, cx).expect("f(gh)");
            assert_eq!(left, right);
            assert_eq!(left.primitive_factors().len(), 3);
            assert_eq!(left.no_equivalence_claims().len(), 3);
            assert_eq!(
                left.evidence().ranks(),
                Some((ColorRank::Verified, ColorRank::Estimated))
            );
        });
    }

    #[test]
    fn factor_order_is_semantic_and_middle_endpoint_is_exact() {
        with_cx(false, |cx| {
            let x = endpoint(30);
            let a = admit_strict(x, x, 31, ColorRank::Validated, ColorRank::Validated, cx);
            let b = admit_strict(x, x, 32, ColorRank::Validated, ColorRank::Validated, cx);
            assert_eq!(
                a.primitive_path(),
                &[AdmittedDerivedPrimitiveV1::Strict {
                    source_geometry: x.id,
                    target_geometry: x.id,
                    witness: DerivedWitnessIdV1::from_bytes([31; 32]),
                }]
            );
            let ab = compose_derived_morphisms_v1(&a, &b, cx).expect("a then b");
            let ba = compose_derived_morphisms_v1(&b, &a, cx).expect("b then a");
            assert_ne!(ab.id(), ba.id());
            assert_eq!(ab.primitive_factors(), &[a.id(), b.id()]);
            assert_eq!(
                ab.no_equivalence_claims(),
                &[a.no_equivalence_claims()[0], b.no_equivalence_claims()[0]]
            );

            let y = endpoint(33);
            let z = endpoint(34);
            let yz = admit_strict(y, z, 35, ColorRank::Validated, ColorRank::Validated, cx);
            assert_eq!(
                compose_derived_morphisms_v1(&a, &yz, cx),
                Err(DerivedMorphismErrorV1::CompositionEndpointMismatch)
            );
        });
    }

    #[test]
    fn evidence_cannot_reverse_variance_or_strengthen() {
        with_cx(false, |cx| {
            let x = endpoint(40);
            let y = endpoint(41);
            let strengthening = strict_ir(x, y, 42, ColorRank::Estimated, ColorRank::Verified);
            assert!(matches!(
                admit_between_endpoints(strengthening, x, y, cx),
                Err(DerivedMorphismErrorV1::EvidenceStrengthening { .. })
            ));

            let mut reversed = restriction_ir(x, y, 43, ColorRank::Verified, ColorRank::Validated);
            reversed.evidence = DerivedEvidenceTransportV1::RestrictionContravariant {
                input_geometry: x.id,
                output_geometry: y.id,
                input_evidence: evidence_id(x.id),
                output_evidence: evidence_id(y.id),
                input_rank: ColorRank::Verified,
                output_rank: ColorRank::Validated,
            };
            assert_eq!(
                admit_between_endpoints(reversed, x, y, cx),
                Err(DerivedMorphismErrorV1::EvidenceOrientationMismatch)
            );

            let mut missing = strict_ir(x, y, 44, ColorRank::Validated, ColorRank::Estimated);
            missing.evidence = DerivedEvidenceTransportV1::BalanceCorestrictionCovariant {
                input_geometry: x.id,
                output_geometry: y.id,
                input_evidence: DerivedEvidenceArtifactIdV1::from_bytes([0; 32]),
                output_evidence: evidence_id(y.id),
                input_rank: ColorRank::Validated,
                output_rank: ColorRank::Estimated,
            };
            assert_eq!(
                admit_between_endpoints(missing, x, y, cx),
                Err(DerivedMorphismErrorV1::MissingIdentity {
                    field: "input-evidence"
                })
            );
        });
    }

    #[test]
    fn contravariant_restrictions_compose_in_evidence_order() {
        with_cx(false, |cx| {
            let x = endpoint(44);
            let y = endpoint(45);
            let z = endpoint(46);
            let first = admit_between_endpoints(
                restriction_ir(x, y, 47, ColorRank::Verified, ColorRank::Validated),
                x,
                y,
                cx,
            )
            .expect("restriction Y to X");
            let second = admit_between_endpoints(
                restriction_ir(y, z, 48, ColorRank::Verified, ColorRank::Verified),
                y,
                z,
                cx,
            )
            .expect("restriction Z to Y");
            let composite =
                compose_derived_morphisms_v1(&first, &second, cx).expect("restriction Z to X");
            assert_eq!(
                composite.evidence(),
                DerivedEvidenceTransportV1::RestrictionContravariant {
                    input_geometry: z.id,
                    output_geometry: x.id,
                    input_evidence: evidence_id(z.id),
                    output_evidence: evidence_id(x.id),
                    input_rank: ColorRank::Verified,
                    output_rank: ColorRank::Validated,
                }
            );
        });
    }

    #[test]
    fn declared_chart_maps_compose_as_associative_typed_paths() {
        with_cx(false, |cx| {
            let x_charts = [chart(90, 2, 2, 9, 1.0)];
            let y_charts = [chart(91, 2, 2, 9, 1.0)];
            let z_charts = [chart(92, 2, 2, 9, 1.0)];
            let w_charts = [chart(93, 2, 2, 9, 1.0)];
            let x = endpoint_with_charts(86, &x_charts);
            let y = endpoint_with_charts(87, &y_charts);
            let z = endpoint_with_charts(88, &z_charts);
            let w = endpoint_with_charts(89, &w_charts);
            let f = admit_chart_map(
                x,
                y,
                x_charts[0].id,
                y_charts[0].id,
                100,
                ColorRank::Verified,
                ColorRank::Validated,
                cx,
            );
            let g = admit_chart_map(
                y,
                z,
                y_charts[0].id,
                z_charts[0].id,
                104,
                ColorRank::Validated,
                ColorRank::Validated,
                cx,
            );
            let h = admit_chart_map(
                z,
                w,
                z_charts[0].id,
                w_charts[0].id,
                108,
                ColorRank::Validated,
                ColorRank::Estimated,
                cx,
            );

            let fg = compose_derived_morphisms_v1(&f, &g, cx).expect("f then g");
            let gh = compose_derived_morphisms_v1(&g, &h, cx).expect("g then h");
            let left = compose_derived_morphisms_v1(&fg, &h, cx).expect("(fg)h");
            let right = compose_derived_morphisms_v1(&f, &gh, cx).expect("f(gh)");
            assert_eq!(left, right);
            assert_eq!(
                left.class(),
                AdmittedDerivedMorphismClassV1::DeclaredChartMapPath
            );
            assert_eq!(
                left.chart_path(),
                Some(DerivedChartPathV1 {
                    source_chart: x_charts[0].id,
                    target_chart: w_charts[0].id,
                })
            );
            assert_eq!(left.primitive_factors(), &[f.id(), g.id(), h.id()]);
            assert_eq!(
                left.primitive_path(),
                &[
                    f.primitive_path()[0],
                    g.primitive_path()[0],
                    h.primitive_path()[0],
                ]
            );
            assert_eq!(
                left.declared_chart_maps(),
                &[
                    f.declared_chart_maps()[0],
                    g.declared_chart_maps()[0],
                    h.declared_chart_maps()[0],
                ]
            );
            let identity = admit_identity(x, cx);
            assert_eq!(
                compose_derived_morphisms_v1(&identity, &f, cx).expect("identity then f"),
                f
            );
        });
    }

    #[test]
    fn mixed_primitive_families_form_associative_typed_paths() {
        with_cx(false, |cx| {
            let y_charts = [chart(160, 2, 2, 16, 1.0)];
            let z_charts = [chart(161, 2, 2, 16, 1.0)];
            let x = endpoint(162);
            let y = endpoint_with_charts(163, &y_charts);
            let z = endpoint_with_charts(164, &z_charts);
            let w = endpoint(165);
            let f = admit_strict(x, y, 166, ColorRank::Verified, ColorRank::Validated, cx);
            let g = admit_chart_map(
                y,
                z,
                y_charts[0].id,
                z_charts[0].id,
                167,
                ColorRank::Validated,
                ColorRank::Validated,
                cx,
            );
            let h = admit_strict(z, w, 170, ColorRank::Validated, ColorRank::Estimated, cx);

            let fg = compose_derived_morphisms_v1(&f, &g, cx).expect("strict then chart");
            let gh = compose_derived_morphisms_v1(&g, &h, cx).expect("chart then strict");
            let left = compose_derived_morphisms_v1(&fg, &h, cx).expect("(fg)h");
            let right = compose_derived_morphisms_v1(&f, &gh, cx).expect("f(gh)");
            assert_eq!(left, right);
            assert_eq!(
                left.class(),
                AdmittedDerivedMorphismClassV1::HeterogeneousPath
            );
            assert_eq!(left.chart_path(), None);
            assert_eq!(left.primitive_factors(), &[f.id(), g.id(), h.id()]);
            assert_eq!(
                left.primitive_path(),
                &[
                    f.primitive_path()[0],
                    g.primitive_path()[0],
                    h.primitive_path()[0],
                ]
            );
            assert!(matches!(
                left.primitive_path(),
                [
                    AdmittedDerivedPrimitiveV1::Strict { .. },
                    AdmittedDerivedPrimitiveV1::DeclaredChartMap(_),
                    AdmittedDerivedPrimitiveV1::Strict { .. }
                ]
            ));
            assert_eq!(left.declared_chart_maps(), g.declared_chart_maps());
            let left_identity = admit_identity(x, cx);
            let right_identity = admit_identity(w, cx);
            assert_eq!(
                compose_derived_morphisms_v1(&left_identity, &left, cx)
                    .expect("heterogeneous left identity"),
                left
            );
            assert_eq!(
                compose_derived_morphisms_v1(&left, &right_identity, cx)
                    .expect("heterogeneous right identity"),
                left
            );
        });
    }

    #[test]
    fn heterogeneous_paths_preserve_adjacent_chart_seam_checks() {
        with_cx(false, |cx| {
            let y_charts = [chart(180, 2, 2, 18, 1.0)];
            let z_charts = [chart(181, 2, 2, 18, 1.0), chart(182, 2, 2, 18, 1.0)];
            let w_charts = [chart(183, 2, 2, 18, 1.0)];
            let x = endpoint(184);
            let y = endpoint_with_charts(185, &y_charts);
            let z = endpoint_with_charts(186, &z_charts);
            let w = endpoint_with_charts(187, &w_charts);
            let strict = admit_strict(x, y, 188, ColorRank::Verified, ColorRank::Validated, cx);
            let first_chart = admit_chart_map(
                y,
                z,
                y_charts[0].id,
                z_charts[0].id,
                189,
                ColorRank::Validated,
                ColorRank::Validated,
                cx,
            );
            let mixed =
                compose_derived_morphisms_v1(&strict, &first_chart, cx).expect("mixed prefix");
            let wrong_chart = admit_chart_map(
                z,
                w,
                z_charts[1].id,
                w_charts[0].id,
                190,
                ColorRank::Validated,
                ColorRank::Estimated,
                cx,
            );
            assert_eq!(
                compose_derived_morphisms_v1(&mixed, &wrong_chart, cx),
                Err(DerivedMorphismErrorV1::CompositionChartMismatch)
            );
        });
    }

    #[test]
    #[allow(clippy::too_many_lines)] // Full three-level refinement law and retained-lineage fixture.
    fn finite_complex_refinements_compose_as_associative_typed_paths() {
        with_cx(false, |cx| {
            let charts = [chart(190, 2, 2, 21, 1.0)];
            let x_complexes = [complex(
                200,
                210,
                charts[0].id,
                DerivedComplexRoleV1::Tangent,
                &[(0, 1, 21), (1, 1, 22)],
                2,
            )];
            let y_complexes = [complex(
                201,
                211,
                charts[0].id,
                DerivedComplexRoleV1::Tangent,
                &[(0, 2, 21), (1, 1, 22)],
                2,
            )];
            let z_complexes = [complex(
                202,
                212,
                charts[0].id,
                DerivedComplexRoleV1::Tangent,
                &[(0, 2, 21), (1, 2, 22)],
                3,
            )];
            let w_complexes = [complex(
                203,
                213,
                charts[0].id,
                DerivedComplexRoleV1::Tangent,
                &[(-1, 1, 23), (0, 2, 21), (1, 2, 22)],
                0,
            )];
            let x = endpoint_with_parts(220, &charts, &x_complexes);
            let y = endpoint_with_parts(221, &charts, &y_complexes);
            let z = endpoint_with_parts(222, &charts, &z_complexes);
            let w = endpoint_with_parts(223, &charts, &w_complexes);
            let f = admit_complex_refinement(
                x,
                y,
                x_complexes[0].id,
                y_complexes[0].id,
                x_complexes[0].resolution.id,
                y_complexes[0].resolution.id,
                120,
                ColorRank::Verified,
                ColorRank::Validated,
                cx,
            );
            let g = admit_complex_refinement(
                y,
                z,
                y_complexes[0].id,
                z_complexes[0].id,
                y_complexes[0].resolution.id,
                z_complexes[0].resolution.id,
                124,
                ColorRank::Validated,
                ColorRank::Validated,
                cx,
            );
            let h = admit_complex_refinement(
                z,
                w,
                z_complexes[0].id,
                w_complexes[0].id,
                z_complexes[0].resolution.id,
                w_complexes[0].resolution.id,
                128,
                ColorRank::Validated,
                ColorRank::Estimated,
                cx,
            );

            let fg = compose_derived_morphisms_v1(&f, &g, cx).expect("f then g");
            let gh = compose_derived_morphisms_v1(&g, &h, cx).expect("g then h");
            let left = compose_derived_morphisms_v1(&fg, &h, cx).expect("(fg)h");
            let right = compose_derived_morphisms_v1(&f, &gh, cx).expect("f(gh)");
            assert_eq!(left, right);
            assert_eq!(
                left.class(),
                AdmittedDerivedMorphismClassV1::DeclaredComplexRefinementPath
            );
            assert_eq!(left.primitive_factors(), &[f.id(), g.id(), h.id()]);
            assert_eq!(
                left.primitive_path(),
                &[
                    f.primitive_path()[0],
                    g.primitive_path()[0],
                    h.primitive_path()[0],
                ]
            );
            assert!(left.primitive_path().iter().all(|primitive| matches!(
                primitive,
                AdmittedDerivedPrimitiveV1::DeclaredComplexRefinement(_)
            )));
            let left_identity = admit_identity(x, cx);
            let right_identity = admit_identity(w, cx);
            assert_eq!(
                compose_derived_morphisms_v1(&left_identity, &left, cx)
                    .expect("refinement left identity"),
                left
            );
            assert_eq!(
                compose_derived_morphisms_v1(&left, &right_identity, cx)
                    .expect("refinement right identity"),
                left
            );
        });
    }

    #[test]
    #[allow(clippy::too_many_lines)] // Replay plus independent movement of every receipt selector.
    fn complex_refinement_receipt_binds_typed_nominal_artifacts() {
        with_cx(false, |cx| {
            let charts = [chart(170, 2, 2, 17, 1.0)];
            let source_complexes = [complex(
                171,
                173,
                charts[0].id,
                DerivedComplexRoleV1::Tangent,
                &[(0, 1, 17), (1, 1, 18)],
                1,
            )];
            let target_complexes = [complex(
                172,
                174,
                charts[0].id,
                DerivedComplexRoleV1::Tangent,
                &[(0, 2, 17), (1, 1, 18)],
                1,
            )];
            let source = endpoint_with_parts(175, &charts, &source_complexes);
            let target = endpoint_with_parts(176, &charts, &target_complexes);
            let base_ir = complex_refinement_ir(
                source,
                target,
                source_complexes[0].id,
                target_complexes[0].id,
                source_complexes[0].resolution.id,
                target_complexes[0].resolution.id,
                50,
                ColorRank::Validated,
                ColorRank::Estimated,
            );
            let base = admit_between_endpoints(base_ir, source, target, cx)
                .expect("base complex refinement");
            let replay = admit_between_endpoints(base_ir, source, target, cx)
                .expect("replayed complex refinement");
            assert_eq!(base, replay);
            assert_eq!(
                base.primitive_path(),
                &[AdmittedDerivedPrimitiveV1::DeclaredComplexRefinement(
                    DeclaredComplexRefinementPrimitiveV1 {
                        source_geometry: source.id,
                        target_geometry: target.id,
                        source_complex: source_complexes[0].id,
                        target_complex: target_complexes[0].id,
                        source_resolution: source_complexes[0].resolution.id,
                        target_resolution: target_complexes[0].resolution.id,
                        prolongation: DerivedComplexRefinementMapIdV1::from_bytes([51; 32]),
                        commutation: DerivedWitnessIdV1::from_bytes([52; 32]),
                    }
                )]
            );

            let mut changed_prolongation = base_ir;
            if let DerivedMorphismKindV1::DeclaredComplexRefinement { prolongation, .. } =
                &mut changed_prolongation.kind
            {
                *prolongation = DerivedComplexRefinementMapIdV1::from_bytes([53; 32]);
            }
            let changed_prolongation =
                admit_between_endpoints(changed_prolongation, source, target, cx)
                    .expect("changed prolongation remains a declaration");
            assert_ne!(base.id(), changed_prolongation.id());

            let mut changed_commutation = base_ir;
            if let DerivedMorphismKindV1::DeclaredComplexRefinement { commutation, .. } =
                &mut changed_commutation.kind
            {
                *commutation = DerivedWitnessIdV1::from_bytes([54; 32]);
            }
            let changed_commutation =
                admit_between_endpoints(changed_commutation, source, target, cx)
                    .expect("changed commutation remains a declaration");
            assert_ne!(base.id(), changed_commutation.id());

            let mut changed_source_complex = source_complexes[0].clone();
            changed_source_complex.id = complex_id(177);
            let changed_source_complexes = [changed_source_complex];
            let changed_source = endpoint_with_parts(175, &charts, &changed_source_complexes);
            let mut changed_source_complex_ir = base_ir;
            if let DerivedMorphismKindV1::DeclaredComplexRefinement { source_complex, .. } =
                &mut changed_source_complex_ir.kind
            {
                *source_complex = changed_source_complexes[0].id;
            }
            let changed_source_complex =
                admit_between_endpoints(changed_source_complex_ir, changed_source, target, cx)
                    .expect("changed source complex remains a declaration");
            assert_ne!(base.id(), changed_source_complex.id());

            let mut changed_target_complex = target_complexes[0].clone();
            changed_target_complex.id = complex_id(178);
            let changed_target_complexes = [changed_target_complex];
            let changed_target = endpoint_with_parts(176, &charts, &changed_target_complexes);
            let mut changed_target_complex_ir = base_ir;
            if let DerivedMorphismKindV1::DeclaredComplexRefinement { target_complex, .. } =
                &mut changed_target_complex_ir.kind
            {
                *target_complex = changed_target_complexes[0].id;
            }
            let changed_target_complex =
                admit_between_endpoints(changed_target_complex_ir, source, changed_target, cx)
                    .expect("changed target complex remains a declaration");
            assert_ne!(base.id(), changed_target_complex.id());

            let mut changed_source_resolution = source_complexes[0].clone();
            changed_source_resolution.resolution.id = resolution_id(179);
            let changed_source_resolutions = [changed_source_resolution];
            let changed_source = endpoint_with_parts(175, &charts, &changed_source_resolutions);
            let mut changed_source_resolution_ir = base_ir;
            if let DerivedMorphismKindV1::DeclaredComplexRefinement {
                source_resolution, ..
            } = &mut changed_source_resolution_ir.kind
            {
                *source_resolution = changed_source_resolutions[0].resolution.id;
            }
            let changed_source_resolution =
                admit_between_endpoints(changed_source_resolution_ir, changed_source, target, cx)
                    .expect("changed source resolution remains a declaration");
            assert_ne!(base.id(), changed_source_resolution.id());

            let mut changed_target_resolution = target_complexes[0].clone();
            changed_target_resolution.resolution.id = resolution_id(180);
            let changed_target_resolutions = [changed_target_resolution];
            let changed_target = endpoint_with_parts(176, &charts, &changed_target_resolutions);
            let mut changed_target_resolution_ir = base_ir;
            if let DerivedMorphismKindV1::DeclaredComplexRefinement {
                target_resolution, ..
            } = &mut changed_target_resolution_ir.kind
            {
                *target_resolution = changed_target_resolutions[0].resolution.id;
            }
            let changed_target_resolution =
                admit_between_endpoints(changed_target_resolution_ir, source, changed_target, cx)
                    .expect("changed target resolution remains a declaration");
            assert_ne!(base.id(), changed_target_resolution.id());
        });
    }

    #[test]
    #[allow(clippy::too_many_lines)] // One mutation table for selector and authority boundaries.
    fn complex_refinement_refuses_missing_selectors_and_authority_laundering() {
        with_cx(false, |cx| {
            let charts = [chart(40, 2, 2, 4, 1.0)];
            let source_complexes = [complex(
                60,
                80,
                charts[0].id,
                DerivedComplexRoleV1::Tangent,
                &[(0, 1, 4), (1, 1, 5)],
                1,
            )];
            let target_complexes = [complex(
                61,
                81,
                charts[0].id,
                DerivedComplexRoleV1::Tangent,
                &[(0, 2, 4), (1, 1, 5)],
                1,
            )];
            let source = endpoint_with_parts(90, &charts, &source_complexes);
            let target = endpoint_with_parts(91, &charts, &target_complexes);
            let base = complex_refinement_ir(
                source,
                target,
                source_complexes[0].id,
                target_complexes[0].id,
                source_complexes[0].resolution.id,
                target_complexes[0].resolution.id,
                100,
                ColorRank::Validated,
                ColorRank::Estimated,
            );

            let mut unknown_complex = base;
            if let DerivedMorphismKindV1::DeclaredComplexRefinement { source_complex, .. } =
                &mut unknown_complex.kind
            {
                *source_complex = complex_id(62);
            }
            assert_eq!(
                admit_between_endpoints(unknown_complex, source, target, cx),
                Err(DerivedMorphismErrorV1::MissingComplex {
                    field: "source-complex"
                })
            );

            let mut unknown_target_complex = base;
            if let DerivedMorphismKindV1::DeclaredComplexRefinement { target_complex, .. } =
                &mut unknown_target_complex.kind
            {
                *target_complex = complex_id(62);
            }
            assert_eq!(
                admit_between_endpoints(unknown_target_complex, source, target, cx),
                Err(DerivedMorphismErrorV1::MissingComplex {
                    field: "target-complex"
                })
            );

            let mut wrong_resolution = base;
            if let DerivedMorphismKindV1::DeclaredComplexRefinement {
                source_resolution, ..
            } = &mut wrong_resolution.kind
            {
                *source_resolution = resolution_id(82);
            }
            assert_eq!(
                admit_between_endpoints(wrong_resolution, source, target, cx),
                Err(DerivedMorphismErrorV1::ComplexRefinementMismatch {
                    field: "source-resolution"
                })
            );

            let mut wrong_target_resolution = base;
            if let DerivedMorphismKindV1::DeclaredComplexRefinement {
                target_resolution, ..
            } = &mut wrong_target_resolution.kind
            {
                *target_resolution = resolution_id(82);
            }
            assert_eq!(
                admit_between_endpoints(wrong_target_resolution, source, target, cx),
                Err(DerivedMorphismErrorV1::ComplexRefinementMismatch {
                    field: "target-resolution"
                })
            );

            let mut zero_complex = base;
            if let DerivedMorphismKindV1::DeclaredComplexRefinement { source_complex, .. } =
                &mut zero_complex.kind
            {
                *source_complex = complex_id(0);
            }
            assert_eq!(
                admit_between_endpoints(zero_complex, source, target, cx),
                Err(DerivedMorphismErrorV1::MissingIdentity {
                    field: "source-complex"
                })
            );

            let mut zero_target_complex = base;
            if let DerivedMorphismKindV1::DeclaredComplexRefinement { target_complex, .. } =
                &mut zero_target_complex.kind
            {
                *target_complex = complex_id(0);
            }
            assert_eq!(
                admit_between_endpoints(zero_target_complex, source, target, cx),
                Err(DerivedMorphismErrorV1::MissingIdentity {
                    field: "target-complex"
                })
            );

            let mut zero_source_resolution = base;
            if let DerivedMorphismKindV1::DeclaredComplexRefinement {
                source_resolution, ..
            } = &mut zero_source_resolution.kind
            {
                *source_resolution = resolution_id(0);
            }
            assert_eq!(
                admit_between_endpoints(zero_source_resolution, source, target, cx),
                Err(DerivedMorphismErrorV1::MissingIdentity {
                    field: "source-resolution"
                })
            );

            let mut zero_target_resolution = base;
            if let DerivedMorphismKindV1::DeclaredComplexRefinement {
                target_resolution, ..
            } = &mut zero_target_resolution.kind
            {
                *target_resolution = resolution_id(0);
            }
            assert_eq!(
                admit_between_endpoints(zero_target_resolution, source, target, cx),
                Err(DerivedMorphismErrorV1::MissingIdentity {
                    field: "target-resolution"
                })
            );

            let mut zero_prolongation = base;
            if let DerivedMorphismKindV1::DeclaredComplexRefinement { prolongation, .. } =
                &mut zero_prolongation.kind
            {
                *prolongation = DerivedComplexRefinementMapIdV1::from_bytes([0; 32]);
            }
            assert_eq!(
                admit_between_endpoints(zero_prolongation, source, target, cx),
                Err(DerivedMorphismErrorV1::MissingIdentity {
                    field: "complex-refinement-prolongation"
                })
            );

            let mut zero_commutation = base;
            if let DerivedMorphismKindV1::DeclaredComplexRefinement { commutation, .. } =
                &mut zero_commutation.kind
            {
                *commutation = DerivedWitnessIdV1::from_bytes([0; 32]);
            }
            assert_eq!(
                admit_between_endpoints(zero_commutation, source, target, cx),
                Err(DerivedMorphismErrorV1::MissingIdentity {
                    field: "complex-refinement-commutation"
                })
            );

            let mut zero_no_claim = base;
            if let DerivedEquivalenceBoundaryV1::NoClaim { artifact } =
                &mut zero_no_claim.equivalence
            {
                *artifact = DerivedNoClaimIdV1::from_bytes([0; 32]);
            }
            assert_eq!(
                admit_between_endpoints(zero_no_claim, source, target, cx),
                Err(DerivedMorphismErrorV1::MissingIdentity {
                    field: "no-equivalence-claim"
                })
            );

            let mut identity_evidence = base;
            identity_evidence.evidence = DerivedEvidenceTransportV1::Identity;
            assert_eq!(
                admit_between_endpoints(identity_evidence, source, target, cx),
                Err(DerivedMorphismErrorV1::EvidenceOrientationMismatch)
            );

            let mut laundering = base;
            laundering.equivalence = DerivedEquivalenceBoundaryV1::IdentityOnly;
            assert_eq!(
                admit_between_endpoints(laundering, source, target, cx),
                Err(DerivedMorphismErrorV1::EquivalenceLaundering)
            );
        });
    }

    #[test]
    #[allow(clippy::too_many_lines)] // Independent mutations for every rank-envelope rule.
    fn complex_refinement_refuses_role_chart_and_shape_regressions() {
        with_cx(false, |cx| {
            let charts = [chart(41, 2, 2, 6, 1.0)];
            let source_complexes = [complex(
                70,
                90,
                charts[0].id,
                DerivedComplexRoleV1::Tangent,
                &[(0, 2, 6), (1, 1, 7)],
                2,
            )];
            let source = endpoint_with_parts(110, &charts, &source_complexes);

            let wrong_role = [complex(
                71,
                91,
                charts[0].id,
                DerivedComplexRoleV1::Cotangent,
                &[(0, 3, 6), (1, 1, 7)],
                2,
            )];
            let wrong_role_target = endpoint_with_parts(111, &charts, &wrong_role);
            assert_eq!(
                try_complex_refinement(
                    source,
                    wrong_role_target,
                    &source_complexes[0],
                    &wrong_role[0],
                    cx,
                ),
                Err(DerivedMorphismErrorV1::ComplexRefinementMismatch {
                    field: "complex-role"
                })
            );

            let other_charts = [chart(42, 2, 2, 6, 1.0)];
            let wrong_chart = [complex(
                72,
                92,
                other_charts[0].id,
                DerivedComplexRoleV1::Tangent,
                &[(0, 3, 6), (1, 1, 7)],
                2,
            )];
            let wrong_chart_target = endpoint_with_parts(112, &other_charts, &wrong_chart);
            assert_eq!(
                try_complex_refinement(
                    source,
                    wrong_chart_target,
                    &source_complexes[0],
                    &wrong_chart[0],
                    cx,
                ),
                Err(DerivedMorphismErrorV1::ComplexRefinementMismatch {
                    field: "complex-chart"
                })
            );

            let mut changed_chart = charts[0].clone();
            changed_chart.coordinates.scale_to_canonical = 2.0;
            let changed_charts = [changed_chart];
            let changed_chart_complexes = [complex(
                73,
                93,
                charts[0].id,
                DerivedComplexRoleV1::Tangent,
                &[(0, 3, 6), (1, 1, 7)],
                2,
            )];
            let changed_chart_target =
                endpoint_with_parts(113, &changed_charts, &changed_chart_complexes);
            assert_eq!(
                try_complex_refinement(
                    source,
                    changed_chart_target,
                    &source_complexes[0],
                    &changed_chart_complexes[0],
                    cx,
                ),
                Err(DerivedMorphismErrorV1::ComplexRefinementMismatch {
                    field: "chart-semantics"
                })
            );

            for (target_complex, expected_field) in [
                (
                    complex(
                        74,
                        94,
                        charts[0].id,
                        DerivedComplexRoleV1::Tangent,
                        &[(0, 1, 6), (1, 1, 7)],
                        2,
                    ),
                    "degree-rank",
                ),
                (
                    complex(
                        75,
                        95,
                        charts[0].id,
                        DerivedComplexRoleV1::Tangent,
                        &[(0, 3, 8), (1, 1, 7)],
                        2,
                    ),
                    "degree-quantity",
                ),
                (
                    complex(
                        76,
                        96,
                        charts[0].id,
                        DerivedComplexRoleV1::Tangent,
                        &[(0, 3, 6)],
                        2,
                    ),
                    "degree-coverage",
                ),
                (
                    complex(
                        67,
                        87,
                        charts[0].id,
                        DerivedComplexRoleV1::Tangent,
                        &[(1, 3, 6), (2, 1, 7)],
                        2,
                    ),
                    "degree-coverage",
                ),
                (
                    complex(
                        77,
                        97,
                        charts[0].id,
                        DerivedComplexRoleV1::Tangent,
                        &[(0, 2, 6), (1, 1, 7)],
                        2,
                    ),
                    "strict-progress",
                ),
                (
                    complex(
                        78,
                        98,
                        charts[0].id,
                        DerivedComplexRoleV1::Tangent,
                        &[(-1, 0, 8), (0, 2, 6), (1, 1, 7)],
                        2,
                    ),
                    "strict-progress",
                ),
                (
                    complex(
                        79,
                        99,
                        charts[0].id,
                        DerivedComplexRoleV1::Tangent,
                        &[(0, 3, 6), (1, 1, 7)],
                        1,
                    ),
                    "truncation-policy",
                ),
            ] {
                let target_complexes = [target_complex];
                let target = endpoint_with_parts(114, &charts, &target_complexes);
                assert_eq!(
                    try_complex_refinement(
                        source,
                        target,
                        &source_complexes[0],
                        &target_complexes[0],
                        cx,
                    ),
                    Err(DerivedMorphismErrorV1::ComplexRefinementMismatch {
                        field: expected_field
                    })
                );
            }

            for (seed, resolution_seed, truncation_order) in [(68, 88, 3), (69, 89, 0)] {
                let target_complexes = [complex(
                    seed,
                    resolution_seed,
                    charts[0].id,
                    DerivedComplexRoleV1::Tangent,
                    &[(0, 2, 6), (1, 1, 7)],
                    truncation_order,
                )];
                let target = endpoint_with_parts(115, &charts, &target_complexes);
                assert!(
                    try_complex_refinement(
                        source,
                        target,
                        &source_complexes[0],
                        &target_complexes[0],
                        cx,
                    )
                    .is_ok(),
                    "truncation-only improvement should admit"
                );
            }

            let trailing_degree = [complex(
                66,
                86,
                charts[0].id,
                DerivedComplexRoleV1::Tangent,
                &[(0, 2, 6), (1, 1, 7), (2, 1, 8)],
                2,
            )];
            let trailing_degree_target = endpoint_with_parts(116, &charts, &trailing_degree);
            assert!(
                try_complex_refinement(
                    source,
                    trailing_degree_target,
                    &source_complexes[0],
                    &trailing_degree[0],
                    cx,
                )
                .is_ok(),
                "positive trailing degree should establish strict progress"
            );
        });
    }

    #[test]
    #[allow(clippy::too_many_lines)] // Exact homogeneous and heterogeneous seam fixture.
    fn complex_refinement_seams_survive_heterogeneous_parenthesization() {
        with_cx(false, |cx| {
            let charts = [chart(150, 2, 2, 15, 1.0)];
            let x_complexes = [complex(
                151,
                160,
                charts[0].id,
                DerivedComplexRoleV1::Tangent,
                &[(0, 1, 15), (1, 1, 16)],
                1,
            )];
            let y_complexes = [
                complex(
                    152,
                    161,
                    charts[0].id,
                    DerivedComplexRoleV1::Tangent,
                    &[(0, 2, 15), (1, 1, 16)],
                    1,
                ),
                complex(
                    153,
                    162,
                    charts[0].id,
                    DerivedComplexRoleV1::Tangent,
                    &[(0, 2, 15), (1, 1, 16)],
                    1,
                ),
            ];
            let z_complexes = [complex(
                154,
                163,
                charts[0].id,
                DerivedComplexRoleV1::Tangent,
                &[(0, 2, 15), (1, 2, 16)],
                1,
            )];
            let w = endpoint(159);
            let x = endpoint_with_parts(155, &charts, &x_complexes);
            let y = endpoint_with_parts(156, &charts, &y_complexes);
            let z = endpoint_with_parts(157, &charts, &z_complexes);
            let strict = admit_strict(w, x, 48, ColorRank::Verified, ColorRank::Verified, cx);
            let f = admit_complex_refinement(
                x,
                y,
                x_complexes[0].id,
                y_complexes[0].id,
                x_complexes[0].resolution.id,
                y_complexes[0].resolution.id,
                40,
                ColorRank::Verified,
                ColorRank::Validated,
                cx,
            );
            let good = admit_complex_refinement(
                y,
                z,
                y_complexes[0].id,
                z_complexes[0].id,
                y_complexes[0].resolution.id,
                z_complexes[0].resolution.id,
                44,
                ColorRank::Validated,
                ColorRank::Estimated,
                cx,
            );
            let wrong = admit_complex_refinement(
                y,
                z,
                y_complexes[1].id,
                z_complexes[0].id,
                y_complexes[1].resolution.id,
                z_complexes[0].resolution.id,
                52,
                ColorRank::Validated,
                ColorRank::Estimated,
                cx,
            );
            assert_eq!(
                compose_derived_morphisms_v1(&f, &wrong, cx),
                Err(DerivedMorphismErrorV1::CompositionRefinementMismatch)
            );

            let mixed_prefix =
                compose_derived_morphisms_v1(&strict, &f, cx).expect("strict then refinement");
            assert_eq!(
                compose_derived_morphisms_v1(&mixed_prefix, &wrong, cx),
                Err(DerivedMorphismErrorV1::CompositionRefinementMismatch)
            );
            let left = compose_derived_morphisms_v1(&mixed_prefix, &good, cx)
                .expect("(strict-refinement)-refinement");
            let refinement_path =
                compose_derived_morphisms_v1(&f, &good, cx).expect("refinement path");
            let right = compose_derived_morphisms_v1(&strict, &refinement_path, cx)
                .expect("strict-(refinement-refinement)");
            assert_eq!(left, right);
            assert_eq!(
                left.class(),
                AdmittedDerivedMorphismClassV1::HeterogeneousPath
            );
            assert!(matches!(
                left.primitive_path(),
                [
                    AdmittedDerivedPrimitiveV1::Strict { .. },
                    AdmittedDerivedPrimitiveV1::DeclaredComplexRefinement(_),
                    AdmittedDerivedPrimitiveV1::DeclaredComplexRefinement(_)
                ]
            ));
        });
    }

    #[test]
    fn declared_inclusion_receipt_binds_nominal_artifacts_without_containment_authority() {
        with_cx(false, |cx| {
            let source = endpoint(180);
            let target = endpoint(181);
            let base_ir = inclusion_ir(
                source,
                target,
                20,
                ColorRank::Validated,
                ColorRank::Estimated,
            );
            let base = admit_between_endpoints(base_ir, source, target, cx)
                .expect("base inclusion declaration");
            let replay = admit_between_endpoints(base_ir, source, target, cx)
                .expect("replayed inclusion declaration");
            assert_eq!(base, replay);
            assert_eq!(
                base.primitive_path(),
                &[AdmittedDerivedPrimitiveV1::DeclaredInclusion(
                    DeclaredInclusionPrimitiveV1 {
                        source_geometry: source.id,
                        target_geometry: target.id,
                        map: DerivedInclusionMapIdV1::from_bytes([21; 32]),
                        containment: DerivedWitnessIdV1::from_bytes([22; 32]),
                    }
                )]
            );

            let mut strict_ir = strict_ir(
                source,
                target,
                21,
                ColorRank::Validated,
                ColorRank::Estimated,
            );
            strict_ir.equivalence = base_ir.equivalence;
            let strict = admit_between_endpoints(strict_ir, source, target, cx)
                .expect("strict comparator with the same no-claim artifact");
            assert_eq!(base.evidence(), strict.evidence());
            assert_eq!(base.no_equivalence_claims(), strict.no_equivalence_claims());
            assert_ne!(base.id(), strict.id(), "map families require distinct tags");

            let mut changed_map = base_ir;
            if let DerivedMorphismKindV1::DeclaredInclusion { map, .. } = &mut changed_map.kind {
                *map = DerivedInclusionMapIdV1::from_bytes([23; 32]);
            }
            let changed_map = admit_between_endpoints(changed_map, source, target, cx)
                .expect("changed inclusion map remains a declaration");
            assert_ne!(base.id(), changed_map.id());

            let mut changed_containment = base_ir;
            if let DerivedMorphismKindV1::DeclaredInclusion { containment, .. } =
                &mut changed_containment.kind
            {
                *containment = DerivedWitnessIdV1::from_bytes([24; 32]);
            }
            let changed_containment =
                admit_between_endpoints(changed_containment, source, target, cx)
                    .expect("changed containment artifact remains a declaration");
            assert_ne!(base.id(), changed_containment.id());
        });
    }

    #[test]
    #[allow(clippy::too_many_lines)] // Independent mutations for every inclusion authority boundary.
    fn declared_inclusion_refuses_missing_ids_and_authority_laundering() {
        with_cx(false, |cx| {
            let source = endpoint(182);
            let target = endpoint(183);
            let base = inclusion_ir(
                source,
                target,
                30,
                ColorRank::Validated,
                ColorRank::Estimated,
            );

            let mut zero_map = base;
            if let DerivedMorphismKindV1::DeclaredInclusion { map, .. } = &mut zero_map.kind {
                *map = DerivedInclusionMapIdV1::from_bytes([0; 32]);
            }
            assert_eq!(
                admit_between_endpoints(zero_map, source, target, cx),
                Err(DerivedMorphismErrorV1::MissingIdentity {
                    field: "inclusion-map"
                })
            );

            let mut zero_containment = base;
            if let DerivedMorphismKindV1::DeclaredInclusion { containment, .. } =
                &mut zero_containment.kind
            {
                *containment = DerivedWitnessIdV1::from_bytes([0; 32]);
            }
            assert_eq!(
                admit_between_endpoints(zero_containment, source, target, cx),
                Err(DerivedMorphismErrorV1::MissingIdentity {
                    field: "inclusion-containment"
                })
            );

            let mut zero_no_claim = base;
            if let DerivedEquivalenceBoundaryV1::NoClaim { artifact } =
                &mut zero_no_claim.equivalence
            {
                *artifact = DerivedNoClaimIdV1::from_bytes([0; 32]);
            }
            assert_eq!(
                admit_between_endpoints(zero_no_claim, source, target, cx),
                Err(DerivedMorphismErrorV1::MissingIdentity {
                    field: "no-equivalence-claim"
                })
            );

            let mut identity_evidence = base;
            identity_evidence.evidence = DerivedEvidenceTransportV1::Identity;
            assert_eq!(
                admit_between_endpoints(identity_evidence, source, target, cx),
                Err(DerivedMorphismErrorV1::EvidenceOrientationMismatch)
            );

            let mut self_identity_evidence = inclusion_ir(
                source,
                source,
                31,
                ColorRank::Validated,
                ColorRank::Estimated,
            );
            self_identity_evidence.evidence = DerivedEvidenceTransportV1::Identity;
            assert_eq!(
                admit_between_endpoints(self_identity_evidence, source, source, cx),
                Err(DerivedMorphismErrorV1::EvidenceOrientationMismatch)
            );

            let mut laundering = base;
            laundering.equivalence = DerivedEquivalenceBoundaryV1::IdentityOnly;
            assert_eq!(
                admit_between_endpoints(laundering, source, target, cx),
                Err(DerivedMorphismErrorV1::EquivalenceLaundering)
            );

            let incompatible_target = GeometryEndpointV1 {
                model_version: DerivedModelVersionIdV1::from_bytes([5; 32]),
                ..target
            };
            assert_eq!(
                admit_between_endpoints(base, source, incompatible_target, cx),
                Err(DerivedMorphismErrorV1::ModelVersionMismatch)
            );
        });
    }

    #[test]
    #[allow(clippy::too_many_lines)] // Homogeneous and mixed inclusion composition laws.
    fn declared_inclusions_compose_associatively_as_typed_paths() {
        with_cx(false, |cx| {
            let x = endpoint(184);
            let y = endpoint(185);
            let z = endpoint(186);
            let w = endpoint(187);
            let f = admit_inclusion(x, y, 40, ColorRank::Verified, ColorRank::Validated, cx);
            let g = admit_inclusion(y, z, 41, ColorRank::Validated, ColorRank::Estimated, cx);
            let h = admit_inclusion(z, w, 42, ColorRank::Estimated, ColorRank::Estimated, cx);
            let fg = compose_derived_morphisms_v1(&f, &g, cx).expect("f then g");
            let gh = compose_derived_morphisms_v1(&g, &h, cx).expect("g then h");
            let left = compose_derived_morphisms_v1(&fg, &h, cx).expect("(fg)h");
            let right = compose_derived_morphisms_v1(&f, &gh, cx).expect("f(gh)");
            assert_eq!(left, right);
            assert_eq!(
                left.class(),
                AdmittedDerivedMorphismClassV1::DeclaredInclusionPath
            );
            assert_eq!(
                left.primitive_path(),
                &[
                    AdmittedDerivedPrimitiveV1::DeclaredInclusion(DeclaredInclusionPrimitiveV1 {
                        source_geometry: x.id,
                        target_geometry: y.id,
                        map: DerivedInclusionMapIdV1::from_bytes([41; 32]),
                        containment: DerivedWitnessIdV1::from_bytes([42; 32]),
                    },),
                    AdmittedDerivedPrimitiveV1::DeclaredInclusion(DeclaredInclusionPrimitiveV1 {
                        source_geometry: y.id,
                        target_geometry: z.id,
                        map: DerivedInclusionMapIdV1::from_bytes([42; 32]),
                        containment: DerivedWitnessIdV1::from_bytes([43; 32]),
                    },),
                    AdmittedDerivedPrimitiveV1::DeclaredInclusion(DeclaredInclusionPrimitiveV1 {
                        source_geometry: z.id,
                        target_geometry: w.id,
                        map: DerivedInclusionMapIdV1::from_bytes([43; 32]),
                        containment: DerivedWitnessIdV1::from_bytes([44; 32]),
                    },),
                ]
            );
            assert_eq!(left.primitive_factors(), &[f.id(), g.id(), h.id()]);
            assert_eq!(
                left.no_equivalence_claims(),
                &[
                    DerivedNoClaimIdV1::from_bytes([104; 32]),
                    DerivedNoClaimIdV1::from_bytes([105; 32]),
                    DerivedNoClaimIdV1::from_bytes([106; 32]),
                ]
            );
            assert_eq!(
                compose_derived_morphisms_v1(&admit_identity(x, cx), &f, cx)
                    .expect("inclusion left identity"),
                f
            );
            assert_eq!(
                compose_derived_morphisms_v1(&h, &admit_identity(w, cx), cx)
                    .expect("inclusion right identity"),
                h
            );
            assert_eq!(
                compose_derived_morphisms_v1(&f, &h, cx),
                Err(DerivedMorphismErrorV1::CompositionEndpointMismatch)
            );

            let strict_xy = admit_strict(x, y, 50, ColorRank::Verified, ColorRank::Validated, cx);
            let inclusion_yz =
                admit_inclusion(y, z, 51, ColorRank::Validated, ColorRank::Estimated, cx);
            let strict_zw = admit_strict(z, w, 52, ColorRank::Estimated, ColorRank::Estimated, cx);
            let strict_inclusion = compose_derived_morphisms_v1(&strict_xy, &inclusion_yz, cx)
                .expect("strict then inclusion");
            let inclusion_strict = compose_derived_morphisms_v1(&inclusion_yz, &strict_zw, cx)
                .expect("inclusion then strict");
            let mixed_left = compose_derived_morphisms_v1(&strict_inclusion, &strict_zw, cx)
                .expect("(strict-inclusion)-strict");
            let mixed_right = compose_derived_morphisms_v1(&strict_xy, &inclusion_strict, cx)
                .expect("strict-(inclusion-strict)");
            assert_eq!(mixed_left, mixed_right);
            assert_eq!(
                mixed_left.class(),
                AdmittedDerivedMorphismClassV1::HeterogeneousPath
            );
            assert!(matches!(
                mixed_left.primitive_path(),
                [
                    AdmittedDerivedPrimitiveV1::Strict { .. },
                    AdmittedDerivedPrimitiveV1::DeclaredInclusion(_),
                    AdmittedDerivedPrimitiveV1::Strict { .. }
                ]
            ));
            assert_eq!(
                mixed_left.primitive_factors(),
                &[strict_xy.id(), inclusion_yz.id(), strict_zw.id()]
            );
            assert_eq!(
                mixed_left.no_equivalence_claims(),
                &[
                    DerivedNoClaimIdV1::from_bytes([114; 32]),
                    DerivedNoClaimIdV1::from_bytes([115; 32]),
                    DerivedNoClaimIdV1::from_bytes([116; 32]),
                ]
            );
        });
    }

    #[test]
    #[allow(clippy::too_many_lines)] // Replay, exact accessors, and ordered-leg identity.
    fn standalone_span_binds_ordered_admitted_legs_without_direct_transport() {
        with_cx(false, |cx| {
            assert_eq!(
                <DerivedSpanCorrespondenceIdentitySchemaV1 as CanonicalSchema>::CONTEXT,
                "exact source, common apex, exact target, ordered admitted legs, and no-authority boundary"
            );
            let source = endpoint(190);
            let apex = endpoint(191);
            let target = endpoint(192);
            let left = admit_strict(
                apex,
                source,
                60,
                ColorRank::Verified,
                ColorRank::Validated,
                cx,
            );
            let right = admit_strict(
                apex,
                target,
                61,
                ColorRank::Verified,
                ColorRank::Estimated,
                cx,
            );
            let ir = span_ir(source, apex, target, &left, &right, 62);
            let admitted = admit_derived_span_correspondence_v1(ir, &left, &right, cx)
                .expect("valid standalone span");
            let replay = admit_derived_span_correspondence_v1(ir, &left, &right, cx)
                .expect("replayed standalone span");
            assert_eq!(admitted, replay);
            assert_eq!(admitted.source(), source.id);
            assert_eq!(admitted.apex(), apex.id);
            assert_eq!(admitted.target(), target.id);
            assert_eq!(admitted.left_leg(), left.id());
            assert_eq!(admitted.right_leg(), right.id());
            assert_eq!(
                admitted.no_claim(),
                DerivedNoClaimIdV1::from_bytes([62; 32])
            );
            assert_eq!(admitted.id(), admitted.identity_receipt().id());

            let mut changed_source = ir;
            changed_source.source = geometry_id(200);
            let mut changed_apex = ir;
            changed_apex.apex = geometry_id(200);
            let mut changed_target = ir;
            changed_target.target = geometry_id(200);
            let mut changed_left_leg = ir;
            changed_left_leg.left_leg = right.id();
            let mut changed_right_leg = ir;
            changed_right_leg.right_leg = left.id();
            for (field, changed) in [
                ("source", changed_source),
                ("apex", changed_apex),
                ("target", changed_target),
                ("left-leg", changed_left_leg),
                ("right-leg", changed_right_leg),
            ] {
                let changed = span_correspondence_receipt(changed, cx)
                    .expect("one-field span receipt mutation");
                assert_ne!(admitted.id(), changed.id(), "{field} must move identity");
            }

            let mut changed_no_claim = ir;
            changed_no_claim.no_claim = DerivedNoClaimIdV1::from_bytes([63; 32]);
            let changed_no_claim =
                admit_derived_span_correspondence_v1(changed_no_claim, &left, &right, cx)
                    .expect("changed no-claim remains a structural span");
            assert_ne!(admitted.id(), changed_no_claim.id());

            let common_target = endpoint(193);
            let first_leg = admit_strict(
                apex,
                common_target,
                64,
                ColorRank::Verified,
                ColorRank::Validated,
                cx,
            );
            let second_leg = admit_strict(
                apex,
                common_target,
                65,
                ColorRank::Verified,
                ColorRank::Estimated,
                cx,
            );
            let ordered_ir = span_ir(
                common_target,
                apex,
                common_target,
                &first_leg,
                &second_leg,
                66,
            );
            let ordered =
                admit_derived_span_correspondence_v1(ordered_ir, &first_leg, &second_leg, cx)
                    .expect("ordered equal-endpoint span");
            let swapped_ir = span_ir(
                common_target,
                apex,
                common_target,
                &second_leg,
                &first_leg,
                66,
            );
            let swapped =
                admit_derived_span_correspondence_v1(swapped_ir, &second_leg, &first_leg, cx)
                    .expect("swapped equal-endpoint span");
            assert_ne!(
                ordered.id(),
                swapped.id(),
                "left/right leg order is semantic"
            );
        });
    }

    #[test]
    #[allow(clippy::too_many_lines)] // Independent raw-binding and four orientation refusals.
    fn standalone_span_refuses_raw_leg_and_orientation_mismatches() {
        with_cx(false, |cx| {
            let source = endpoint(194);
            let apex = endpoint(195);
            let target = endpoint(196);
            let other = endpoint(197);
            let left = admit_strict(
                apex,
                source,
                70,
                ColorRank::Verified,
                ColorRank::Validated,
                cx,
            );
            let right = admit_strict(
                apex,
                target,
                71,
                ColorRank::Verified,
                ColorRank::Estimated,
                cx,
            );
            let base = span_ir(source, apex, target, &left, &right, 72);

            let mut wrong_schema = base;
            wrong_schema.schema_version += 1;
            assert_eq!(
                admit_derived_span_correspondence_v1(wrong_schema, &left, &right, cx),
                Err(DerivedSpanCorrespondenceErrorV1::UnsupportedSchemaVersion {
                    found: 2,
                    supported: DERIVED_SPAN_CORRESPONDENCE_SCHEMA_VERSION_V1,
                })
            );

            let mut zero_no_claim = base;
            zero_no_claim.no_claim = DerivedNoClaimIdV1::from_bytes([0; 32]);
            assert_eq!(
                admit_derived_span_correspondence_v1(zero_no_claim, &left, &right, cx),
                Err(DerivedSpanCorrespondenceErrorV1::MissingIdentity {
                    field: "no-correspondence-claim"
                })
            );

            let mut wrong_left_id = base;
            wrong_left_id.left_leg = right.id();
            assert_eq!(
                admit_derived_span_correspondence_v1(wrong_left_id, &left, &right, cx),
                Err(DerivedSpanCorrespondenceErrorV1::LegIdentityMismatch { field: "left-leg" })
            );

            let mut wrong_right_id = base;
            wrong_right_id.right_leg = left.id();
            assert_eq!(
                admit_derived_span_correspondence_v1(wrong_right_id, &left, &right, cx),
                Err(DerivedSpanCorrespondenceErrorV1::LegIdentityMismatch { field: "right-leg" })
            );

            let wrong_left_source = admit_strict(
                other,
                source,
                73,
                ColorRank::Verified,
                ColorRank::Validated,
                cx,
            );
            let wrong_left_source_ir =
                span_ir(source, apex, target, &wrong_left_source, &right, 74);
            assert_eq!(
                admit_derived_span_correspondence_v1(
                    wrong_left_source_ir,
                    &wrong_left_source,
                    &right,
                    cx,
                ),
                Err(DerivedSpanCorrespondenceErrorV1::LegOrientationMismatch {
                    field: "left-source-apex"
                })
            );

            let wrong_left_target = admit_strict(
                apex,
                other,
                75,
                ColorRank::Verified,
                ColorRank::Validated,
                cx,
            );
            let wrong_left_target_ir =
                span_ir(source, apex, target, &wrong_left_target, &right, 76);
            assert_eq!(
                admit_derived_span_correspondence_v1(
                    wrong_left_target_ir,
                    &wrong_left_target,
                    &right,
                    cx,
                ),
                Err(DerivedSpanCorrespondenceErrorV1::LegOrientationMismatch {
                    field: "left-target-source"
                })
            );

            let wrong_right_source = admit_strict(
                other,
                target,
                77,
                ColorRank::Verified,
                ColorRank::Estimated,
                cx,
            );
            let wrong_right_source_ir =
                span_ir(source, apex, target, &left, &wrong_right_source, 78);
            assert_eq!(
                admit_derived_span_correspondence_v1(
                    wrong_right_source_ir,
                    &left,
                    &wrong_right_source,
                    cx,
                ),
                Err(DerivedSpanCorrespondenceErrorV1::LegOrientationMismatch {
                    field: "right-source-apex"
                })
            );

            let wrong_right_target = admit_strict(
                apex,
                other,
                79,
                ColorRank::Verified,
                ColorRank::Estimated,
                cx,
            );
            let wrong_right_target_ir =
                span_ir(source, apex, target, &left, &wrong_right_target, 80);
            assert_eq!(
                admit_derived_span_correspondence_v1(
                    wrong_right_target_ir,
                    &left,
                    &wrong_right_target,
                    cx,
                ),
                Err(DerivedSpanCorrespondenceErrorV1::LegOrientationMismatch {
                    field: "right-target-target"
                })
            );
        });
    }

    #[test]
    fn standalone_span_admits_graph_shape_and_refuses_entry_cancellation() {
        with_cx(false, |cx| {
            let source_and_apex = endpoint(198);
            let target = endpoint(199);
            let left_identity = admit_identity(source_and_apex, cx);
            let right = admit_strict(
                source_and_apex,
                target,
                81,
                ColorRank::Verified,
                ColorRank::Validated,
                cx,
            );
            let ir = span_ir(
                source_and_apex,
                source_and_apex,
                target,
                &left_identity,
                &right,
                82,
            );
            let graph = admit_derived_span_correspondence_v1(ir, &left_identity, &right, cx)
                .expect("identity-left graph shape is structurally admissible");
            assert_eq!(graph.source(), graph.apex());
            assert_eq!(graph.left_leg(), left_identity.id());
            assert_eq!(graph.right_leg(), right.id());

            with_cx(true, |cancelled_cx| {
                assert_eq!(
                    admit_derived_span_correspondence_v1(ir, &left_identity, &right, cancelled_cx,),
                    Err(DerivedSpanCorrespondenceErrorV1::Cancelled {
                        stage: "span-admission-entry"
                    })
                );
            });
        });
    }

    #[test]
    fn declared_chart_map_receipt_binds_family_specific_ids() {
        with_cx(false, |cx| {
            let x_charts = [chart(110, 2, 2, 10, 1.0), chart(109, 2, 2, 10, 1.0)];
            let y_charts = [chart(111, 2, 2, 10, 1.0), chart(117, 2, 2, 10, 1.0)];
            let x = endpoint_with_charts(112, &x_charts);
            let y = endpoint_with_charts(113, &y_charts);
            let base_ir = chart_map_ir(
                x,
                y,
                x_charts[0].id,
                y_charts[0].id,
                114,
                ColorRank::Validated,
                ColorRank::Estimated,
            );
            let base = admit_between_endpoints(base_ir, x, y, cx).expect("base chart map");
            let replay = admit_between_endpoints(base_ir, x, y, cx).expect("replayed chart map");
            assert_eq!(base, replay);
            assert_eq!(
                base.declared_chart_maps(),
                &[DeclaredChartMapPrimitiveV1 {
                    source_geometry: x.id,
                    target_geometry: y.id,
                    source_chart: x_charts[0].id,
                    target_chart: y_charts[0].id,
                    overlap: DerivedChartOverlapIdV1::from_bytes([115; 32]),
                    map: DerivedChartMapIdV1::from_bytes([116; 32]),
                }]
            );

            let mut changed_source = base_ir;
            if let DerivedMorphismKindV1::DeclaredChartMap { source_chart, .. } =
                &mut changed_source.kind
            {
                *source_chart = x_charts[1].id;
            }
            let changed_source = admit_between_endpoints(changed_source, x, y, cx)
                .expect("changed source chart remains structural");
            assert_ne!(base.id(), changed_source.id());

            let mut changed_target = base_ir;
            if let DerivedMorphismKindV1::DeclaredChartMap { target_chart, .. } =
                &mut changed_target.kind
            {
                *target_chart = y_charts[1].id;
            }
            let changed_target = admit_between_endpoints(changed_target, x, y, cx)
                .expect("changed target chart remains structural");
            assert_ne!(base.id(), changed_target.id());

            let mut changed_overlap = base_ir;
            if let DerivedMorphismKindV1::DeclaredChartMap { overlap, .. } =
                &mut changed_overlap.kind
            {
                *overlap = DerivedChartOverlapIdV1::from_bytes([118; 32]);
            }
            let changed_overlap = admit_between_endpoints(changed_overlap, x, y, cx)
                .expect("changed overlap remains structural");
            assert_ne!(base.id(), changed_overlap.id());

            let mut changed_map = base_ir;
            if let DerivedMorphismKindV1::DeclaredChartMap { map, .. } = &mut changed_map.kind {
                *map = DerivedChartMapIdV1::from_bytes([119; 32]);
            }
            let changed_map = admit_between_endpoints(changed_map, x, y, cx)
                .expect("changed map remains structural");
            assert_ne!(base.id(), changed_map.id());

            let generic_strict =
                admit_strict(x, y, 114, ColorRank::Validated, ColorRank::Estimated, cx);
            assert_ne!(base.id(), generic_strict.id());
        });
    }

    #[test]
    #[allow(clippy::too_many_lines)] // One mutation table for every chart-map admission boundary.
    fn declared_chart_map_admission_refuses_missing_or_incompatible_charts() {
        with_cx(false, |cx| {
            let x_charts = [chart(120, 2, 2, 11, 1.0)];
            let y_charts = [chart(121, 2, 2, 11, 1.0)];
            let x = endpoint_with_charts(122, &x_charts);
            let y = endpoint_with_charts(123, &y_charts);
            let base = chart_map_ir(
                x,
                y,
                x_charts[0].id,
                y_charts[0].id,
                124,
                ColorRank::Validated,
                ColorRank::Estimated,
            );

            let mut missing = base;
            if let DerivedMorphismKindV1::DeclaredChartMap { source_chart, .. } = &mut missing.kind
            {
                *source_chart = chart_id(125);
            }
            assert_eq!(
                admit_between_endpoints(missing, x, y, cx),
                Err(DerivedMorphismErrorV1::MissingChart {
                    field: "source-chart"
                })
            );

            let mut missing_target = base;
            if let DerivedMorphismKindV1::DeclaredChartMap { target_chart, .. } =
                &mut missing_target.kind
            {
                *target_chart = chart_id(125);
            }
            assert_eq!(
                admit_between_endpoints(missing_target, x, y, cx),
                Err(DerivedMorphismErrorV1::MissingChart {
                    field: "target-chart"
                })
            );

            let mut zero_source = base;
            if let DerivedMorphismKindV1::DeclaredChartMap { source_chart, .. } =
                &mut zero_source.kind
            {
                *source_chart = chart_id(0);
            }
            assert_eq!(
                admit_between_endpoints(zero_source, x, y, cx),
                Err(DerivedMorphismErrorV1::MissingIdentity {
                    field: "source-chart"
                })
            );

            let mut zero_target = base;
            if let DerivedMorphismKindV1::DeclaredChartMap { target_chart, .. } =
                &mut zero_target.kind
            {
                *target_chart = chart_id(0);
            }
            assert_eq!(
                admit_between_endpoints(zero_target, x, y, cx),
                Err(DerivedMorphismErrorV1::MissingIdentity {
                    field: "target-chart"
                })
            );

            let mut zero_overlap = base;
            if let DerivedMorphismKindV1::DeclaredChartMap { overlap, .. } = &mut zero_overlap.kind
            {
                *overlap = DerivedChartOverlapIdV1::from_bytes([0; 32]);
            }
            assert_eq!(
                admit_between_endpoints(zero_overlap, x, y, cx),
                Err(DerivedMorphismErrorV1::MissingIdentity {
                    field: "chart-overlap"
                })
            );

            let mut zero_map = base;
            if let DerivedMorphismKindV1::DeclaredChartMap { map, .. } = &mut zero_map.kind {
                *map = DerivedChartMapIdV1::from_bytes([0; 32]);
            }
            assert_eq!(
                admit_between_endpoints(zero_map, x, y, cx),
                Err(DerivedMorphismErrorV1::MissingIdentity { field: "chart-map" })
            );

            let mut laundering = base;
            laundering.equivalence = DerivedEquivalenceBoundaryV1::IdentityOnly;
            assert_eq!(
                admit_between_endpoints(laundering, x, y, cx),
                Err(DerivedMorphismErrorV1::EquivalenceLaundering)
            );

            let dimension_chart = [chart(126, 3, 3, 11, 1.0)];
            let dimension_target = endpoint_with_charts(123, &dimension_chart);
            assert_eq!(
                admit_between_endpoints(
                    chart_map_ir(
                        x,
                        dimension_target,
                        x_charts[0].id,
                        dimension_chart[0].id,
                        127,
                        ColorRank::Validated,
                        ColorRank::Estimated,
                    ),
                    x,
                    dimension_target,
                    cx,
                ),
                Err(DerivedMorphismErrorV1::ChartDimensionMismatch)
            );

            let quantity_chart = [chart(128, 2, 2, 12, 1.0)];
            let quantity_target = endpoint_with_charts(123, &quantity_chart);
            assert_eq!(
                admit_between_endpoints(
                    chart_map_ir(
                        x,
                        quantity_target,
                        x_charts[0].id,
                        quantity_chart[0].id,
                        129,
                        ColorRank::Validated,
                        ColorRank::Estimated,
                    ),
                    x,
                    quantity_target,
                    cx,
                ),
                Err(DerivedMorphismErrorV1::ChartConventionMismatch)
            );
        });
    }

    #[test]
    fn declared_chart_map_admission_checks_each_chart_convention_field() {
        with_cx(false, |cx| {
            let source_charts = [chart(150, 2, 2, 15, 1.0)];
            let source = endpoint_with_charts(151, &source_charts);

            let mut wrong_frame = chart(152, 2, 2, 15, 1.0);
            wrong_frame.frame = DerivedFrameIdV1::from_bytes([99; 32]);
            let wrong_frame_charts = [wrong_frame];
            let frame_target = endpoint_with_charts(153, &wrong_frame_charts);
            assert_eq!(
                admit_between_endpoints(
                    chart_map_ir(
                        source,
                        frame_target,
                        source_charts[0].id,
                        wrong_frame_charts[0].id,
                        154,
                        ColorRank::Validated,
                        ColorRank::Estimated,
                    ),
                    source,
                    frame_target,
                    cx,
                ),
                Err(DerivedMorphismErrorV1::ChartConventionMismatch)
            );

            let mut wrong_system = chart(155, 2, 2, 15, 1.0);
            wrong_system.coordinates.system = DerivedUnitSystemIdV1::from_bytes([98; 32]);
            let wrong_system_charts = [wrong_system];
            let system_target = endpoint_with_charts(153, &wrong_system_charts);
            assert_eq!(
                admit_between_endpoints(
                    chart_map_ir(
                        source,
                        system_target,
                        source_charts[0].id,
                        wrong_system_charts[0].id,
                        158,
                        ColorRank::Validated,
                        ColorRank::Estimated,
                    ),
                    source,
                    system_target,
                    cx,
                ),
                Err(DerivedMorphismErrorV1::ChartConventionMismatch)
            );

            let wrong_scale_charts = [chart(159, 2, 2, 15, f64::from_bits(1.0_f64.to_bits() + 1))];
            let scale_target = endpoint_with_charts(153, &wrong_scale_charts);
            assert_eq!(
                admit_between_endpoints(
                    chart_map_ir(
                        source,
                        scale_target,
                        source_charts[0].id,
                        wrong_scale_charts[0].id,
                        162,
                        ColorRank::Validated,
                        ColorRank::Estimated,
                    ),
                    source,
                    scale_target,
                    cx,
                ),
                Err(DerivedMorphismErrorV1::ChartConventionMismatch)
            );
        });
    }

    #[test]
    fn declared_chart_map_self_map_cannot_use_identity_evidence_transport() {
        with_cx(false, |cx| {
            let charts = [chart(146, 2, 2, 14, 1.0), chart(147, 2, 2, 14, 1.0)];
            let object = endpoint_with_charts(148, &charts);
            let mut ir = chart_map_ir(
                object,
                object,
                charts[0].id,
                charts[1].id,
                149,
                ColorRank::Validated,
                ColorRank::Validated,
            );
            ir.evidence = DerivedEvidenceTransportV1::Identity;
            assert_eq!(
                admit_between_endpoints(ir, object, object, cx),
                Err(DerivedMorphismErrorV1::EvidenceOrientationMismatch)
            );
        });
    }

    #[test]
    fn chart_composition_refuses_wrong_chart_but_accepts_a_strict_boundary() {
        with_cx(false, |cx| {
            let x_charts = [chart(130, 2, 2, 13, 1.0)];
            let y_charts = [chart(131, 2, 2, 13, 1.0), chart(132, 2, 2, 13, 1.0)];
            let z_charts = [chart(133, 2, 2, 13, 1.0)];
            let x = endpoint_with_charts(134, &x_charts);
            let y = endpoint_with_charts(135, &y_charts);
            let z = endpoint_with_charts(136, &z_charts);
            let f = admit_chart_map(
                x,
                y,
                x_charts[0].id,
                y_charts[0].id,
                137,
                ColorRank::Validated,
                ColorRank::Validated,
                cx,
            );
            let wrong_middle = admit_chart_map(
                y,
                z,
                y_charts[1].id,
                z_charts[0].id,
                140,
                ColorRank::Validated,
                ColorRank::Estimated,
                cx,
            );
            assert_eq!(
                compose_derived_morphisms_v1(&f, &wrong_middle, cx),
                Err(DerivedMorphismErrorV1::CompositionChartMismatch)
            );

            let strict = admit_strict(x, y, 143, ColorRank::Validated, ColorRank::Validated, cx);
            let mixed = compose_derived_morphisms_v1(&strict, &wrong_middle, cx)
                .expect("strict boundary carries no chart seam");
            assert_eq!(
                mixed.class(),
                AdmittedDerivedMorphismClassV1::HeterogeneousPath
            );
        });
    }

    #[test]
    fn composition_refuses_nominal_evidence_seam_and_variance_mismatch() {
        with_cx(false, |cx| {
            let x = endpoint(80);
            let y = endpoint(81);
            let z = endpoint(82);
            let first = admit_strict(x, y, 83, ColorRank::Validated, ColorRank::Validated, cx);

            let mut mismatched_ir = strict_ir(y, z, 84, ColorRank::Validated, ColorRank::Estimated);
            mismatched_ir.evidence = DerivedEvidenceTransportV1::BalanceCorestrictionCovariant {
                input_geometry: y.id,
                output_geometry: z.id,
                input_evidence: DerivedEvidenceArtifactIdV1::from_bytes([99; 32]),
                output_evidence: evidence_id(z.id),
                input_rank: ColorRank::Validated,
                output_rank: ColorRank::Estimated,
            };
            let mismatched = admit_between_endpoints(mismatched_ir, y, z, cx)
                .expect("nominal evidence IDs are structural declarations");
            assert_eq!(
                compose_derived_morphisms_v1(&first, &mismatched, cx),
                Err(DerivedMorphismErrorV1::CompositionEvidenceMismatch)
            );

            let rank_mismatched =
                admit_strict(y, z, 86, ColorRank::Verified, ColorRank::Estimated, cx);
            assert_eq!(
                compose_derived_morphisms_v1(&first, &rank_mismatched, cx),
                Err(DerivedMorphismErrorV1::CompositionEvidenceMismatch)
            );

            let restriction = admit_between_endpoints(
                restriction_ir(y, z, 85, ColorRank::Verified, ColorRank::Validated),
                y,
                z,
                cx,
            )
            .expect("valid restriction");
            assert_eq!(
                compose_derived_morphisms_v1(&first, &restriction, cx),
                Err(DerivedMorphismErrorV1::CompositionVarianceMismatch)
            );
        });
    }

    #[test]
    fn strict_maps_cannot_launder_equivalence_or_conventions() {
        with_cx(false, |cx| {
            let x = endpoint(50);
            let mut y = endpoint(51);
            let mut equivalence = strict_ir(x, y, 52, ColorRank::Validated, ColorRank::Validated);
            equivalence.equivalence = DerivedEquivalenceBoundaryV1::IdentityOnly;
            assert_eq!(
                admit_between_endpoints(equivalence, x, y, cx),
                Err(DerivedMorphismErrorV1::EquivalenceLaundering)
            );

            let mut version_changed = y;
            version_changed.model_version = DerivedModelVersionIdV1::from_bytes([98; 32]);
            assert_eq!(
                admit_between_endpoints(
                    strict_ir(
                        x,
                        version_changed,
                        53,
                        ColorRank::Validated,
                        ColorRank::Validated,
                    ),
                    x,
                    version_changed,
                    cx,
                ),
                Err(DerivedMorphismErrorV1::ModelVersionMismatch)
            );

            y.frame = DerivedFrameIdV1::from_bytes([99; 32]);
            assert_eq!(
                admit_between_endpoints(
                    strict_ir(x, y, 54, ColorRank::Validated, ColorRank::Validated,),
                    x,
                    y,
                    cx,
                ),
                Err(DerivedMorphismErrorV1::FrameMismatch)
            );
        });
    }

    #[test]
    fn replay_is_deterministic_and_entry_cancellation_fails_closed() {
        with_cx(false, |cx| {
            let x = endpoint(60);
            let y = endpoint(61);
            let first = admit_strict(x, y, 62, ColorRank::Verified, ColorRank::Validated, cx);
            let replay = admit_strict(x, y, 62, ColorRank::Verified, ColorRank::Validated, cx);
            assert_eq!(first, replay);
        });
        with_cx(true, |cx| {
            let x = endpoint(60);
            let y = endpoint(61);
            assert!(matches!(
                admit_between_endpoints(
                    strict_ir(x, y, 62, ColorRank::Verified, ColorRank::Validated,),
                    x,
                    y,
                    cx,
                ),
                Err(DerivedMorphismErrorV1::Cancelled { .. })
            ));
        });
    }

    #[test]
    fn lineage_cap_and_composition_entry_cancellation_fail_closed() {
        assert_eq!(
            checked_combined_len("primitive-lineage", DERIVED_MORPHISM_MAX_FACTORS_V1 - 1, 1,),
            Ok(DERIVED_MORPHISM_MAX_FACTORS_V1)
        );
        assert_eq!(
            checked_combined_len("primitive-lineage", DERIVED_MORPHISM_MAX_FACTORS_V1, 1,),
            Err(DerivedMorphismErrorV1::ResourceLimit {
                field: "primitive-lineage",
                requested: DERIVED_MORPHISM_MAX_FACTORS_V1 + 1,
                limit: DERIVED_MORPHISM_MAX_FACTORS_V1,
            })
        );

        let (identity, strict) = with_cx(false, |cx| {
            let x = endpoint(70);
            let y = endpoint(71);
            (
                admit_identity(x, cx),
                admit_strict(x, y, 72, ColorRank::Validated, ColorRank::Estimated, cx),
            )
        });
        with_cx(true, |cx| {
            assert!(matches!(
                compose_derived_morphisms_v1(&identity, &strict, cx),
                Err(DerivedMorphismErrorV1::Cancelled { .. })
            ));
        });
    }

    #[test]
    fn fixed_resolution_candidate_is_domain_separate_replayable_and_no_authority() {
        assert_ne!(
            <DerivedFixedResolutionQuasiIsomorphismCandidateIdentitySchemaV1 as CanonicalSchema>::DOMAIN,
            <DerivedMorphismIdentitySchemaV1 as CanonicalSchema>::DOMAIN
        );
        assert_ne!(
            <DerivedFixedResolutionQuasiIsomorphismCandidateIdentitySchemaV1 as CanonicalSchema>::DOMAIN,
            <DerivedSpanCorrespondenceIdentitySchemaV1 as CanonicalSchema>::DOMAIN
        );
        assert_eq!(
            <DerivedFixedResolutionQuasiIsomorphismCandidateIdentitySchemaV1 as CanonicalSchema>::FIELDS
                .len(),
            16
        );

        with_cx(false, |cx| {
            let (source, target, path, ir) = fixed_resolution_candidate_fixture(cx);
            let first = admit_derived_fixed_resolution_quasi_isomorphism_candidate_v1(
                &ir, &source, &target, &path, cx,
            )
            .expect("valid structural quasi-isomorphism candidate");
            let replay = admit_derived_fixed_resolution_quasi_isomorphism_candidate_v1(
                &ir, &source, &target, &path, cx,
            )
            .expect("deterministic structural candidate replay");

            assert_eq!(first, replay);
            assert_eq!(first.refinement_path(), path.id());
            assert_eq!(first.source_geometry(), source.id());
            assert_eq!(first.target_geometry(), target.id());
            assert_eq!(first.source_local_model(), ir.source_local_model);
            assert_eq!(first.target_local_model(), ir.target_local_model);
            assert_eq!(first.complex_role(), DerivedComplexRoleV1::Tangent);
            assert_eq!(first.source_complex(), complex_id(70));
            assert_eq!(first.target_complex(), complex_id(73));
            assert_eq!(first.source_resolution(), resolution_id(80));
            assert_eq!(first.target_resolution(), resolution_id(83));
            assert_eq!(
                first.source_scope_witness(),
                DerivedWitnessIdV1::from_bytes([112; 32])
            );
            assert_eq!(
                first.target_scope_witness(),
                DerivedWitnessIdV1::from_bytes([115; 32])
            );
            assert_eq!(first.nominal_theorem(), ir.nominal_theorem);
            assert_eq!(first.nominal_checker(), ir.nominal_checker);
            assert_eq!(first.nominal_check_receipt(), ir.nominal_check_receipt);
            assert_eq!(first.no_authority(), ir.no_authority);
            assert_eq!(
                path.class(),
                AdmittedDerivedMorphismClassV1::DeclaredComplexRefinementPath
            );
            assert_eq!(path.no_equivalence_claims().len(), 1);
        });
    }

    #[test]
    fn fixed_resolution_candidate_receipt_binds_every_ordered_field() {
        with_cx(false, |cx| {
            let (source, target, path, ir) = fixed_resolution_candidate_fixture(cx);
            let binding = fixed_resolution_quasi_isomorphism_candidate_binding(
                &ir, &source, &target, &path, cx,
            )
            .expect("valid candidate binding");
            let baseline = fixed_resolution_quasi_isomorphism_candidate_receipt(&binding, cx)
                .expect("candidate receipt")
                .id();

            macro_rules! assert_field_moves_identity {
                ($field:ident, $value:expr) => {{
                    let mut changed = binding;
                    changed.$field = $value;
                    let changed =
                        fixed_resolution_quasi_isomorphism_candidate_receipt(&changed, cx)
                            .expect("mutated candidate receipt")
                            .id();
                    assert_ne!(baseline, changed, stringify!($field));
                }};
            }

            assert_field_moves_identity!(
                refinement_path,
                DerivedMorphismIdV1::parse_slice(&[201; 32]).expect("nonzero path identity")
            );
            assert_field_moves_identity!(source_geometry, geometry_id(202));
            assert_field_moves_identity!(target_geometry, geometry_id(203));
            assert_field_moves_identity!(
                source_local_model,
                DerivedLocalModelIdV1::from_bytes([204; 32])
            );
            assert_field_moves_identity!(
                target_local_model,
                DerivedLocalModelIdV1::from_bytes([205; 32])
            );
            assert_field_moves_identity!(complex_role, DerivedComplexRoleV1::Cotangent);
            assert_field_moves_identity!(source_complex, complex_id(206));
            assert_field_moves_identity!(target_complex, complex_id(207));
            assert_field_moves_identity!(source_resolution, resolution_id(208));
            assert_field_moves_identity!(target_resolution, resolution_id(209));
            assert_field_moves_identity!(
                source_scope_witness,
                DerivedWitnessIdV1::from_bytes([210; 32])
            );
            assert_field_moves_identity!(
                target_scope_witness,
                DerivedWitnessIdV1::from_bytes([211; 32])
            );
            assert_field_moves_identity!(
                nominal_theorem,
                DerivedTheoremIdV1::from_bytes([212; 32])
            );
            assert_field_moves_identity!(
                nominal_checker,
                DerivedCheckerIdV1::from_bytes([213; 32])
            );
            assert_field_moves_identity!(
                nominal_check_receipt,
                DerivedWitnessIdV1::from_bytes([214; 32])
            );
            assert_field_moves_identity!(no_authority, DerivedNoClaimIdV1::from_bytes([215; 32]));
        });
    }

    #[test]
    fn fixed_resolution_candidate_refuses_path_role_scope_and_opaque_identity_defects() {
        with_cx(false, |cx| {
            let (source, target, path, ir) = fixed_resolution_candidate_fixture(cx);

            let mut bad_schema = ir;
            bad_schema.schema_version = 2;
            assert!(matches!(
                admit_derived_fixed_resolution_quasi_isomorphism_candidate_v1(
                    &bad_schema,
                    &source,
                    &target,
                    &path,
                    cx,
                ),
                Err(DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1::UnsupportedSchemaVersion {
                    found: 2,
                    supported: DERIVED_FIXED_RESOLUTION_QUASI_ISOMORPHISM_CANDIDATE_SCHEMA_VERSION_V1,
                })
            ));

            let mut wrong_path = ir;
            wrong_path.refinement_path =
                DerivedMorphismIdV1::parse_slice(&[216; 32]).expect("nonzero path identity");
            assert_eq!(
                admit_derived_fixed_resolution_quasi_isomorphism_candidate_v1(
                    &wrong_path,
                    &source,
                    &target,
                    &path,
                    cx,
                ),
                Err(DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1::PathIdentityMismatch)
            );

            let strict_path = admit_derived_morphism_v1(
                DerivedMorphismIrV1 {
                    schema_version: DERIVED_MORPHISM_SCHEMA_VERSION_V1,
                    source: source.id(),
                    target: target.id(),
                    kind: DerivedMorphismKindV1::Strict {
                        witness: DerivedWitnessIdV1::from_bytes([220; 32]),
                    },
                    evidence: DerivedEvidenceTransportV1::BalanceCorestrictionCovariant {
                        input_geometry: source.id(),
                        output_geometry: target.id(),
                        input_evidence: evidence_id(source.id()),
                        output_evidence: evidence_id(target.id()),
                        input_rank: ColorRank::Validated,
                        output_rank: ColorRank::Estimated,
                    },
                    equivalence: DerivedEquivalenceBoundaryV1::NoClaim {
                        artifact: DerivedNoClaimIdV1::from_bytes([221; 32]),
                    },
                },
                &source,
                &target,
                cx,
            )
            .expect("valid non-refinement control path");
            let mut wrong_class = ir;
            wrong_class.refinement_path = strict_path.id();
            assert_eq!(
                admit_derived_fixed_resolution_quasi_isomorphism_candidate_v1(
                    &wrong_class,
                    &source,
                    &target,
                    &strict_path,
                    cx,
                ),
                Err(
                    DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1::PathClassMismatch {
                        found: AdmittedDerivedMorphismClassV1::Strict,
                    }
                )
            );

            let mut wrong_selector = ir;
            wrong_selector.source_resolution = resolution_id(222);
            assert_eq!(
                admit_derived_fixed_resolution_quasi_isomorphism_candidate_v1(
                    &wrong_selector,
                    &source,
                    &target,
                    &path,
                    cx,
                ),
                Err(
                    DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1::PathShapeMismatch {
                        field: "source-selector",
                        index: 0,
                    }
                )
            );

            let mut wrong_role = ir;
            wrong_role.complex_role = DerivedComplexRoleV1::Cotangent;
            assert!(matches!(
                admit_derived_fixed_resolution_quasi_isomorphism_candidate_v1(
                    &wrong_role,
                    &source,
                    &target,
                    &path,
                    cx,
                ),
                Err(
                    DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1::ComplexRoleMismatch {
                        field: "source-complex-role",
                    }
                )
            ));

            let mut missing_model = ir;
            missing_model.source_local_model = DerivedLocalModelIdV1::from_bytes([217; 32]);
            assert_eq!(
                admit_derived_fixed_resolution_quasi_isomorphism_candidate_v1(
                    &missing_model,
                    &source,
                    &target,
                    &path,
                    cx,
                ),
                Err(
                    DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1::MissingLocalModel {
                        field: "source-local-model",
                    }
                )
            );

            for (field, mut malformed) in [
                ("nominal-theorem", ir),
                ("nominal-checker", ir),
                ("nominal-check-receipt", ir),
                ("no-authority", ir),
            ] {
                match field {
                    "nominal-theorem" => {
                        malformed.nominal_theorem = DerivedTheoremIdV1::from_bytes([0; 32]);
                    }
                    "nominal-checker" => {
                        malformed.nominal_checker = DerivedCheckerIdV1::from_bytes([0; 32]);
                    }
                    "nominal-check-receipt" => {
                        malformed.nominal_check_receipt = DerivedWitnessIdV1::from_bytes([0; 32]);
                    }
                    "no-authority" => {
                        malformed.no_authority = DerivedNoClaimIdV1::from_bytes([0; 32]);
                    }
                    _ => unreachable!(),
                }
                assert_eq!(
                    admit_derived_fixed_resolution_quasi_isomorphism_candidate_v1(
                        &malformed, &source, &target, &path, cx,
                    ),
                    Err(
                        DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1::MissingIdentity {
                            field,
                        }
                    )
                );
            }

            let literal_model = DerivedLocalModelV1 {
                presentation: PresentationScopeV1::Literal {
                    no_claim: DerivedNoClaimIdV1::from_bytes([218; 32]),
                },
                ..source.ir().local_models[0].clone()
            };
            assert_eq!(
                fixed_resolution_scope_witness(
                    &literal_model,
                    ir.source_resolution,
                    "source-presentation",
                    "source-presentation-resolution",
                    "source-scope-witness",
                ),
                Err(DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1::PresentationScopeMismatch {
                    field: "source-presentation",
                })
            );
            let external_model = DerivedLocalModelV1 {
                presentation: PresentationScopeV1::ExternallyChecked {
                    witness: DerivedWitnessIdV1::from_bytes([219; 32]),
                },
                ..source.ir().local_models[0].clone()
            };
            assert_eq!(
                fixed_resolution_scope_witness(
                    &external_model,
                    ir.source_resolution,
                    "source-presentation",
                    "source-presentation-resolution",
                    "source-scope-witness",
                ),
                Err(DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1::PresentationScopeMismatch {
                    field: "source-presentation",
                })
            );
        });
    }

    #[test]
    fn fixed_resolution_candidate_entry_cancellation_fails_closed() {
        let (source, target, path, ir) = with_cx(false, fixed_resolution_candidate_fixture);
        with_cx(true, |cx| {
            assert!(matches!(
                admit_derived_fixed_resolution_quasi_isomorphism_candidate_v1(
                    &ir, &source, &target, &path, cx,
                ),
                Err(
                    DerivedFixedResolutionQuasiIsomorphismCandidateErrorV1::Cancelled {
                        stage: "candidate-admission-entry",
                    }
                )
            ));
        });
    }
}
