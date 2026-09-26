//! The versioned `.fsim` project schema (bead f85xj.6.1): the user-facing
//! contract for the ratified thermal-design-assurance vertical
//! (`frankensim-vertical-ratification-v1`).
//!
//! One semantic model, two spellings: a canonical s-expression grammar and an
//! isomorphic JSON rendering, both over `fs-ir`'s typed AST. Canonical bytes
//! are the checked s-expression render, hashed under
//! [`wire::FSIM_CANONICAL_DOMAIN`]; the JSON spelling reaches the same AST
//! and therefore the same canonical hash. The Five Explicits — units, seeds,
//! budgets, versions, capabilities — are mandatory sections; every omission
//! is a named [`fs_scenario::Violation`] with a fix, unknown fields are
//! refused, and the only defaults are receipted, never silent. Version bumps
//! travel through explicit [`migration`] receipts.

pub mod assignment;
pub mod bind;
pub mod decision;
pub mod fansystem;
pub mod migration;
pub mod spec;
pub mod study;
pub mod wire;

/// The current `.fsim` study schema version.
pub const STUDY_FSIM_VERSION: u32 = 1;

/// The current `.fsim` project schema version. Readers admit exactly this version;
/// older envelopes must pass through [`migration::migrate_envelope`].
///
/// Version 5 adds optional ambient radiation on convective exterior surfaces
/// (bead q61wp.74), with an immutable emissivity card, explicit reservoir and
/// query temperatures, and declared coupling controls. Version-4 documents
/// carry no radiation declaration and migrate with a receipted envelope and
/// schema rewrite; no radiative exchange is inferred. Version 4 added
/// airflow-convection, version 3 conduction, and version 2 fan-system inputs.
/// Version 6 adds an optional declared relative conductivity tolerance on a
/// material binding, with a mandatory basis and source (bead q61wp.72); the
/// solve stage propagates it into the Parameters budget term. Version-5
/// documents declare none and migrate without inventing one.
/// Version 7 adds an optional declared surface-offset band per geometry
/// artifact, with a mandatory basis and source (bead q61wp.79); the solve stage
/// propagates it into the Geometry budget term. Version-6 documents declare
/// none and migrate without inventing one.
pub const FSIM_VERSION: u32 = 7;

pub use assignment::{
    ConductionInterfaceLimits, ConductionInterfaceResolution, ConductionSourceFace,
    GEOMETRY_ASSIGNMENT_REPORT_DOMAIN, GEOMETRY_SOURCE_IDENTITY_DOMAIN, GeometryResolution,
    ImportedMeshLibrary, ResolvedConductionInterfacePair, ResolvedGeometryArtifact,
    ResolvedProjectAssignment, geometry_source_identity, resolve_conduction_interface_pairs,
    resolve_geometry_assignments,
};
pub use bind::{
    Advisory, BindingRequirements, BindingTarget, CONTACT_RESISTANCE_DIMS,
    CONTACT_RESISTANCE_PROPERTY, CardLibrary, MaterialResolution, RequiredProperty,
    ResolvedBinding, ResolvedProperty, RetainedReceipt, TEMPERATURE_AXIS,
    THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY, lower_card_backed_thermal_interfaces,
    resolve_bindings,
};
pub use decision::{
    PROJECT_DECISION_CONTEXT_IDENTITY_DOMAIN, PROJECT_REQUIREMENT_IDENTITY_DOMAIN,
    PROJECT_SAFETY_FACTOR_IDENTITY_DOMAIN, ProjectDecisionAuthority, ProjectDecisionContext,
    ProjectDecisionError, project_decision_authorities, project_decision_authority,
};
pub use fs_io::{HalfSpaceSide, MeshSelector};
pub use migration::{
    MigratedOrNative, MigratedProject, MigrationRule, ProjectMigrationReceipt, migrate_envelope,
    parse_sexpr_migrating,
};
pub use spec::{
    AirflowLeakage, Budgets, ConductionRadiation, ConductionRegion, ConductionSetup,
    ConsequenceClass, Cooling, DecisionGate, DefaultReceipt, EntityDecl, Envelope, Fan,
    FanCurveDecl, FanCurvePoint, FanToleranceBasis, GeometryArtifact, GeometryAssignment,
    InterfaceCardBinding, InterfaceState, MaterialBinding, MaterialTolerance, Metadata, SurfaceOffset,
    OutputRequest, PerfectContactBinding, PowerDissipation, ProjectSpec, RadiatingSurface,
    RequirementDirection, RequirementSeverity, RequirementSource, RequirementSourceKind,
    RequirementSourceReview, SafetyFactorPolicy, Seeds, SolverSettings, ThermalBoundary,
    ThermalBoundaryCondition, ThermalLimit, UnitsDoctrine, Vent, Versions,
    requirement_source_reviews,
};
pub use study::{
    StudyBudgets, StudyConstraints, StudyDomain, StudyHole, StudyObjective, StudyOptimizer,
    StudyPhysics, StudyScenario, StudySpec, canonical_study_hash, parse_study_json,
    parse_study_sexpr, print_study_json, print_study_sexpr,
};
pub use wire::{
    CanonicalizationReceipt, DecodedProject, ProjectError, canonical_hash, lower, parse_json,
    parse_sexpr, parse_sexpr_lenient, print_json, print_sexpr, recognize,
};
