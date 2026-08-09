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
pub mod wire;

/// The current `.fsim` schema version. Readers admit exactly this version;
/// older envelopes must pass through [`migration::migrate_envelope`].
///
/// Version 2 adds the optional `(fan-system ...)` cooling subsection (bead
/// frn2i.1); version-1 documents are valid version-2 documents unchanged,
/// and the migration is the receipted envelope rewrite.
pub const FSIM_VERSION: u32 = 2;

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
    THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY, resolve_bindings,
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
    AirflowLeakage, Budgets, ConsequenceClass, Cooling, DecisionGate, DefaultReceipt, EntityDecl,
    Envelope, Fan, FanCurveDecl, FanCurvePoint, FanToleranceBasis, GeometryArtifact,
    GeometryAssignment, InterfaceCardBinding, InterfaceState, MaterialBinding, Metadata,
    OutputRequest, PerfectContactBinding, PowerDissipation, ProjectSpec, RequirementDirection,
    RequirementSeverity, RequirementSource, RequirementSourceKind, RequirementSourceReview,
    SafetyFactorPolicy, Seeds, SolverSettings, ThermalLimit, UnitsDoctrine, Vent, Versions,
    requirement_source_reviews,
};
pub use wire::{
    CanonicalizationReceipt, DecodedProject, ProjectError, canonical_hash, lower, parse_json,
    parse_sexpr, parse_sexpr_lenient, print_json, print_sexpr, recognize,
};
