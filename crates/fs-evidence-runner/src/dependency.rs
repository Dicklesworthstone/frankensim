//! Frozen owner-scoped direct-dependency evolution policy.
//!
//! This module is data and pure validation only. It does not resolve Cargo
//! metadata, fetch a registry, move a constellation pin, or add a dependency
//! to the package manifest.

use crate::construction::{ConstructionErrorKindV2, ConstructionErrorV2};
use crate::value::StableTokenV2;
use fs_blake3::{ContentHash, hash_domain};
use std::collections::BTreeSet;

/// Domain for the exact current direct-dependency declaration root.
pub const CURRENT_DIRECT_DEPENDENCY_DECLARATION_DOMAIN_V1: &str =
    "org.frankensim.fs-evidence-runner.current-direct-dependency-declaration.v1";

/// Closed source route admitted by the eventual package dependency policy.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum DependencySourceRouteV1 {
    /// A path-owned FrankenSim workspace crate.
    WorkspacePath = 1,
    /// A source and revision pinned by the constellation lock.
    PinnedConstellationSibling = 2,
    /// An ambient package registry; represented only so validation can refuse
    /// it deterministically.
    AmbientRegistry = 3,
}

/// Exact declaration-time identity of one dependency source.
///
/// Pinned rows mirror the phase-owned constellation identity tuple, but this
/// data does not inspect or attest the live checkout or `constellation.lock`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum DependencySourceIdentityV1 {
    /// A workspace package at one exact path relative to this package.
    WorkspacePath {
        /// Exact package-relative Cargo path.
        manifest_relative_path: &'static str,
    },
    /// A sibling package identified by an exact constellation tuple.
    PinnedConstellationSibling {
        /// Exact package-relative Cargo path.
        manifest_relative_path: &'static str,
        /// Exact constellation library key.
        lock_library: &'static str,
        /// Exact constellation version identity.
        lock_version: &'static str,
        /// Exact lowercase hexadecimal Git revision.
        git_head: &'static str,
    },
    /// An ambient registry route, represented only for deterministic refusal.
    AmbientRegistry,
}

impl DependencySourceIdentityV1 {
    /// Route class implied by this exact source identity.
    #[must_use]
    pub const fn route(self) -> DependencySourceRouteV1 {
        match self {
            Self::WorkspacePath { .. } => DependencySourceRouteV1::WorkspacePath,
            Self::PinnedConstellationSibling { .. } => {
                DependencySourceRouteV1::PinnedConstellationSibling
            }
            Self::AmbientRegistry => DependencySourceRouteV1::AmbientRegistry,
        }
    }

    /// Exact Cargo path relative to this package, when the source is a path.
    #[must_use]
    pub const fn manifest_relative_path(self) -> Option<&'static str> {
        match self {
            Self::WorkspacePath {
                manifest_relative_path,
            }
            | Self::PinnedConstellationSibling {
                manifest_relative_path,
                ..
            } => Some(manifest_relative_path),
            Self::AmbientRegistry => None,
        }
    }
}

/// Sole phase allowed to introduce one direct normal dependency.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum DependencyOwnerPhaseV1 {
    /// Phase-1 pure base schema.
    BaseSchema = 1,
    /// Later capability, cancellation, process, and persistence phase.
    CapabilityCancellationProcessPersistence = 2,
}

/// Cargo-route qualifiers forbidden from the direct normal dependency set.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum DependencyRouteQualifierV1 {
    /// Cargo `optional = true`.
    Optional = 1,
    /// A renamed dependency package.
    Renamed = 2,
    /// A build-dependency route.
    BuildDependency = 3,
    /// A proc-macro dependency route.
    ProcMacro = 4,
    /// A target-specific dependency route.
    TargetSpecific = 5,
}

/// One immutable package-policy row.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct DependencyPolicyRowV1 {
    package: &'static str,
    source_identity: DependencySourceIdentityV1,
    required_features: &'static [&'static str],
    owner: DependencyOwnerPhaseV1,
}

impl DependencyPolicyRowV1 {
    const fn new(
        package: &'static str,
        source_identity: DependencySourceIdentityV1,
        required_features: &'static [&'static str],
        owner: DependencyOwnerPhaseV1,
    ) -> Self {
        Self {
            package,
            source_identity,
            required_features,
            owner,
        }
    }

    /// Exact Cargo package name.
    #[must_use]
    pub const fn package(self) -> &'static str {
        self.package
    }

    /// Required source class.
    #[must_use]
    pub const fn source(self) -> DependencySourceRouteV1 {
        self.source_identity.route()
    }

    /// Exact package-relative path and, for siblings, constellation identity.
    #[must_use]
    pub const fn source_identity(self) -> DependencySourceIdentityV1 {
        self.source_identity
    }

    /// Exact sorted feature set.
    #[must_use]
    pub const fn required_features(self) -> &'static [&'static str] {
        self.required_features
    }

    /// Sole phase owner.
    #[must_use]
    pub const fn owner(self) -> DependencyOwnerPhaseV1 {
        self.owner
    }
}

/// Exact current Phase-1 direct normal-dependency set.
pub const CURRENT_DIRECT_NORMAL_DEPENDENCIES_V1: [DependencyPolicyRowV1; 1] =
    [DependencyPolicyRowV1::new(
        "fs-blake3",
        DependencySourceIdentityV1::WorkspacePath {
            manifest_relative_path: "../fs-blake3",
        },
        &[],
        DependencyOwnerPhaseV1::BaseSchema,
    )];

/// Exact eventual complete-package direct normal-dependency allowlist.
///
/// The latter two rows remain data until their named later phase changes the
/// manifest and lockfile. Their presence here does not make them current.
pub const EVENTUAL_DIRECT_NORMAL_DEPENDENCY_ALLOWLIST_V1: [DependencyPolicyRowV1; 3] = [
    DependencyPolicyRowV1::new(
        "fs-blake3",
        DependencySourceIdentityV1::WorkspacePath {
            manifest_relative_path: "../fs-blake3",
        },
        &[],
        DependencyOwnerPhaseV1::BaseSchema,
    ),
    DependencyPolicyRowV1::new(
        "asupersync",
        DependencySourceIdentityV1::PinnedConstellationSibling {
            manifest_relative_path: "../../../asupersync",
            lock_library: "asupersync",
            lock_version: "0.3.8",
            git_head: "5973b0ff31f405ae90fa9e6e2ef5f61a75c5b78b",
        },
        &[],
        DependencyOwnerPhaseV1::CapabilityCancellationProcessPersistence,
    ),
    DependencyPolicyRowV1::new(
        "fsqlite",
        DependencySourceIdentityV1::PinnedConstellationSibling {
            manifest_relative_path: "../../../frankensqlite/crates/fsqlite",
            lock_library: "frankensqlite",
            lock_version: "unversioned-workspace",
            git_head: "987cfb97f86d3fca4d9b44e7871f427636b10126",
        },
        &["async-api"],
        DependencyOwnerPhaseV1::CapabilityCancellationProcessPersistence,
    ),
];

/// Canonically identify the exact current declaration-time dependency rows.
///
/// This root binds static policy data only. It does not inspect Cargo
/// metadata, the lockfile, the constellation, or a sibling checkout, and
/// therefore is not live supply-chain proof.
#[must_use]
pub fn current_direct_dependency_declaration_root_v1() -> ContentHash {
    dependency_declaration_root(&CURRENT_DIRECT_NORMAL_DEPENDENCIES_V1)
}

fn dependency_declaration_root(rows: &[DependencyPolicyRowV1]) -> ContentHash {
    let mut bytes = Vec::new();
    bytes.extend_from_slice(b"FSCURRENTDIRECTDEPENDENCY\x01");
    bytes.extend_from_slice(
        &u32::try_from(rows.len())
            .expect("the static dependency declaration count fits u32")
            .to_be_bytes(),
    );
    for row in rows {
        push_static_str(&mut bytes, row.package());
        bytes.extend_from_slice(&(row.source() as u16).to_be_bytes());
        match row.source_identity() {
            DependencySourceIdentityV1::WorkspacePath {
                manifest_relative_path,
            } => {
                bytes.extend_from_slice(&1_u16.to_be_bytes());
                push_static_str(&mut bytes, manifest_relative_path);
            }
            DependencySourceIdentityV1::PinnedConstellationSibling {
                manifest_relative_path,
                lock_library,
                lock_version,
                git_head,
            } => {
                bytes.extend_from_slice(&2_u16.to_be_bytes());
                push_static_str(&mut bytes, manifest_relative_path);
                push_static_str(&mut bytes, lock_library);
                push_static_str(&mut bytes, lock_version);
                push_static_str(&mut bytes, git_head);
            }
            DependencySourceIdentityV1::AmbientRegistry => {
                bytes.extend_from_slice(&3_u16.to_be_bytes());
            }
        }
        bytes.extend_from_slice(
            &u32::try_from(row.required_features().len())
                .expect("the static feature count fits u32")
                .to_be_bytes(),
        );
        for feature in row.required_features() {
            push_static_str(&mut bytes, feature);
        }
        bytes.extend_from_slice(&(row.owner() as u16).to_be_bytes());
    }
    hash_domain(CURRENT_DIRECT_DEPENDENCY_DECLARATION_DOMAIN_V1, &bytes)
}

fn push_static_str(bytes: &mut Vec<u8>, value: &'static str) {
    bytes.extend_from_slice(
        &u32::try_from(value.len())
            .expect("the static dependency string length fits u32")
            .to_be_bytes(),
    );
    bytes.extend_from_slice(value.as_bytes());
}

/// Untrusted description of one direct Cargo route.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PresentedDependencyRouteV1 {
    package: StableTokenV2,
    source_identity: DependencySourceIdentityV1,
    features: Box<[StableTokenV2]>,
    owner: DependencyOwnerPhaseV1,
    qualifiers: Box<[DependencyRouteQualifierV1]>,
}

impl PresentedDependencyRouteV1 {
    /// Assemble an untrusted route while validating token syntax and
    /// canonicalizing its duplicate-free feature set.
    pub fn presented(
        package: impl Into<String>,
        source_identity: DependencySourceIdentityV1,
        features: Vec<String>,
        owner: DependencyOwnerPhaseV1,
        mut qualifiers: Vec<DependencyRouteQualifierV1>,
    ) -> Result<Self, ConstructionErrorV2> {
        let package = StableTokenV2::new(package.into()).map_err(|error| {
            ConstructionErrorV2::new(
                ConstructionErrorKindV2::Incompatible,
                "dependency.package",
                "one validated stable package token",
                format_args!("{error:?}"),
            )
        })?;
        let mut features = features
            .into_iter()
            .map(|feature| {
                StableTokenV2::new(feature).map_err(|error| {
                    ConstructionErrorV2::new(
                        ConstructionErrorKindV2::Incompatible,
                        "dependency.feature",
                        "one validated stable feature token",
                        format_args!("{error:?}"),
                    )
                })
            })
            .collect::<Result<Vec<_>, _>>()?;
        let mut seen = BTreeSet::new();
        for feature in &features {
            if !seen.insert(feature.as_str()) {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Duplicate,
                    "dependency.features",
                    "one duplicate-free feature set",
                    feature.as_str(),
                ));
            }
        }
        features.sort();
        let mut seen_qualifiers = BTreeSet::new();
        for qualifier in &qualifiers {
            if !seen_qualifiers.insert(*qualifier) {
                return Err(ConstructionErrorV2::new(
                    ConstructionErrorKindV2::Duplicate,
                    "dependency.qualifiers",
                    "one duplicate-free route-qualifier set",
                    *qualifier as u16,
                ));
            }
        }
        qualifiers.sort_unstable();
        Ok(Self {
            package,
            source_identity,
            features: features.into_boxed_slice(),
            owner,
            qualifiers: qualifiers.into_boxed_slice(),
        })
    }

    /// Presented package name.
    #[must_use]
    pub const fn package(&self) -> &StableTokenV2 {
        &self.package
    }

    /// Presented exact source identity.
    #[must_use]
    pub const fn source_identity(&self) -> DependencySourceIdentityV1 {
        self.source_identity
    }

    /// Presented feature set in canonical order.
    #[must_use]
    pub fn features(&self) -> &[StableTokenV2] {
        &self.features
    }

    /// Presented forbidden route qualifiers in canonical order.
    #[must_use]
    pub fn qualifiers(&self) -> &[DependencyRouteQualifierV1] {
        &self.qualifiers
    }
}

/// Validate the one exact current Phase-1 dependency set.
pub fn validate_current_direct_dependencies_v1(
    presented: &[PresentedDependencyRouteV1],
) -> Result<(), ConstructionErrorV2> {
    validate_exact_dependency_set(presented, &CURRENT_DIRECT_NORMAL_DEPENDENCIES_V1)
}

/// Validate an eventual manifest proposal against the closed complete-package
/// allowlist and exact owner/feature/source rows.
pub fn validate_eventual_direct_dependencies_v1(
    presented: &[PresentedDependencyRouteV1],
) -> Result<(), ConstructionErrorV2> {
    validate_exact_dependency_set(presented, &EVENTUAL_DIRECT_NORMAL_DEPENDENCY_ALLOWLIST_V1)
}

fn validate_exact_dependency_set(
    presented: &[PresentedDependencyRouteV1],
    expected: &[DependencyPolicyRowV1],
) -> Result<(), ConstructionErrorV2> {
    let mut seen = BTreeSet::new();
    for route in presented {
        if !seen.insert(route.package.as_str()) {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Duplicate,
                "dependency.package",
                "one exact route per package",
                route.package.as_str(),
            ));
        }
    }
    if presented.len() < expected.len() {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Missing,
            "dependency.routes",
            "the complete exact dependency set",
            presented.len(),
        ));
    }
    if presented.len() > expected.len() {
        return Err(ConstructionErrorV2::new(
            ConstructionErrorKindV2::Unexpected,
            "dependency.routes",
            "no route beyond the complete exact dependency set",
            presented.len(),
        ));
    }

    for (ordinal, (route, policy)) in presented.iter().zip(expected).enumerate() {
        let actual_features = route
            .features
            .iter()
            .map(StableTokenV2::as_str)
            .collect::<Vec<_>>();
        let exact = route.package.as_str() == policy.package
            && route.source_identity == policy.source_identity
            && actual_features == policy.required_features
            && route.owner == policy.owner
            && route.qualifiers.is_empty()
            && route.source_identity.route() != DependencySourceRouteV1::AmbientRegistry;
        if !exact {
            return Err(ConstructionErrorV2::new(
                ConstructionErrorKindV2::Incompatible,
                "dependency.route",
                "the exact ordered package/source/features/owner normal route",
                format_args!("{ordinal}:{}", route.package.as_str()),
            ));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn route(policy: DependencyPolicyRowV1) -> PresentedDependencyRouteV1 {
        PresentedDependencyRouteV1::presented(
            policy.package(),
            policy.source_identity(),
            policy
                .required_features()
                .iter()
                .map(|feature| (*feature).to_owned())
                .collect(),
            policy.owner(),
            Vec::new(),
        )
        .expect("policy fixture")
    }

    #[test]
    fn current_and_eventual_literal_oracles_are_exact() {
        assert_eq!(CURRENT_DIRECT_NORMAL_DEPENDENCIES_V1.len(), 1);
        assert_eq!(
            CURRENT_DIRECT_NORMAL_DEPENDENCIES_V1[0],
            DependencyPolicyRowV1::new(
                "fs-blake3",
                DependencySourceIdentityV1::WorkspacePath {
                    manifest_relative_path: "../fs-blake3",
                },
                &[],
                DependencyOwnerPhaseV1::BaseSchema,
            )
        );
        assert_eq!(
            CURRENT_DIRECT_NORMAL_DEPENDENCIES_V1[0]
                .source_identity()
                .manifest_relative_path(),
            Some("../fs-blake3")
        );
        assert_eq!(
            EVENTUAL_DIRECT_NORMAL_DEPENDENCY_ALLOWLIST_V1,
            [
                DependencyPolicyRowV1::new(
                    "fs-blake3",
                    DependencySourceIdentityV1::WorkspacePath {
                        manifest_relative_path: "../fs-blake3",
                    },
                    &[],
                    DependencyOwnerPhaseV1::BaseSchema,
                ),
                DependencyPolicyRowV1::new(
                    "asupersync",
                    DependencySourceIdentityV1::PinnedConstellationSibling {
                        manifest_relative_path: "../../../asupersync",
                        lock_library: "asupersync",
                        lock_version: "0.3.8",
                        git_head: "5973b0ff31f405ae90fa9e6e2ef5f61a75c5b78b",
                    },
                    &[],
                    DependencyOwnerPhaseV1::CapabilityCancellationProcessPersistence,
                ),
                DependencyPolicyRowV1::new(
                    "fsqlite",
                    DependencySourceIdentityV1::PinnedConstellationSibling {
                        manifest_relative_path: "../../../frankensqlite/crates/fsqlite",
                        lock_library: "frankensqlite",
                        lock_version: "unversioned-workspace",
                        git_head: "987cfb97f86d3fca4d9b44e7871f427636b10126",
                    },
                    &["async-api"],
                    DependencyOwnerPhaseV1::CapabilityCancellationProcessPersistence,
                ),
            ]
        );
        assert!(
            validate_current_direct_dependencies_v1(&[route(
                CURRENT_DIRECT_NORMAL_DEPENDENCIES_V1[0]
            )])
            .is_ok()
        );
        assert!(
            validate_eventual_direct_dependencies_v1(
                &EVENTUAL_DIRECT_NORMAL_DEPENDENCY_ALLOWLIST_V1.map(route)
            )
            .is_ok()
        );
        let declaration_root = current_direct_dependency_declaration_root_v1();
        assert_eq!(
            declaration_root,
            current_direct_dependency_declaration_root_v1()
        );
        for mutant in [
            DependencyPolicyRowV1::new(
                "fs-blake4",
                CURRENT_DIRECT_NORMAL_DEPENDENCIES_V1[0].source_identity(),
                &[],
                DependencyOwnerPhaseV1::BaseSchema,
            ),
            DependencyPolicyRowV1::new(
                "fs-blake3",
                DependencySourceIdentityV1::WorkspacePath {
                    manifest_relative_path: "../other",
                },
                &[],
                DependencyOwnerPhaseV1::BaseSchema,
            ),
            DependencyPolicyRowV1::new(
                "fs-blake3",
                CURRENT_DIRECT_NORMAL_DEPENDENCIES_V1[0].source_identity(),
                &["unexpected-feature"],
                DependencyOwnerPhaseV1::BaseSchema,
            ),
            DependencyPolicyRowV1::new(
                "fs-blake3",
                CURRENT_DIRECT_NORMAL_DEPENDENCIES_V1[0].source_identity(),
                &[],
                DependencyOwnerPhaseV1::CapabilityCancellationProcessPersistence,
            ),
        ] {
            assert_ne!(dependency_declaration_root(&[mutant]), declaration_root);
        }
    }

    #[test]
    fn every_absent_extra_or_order_mutant_refuses() {
        let exact = EVENTUAL_DIRECT_NORMAL_DEPENDENCY_ALLOWLIST_V1.map(route);
        assert!(validate_eventual_direct_dependencies_v1(&exact[..2]).is_err());
        let mut extra = exact.to_vec();
        extra.push(
            PresentedDependencyRouteV1::presented(
                "other",
                DependencySourceIdentityV1::WorkspacePath {
                    manifest_relative_path: "../other",
                },
                Vec::new(),
                DependencyOwnerPhaseV1::BaseSchema,
                Vec::new(),
            )
            .expect("extra route"),
        );
        assert!(validate_eventual_direct_dependencies_v1(&extra).is_err());
        let mut reordered = exact.to_vec();
        reordered.swap(0, 1);
        assert!(validate_eventual_direct_dependencies_v1(&reordered).is_err());
        let mut duplicate = exact.to_vec();
        duplicate[1] = duplicate[0].clone();
        assert!(validate_eventual_direct_dependencies_v1(&duplicate).is_err());
    }

    #[test]
    fn every_forbidden_route_or_owner_feature_source_mutant_refuses() {
        let expected = CURRENT_DIRECT_NORMAL_DEPENDENCIES_V1[0];
        for qualifier in [
            DependencyRouteQualifierV1::Optional,
            DependencyRouteQualifierV1::Renamed,
            DependencyRouteQualifierV1::BuildDependency,
            DependencyRouteQualifierV1::ProcMacro,
            DependencyRouteQualifierV1::TargetSpecific,
        ] {
            let mutant = PresentedDependencyRouteV1::presented(
                expected.package(),
                expected.source_identity(),
                Vec::new(),
                expected.owner(),
                vec![qualifier],
            )
            .expect("flag mutant");
            assert!(validate_current_direct_dependencies_v1(&[mutant]).is_err());
        }
        for mutant in [
            PresentedDependencyRouteV1::presented(
                expected.package(),
                DependencySourceIdentityV1::AmbientRegistry,
                Vec::new(),
                expected.owner(),
                Vec::new(),
            )
            .expect("registry mutant"),
            PresentedDependencyRouteV1::presented(
                expected.package(),
                DependencySourceIdentityV1::WorkspacePath {
                    manifest_relative_path: "../wrong-package",
                },
                Vec::new(),
                expected.owner(),
                Vec::new(),
            )
            .expect("workspace-path mutant"),
            PresentedDependencyRouteV1::presented(
                expected.package(),
                expected.source_identity(),
                vec!["unexpected".to_owned()],
                expected.owner(),
                Vec::new(),
            )
            .expect("feature mutant"),
            PresentedDependencyRouteV1::presented(
                expected.package(),
                expected.source_identity(),
                Vec::new(),
                DependencyOwnerPhaseV1::CapabilityCancellationProcessPersistence,
                Vec::new(),
            )
            .expect("owner mutant"),
        ] {
            assert!(validate_current_direct_dependencies_v1(&[mutant]).is_err());
        }
    }

    #[test]
    fn every_pinned_source_identity_component_is_exact() {
        let exact = EVENTUAL_DIRECT_NORMAL_DEPENDENCY_ALLOWLIST_V1.map(route);
        let policy = EVENTUAL_DIRECT_NORMAL_DEPENDENCY_ALLOWLIST_V1[1];
        for source_identity in [
            DependencySourceIdentityV1::PinnedConstellationSibling {
                manifest_relative_path: "../../../wrong",
                lock_library: "asupersync",
                lock_version: "0.3.8",
                git_head: "5973b0ff31f405ae90fa9e6e2ef5f61a75c5b78b",
            },
            DependencySourceIdentityV1::PinnedConstellationSibling {
                manifest_relative_path: "../../../asupersync",
                lock_library: "wrong",
                lock_version: "0.3.8",
                git_head: "5973b0ff31f405ae90fa9e6e2ef5f61a75c5b78b",
            },
            DependencySourceIdentityV1::PinnedConstellationSibling {
                manifest_relative_path: "../../../asupersync",
                lock_library: "asupersync",
                lock_version: "0.3.9",
                git_head: "5973b0ff31f405ae90fa9e6e2ef5f61a75c5b78b",
            },
            DependencySourceIdentityV1::PinnedConstellationSibling {
                manifest_relative_path: "../../../asupersync",
                lock_library: "asupersync",
                lock_version: "0.3.8",
                git_head: "6973b0ff31f405ae90fa9e6e2ef5f61a75c5b78b",
            },
        ] {
            let mut mutant = exact.clone();
            mutant[1] = PresentedDependencyRouteV1::presented(
                policy.package(),
                source_identity,
                Vec::new(),
                policy.owner(),
                Vec::new(),
            )
            .expect("source identity mutant");
            let error = validate_eventual_direct_dependencies_v1(&mutant)
                .expect_err("every identity mutation must refuse");
            assert_eq!(error.kind(), ConstructionErrorKindV2::Incompatible);
            assert_eq!(error.field(), "dependency.route");
        }
    }

    #[test]
    fn compiled_phase_one_manifest_has_only_the_owned_normal_row() {
        let manifest = include_str!("../Cargo.toml");
        let dependency_section = manifest
            .split("[dependencies]\n")
            .nth(1)
            .expect("dependencies section")
            .split("\n[")
            .next()
            .expect("section body");
        let rows = dependency_section
            .lines()
            .map(str::trim)
            .filter(|line| !line.is_empty() && !line.starts_with('#'))
            .collect::<Vec<_>>();
        assert_eq!(rows, [r#"fs-blake3 = { path = "../fs-blake3" }"#]);
        assert!(!manifest.contains("[build-dependencies]"));
        assert!(!manifest.contains("[dev-dependencies]"));
        assert!(!manifest.contains("[target."));
        assert!(!manifest.contains("optional = true"));
        assert!(!manifest.contains("package = "));
    }

    #[test]
    fn current_package_is_a_root_member_and_imports_only_its_owned_dependency() {
        let root_manifest = include_str!("../../../Cargo.toml");
        assert!(
            root_manifest
                .lines()
                .any(|line| line.trim() == r#""crates/fs-evidence-runner","#)
        );

        let source_units = [
            ("budget", include_str!("budget.rs")),
            ("canonical", include_str!("canonical.rs")),
            ("capability", include_str!("capability.rs")),
            ("catalog", include_str!("catalog.rs")),
            ("command", include_str!("command.rs")),
            ("construction", include_str!("construction.rs")),
            ("coverage", include_str!("coverage.rs")),
            ("dependency", include_str!("dependency.rs")),
            ("diagnostic", include_str!("diagnostic.rs")),
            ("identity", include_str!("identity.rs")),
            ("lib", include_str!("lib.rs")),
            ("limits", include_str!("limits.rs")),
            ("logging", include_str!("logging.rs")),
            ("path", include_str!("path.rs")),
            ("projection", include_str!("projection.rs")),
            ("publication", include_str!("publication.rs")),
            ("state", include_str!("state.rs")),
            ("value", include_str!("value.rs")),
        ];
        let blake3_import_prefix = ["use ", "fs_", "blake3", "::"].concat();
        let importing_modules = source_units
            .iter()
            .filter_map(|(module, source)| {
                source
                    .lines()
                    .any(|line| line.trim_start().starts_with(&blake3_import_prefix))
                    .then_some(*module)
            })
            .collect::<std::collections::BTreeSet<_>>();
        assert_eq!(
            importing_modules,
            [
                "budget",
                "canonical",
                "catalog",
                "coverage",
                "dependency",
                "diagnostic",
                "identity",
                "limits",
                "logging",
                "projection",
                "publication",
            ]
            .into_iter()
            .collect(),
            "the exact handwritten modules that import the sole runtime dependency must remain closed"
        );
        let source = source_units
            .iter()
            .map(|(_, source)| *source)
            .collect::<String>();
        for forbidden in [
            ["use ", "asupersync"].concat(),
            ["asupersync", "::"].concat(),
            ["use ", "fsqlite"].concat(),
            ["fsqlite", "::"].concat(),
        ] {
            assert!(!source.contains(&forbidden), "{forbidden}");
        }
        assert!(!source.contains(&["extern", " crate "].concat()));
    }
}
