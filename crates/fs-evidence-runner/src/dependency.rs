//! Frozen owner-scoped direct-dependency evolution policy.
//!
//! This module is data and pure validation only. It does not resolve Cargo
//! metadata, fetch a registry, move a constellation pin, or add a dependency
//! to the package manifest.

use crate::construction::{ConstructionErrorKindV2, ConstructionErrorV2};
use crate::value::StableTokenV2;
use std::collections::BTreeSet;

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

/// Sole phase allowed to introduce one direct normal dependency.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u16)]
pub enum DependencyOwnerPhaseV1 {
    /// Phase-1 pure base schema.
    BaseSchema = 1,
    /// Later capability, cancellation, process, and persistence phase.
    CapabilityCancellationProcessPersistence = 2,
}

/// One immutable package-policy row.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct DependencyPolicyRowV1 {
    package: &'static str,
    source: DependencySourceRouteV1,
    required_features: &'static [&'static str],
    owner: DependencyOwnerPhaseV1,
}

impl DependencyPolicyRowV1 {
    const fn new(
        package: &'static str,
        source: DependencySourceRouteV1,
        required_features: &'static [&'static str],
        owner: DependencyOwnerPhaseV1,
    ) -> Self {
        Self {
            package,
            source,
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
        self.source
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
        DependencySourceRouteV1::WorkspacePath,
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
        DependencySourceRouteV1::WorkspacePath,
        &[],
        DependencyOwnerPhaseV1::BaseSchema,
    ),
    DependencyPolicyRowV1::new(
        "asupersync",
        DependencySourceRouteV1::PinnedConstellationSibling,
        &[],
        DependencyOwnerPhaseV1::CapabilityCancellationProcessPersistence,
    ),
    DependencyPolicyRowV1::new(
        "fsqlite",
        DependencySourceRouteV1::PinnedConstellationSibling,
        &["async-api"],
        DependencyOwnerPhaseV1::CapabilityCancellationProcessPersistence,
    ),
];

/// Untrusted description of one direct Cargo route.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PresentedDependencyRouteV1 {
    package: StableTokenV2,
    source: DependencySourceRouteV1,
    features: Box<[StableTokenV2]>,
    owner: DependencyOwnerPhaseV1,
    optional: bool,
    renamed: bool,
    build_dependency: bool,
    proc_macro: bool,
    target_specific: bool,
}

impl PresentedDependencyRouteV1 {
    /// Assemble an untrusted route while validating token syntax and
    /// canonicalizing its duplicate-free feature set.
    #[allow(clippy::too_many_arguments)]
    pub fn presented(
        package: impl Into<String>,
        source: DependencySourceRouteV1,
        features: Vec<String>,
        owner: DependencyOwnerPhaseV1,
        optional: bool,
        renamed: bool,
        build_dependency: bool,
        proc_macro: bool,
        target_specific: bool,
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
        Ok(Self {
            package,
            source,
            features: features.into_boxed_slice(),
            owner,
            optional,
            renamed,
            build_dependency,
            proc_macro,
            target_specific,
        })
    }

    /// Presented package name.
    #[must_use]
    pub const fn package(&self) -> &StableTokenV2 {
        &self.package
    }

    /// Presented feature set in canonical order.
    #[must_use]
    pub fn features(&self) -> &[StableTokenV2] {
        &self.features
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
            && route.source == policy.source
            && actual_features == policy.required_features
            && route.owner == policy.owner
            && !route.optional
            && !route.renamed
            && !route.build_dependency
            && !route.proc_macro
            && !route.target_specific
            && route.source != DependencySourceRouteV1::AmbientRegistry;
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
            policy.source(),
            policy
                .required_features()
                .iter()
                .map(|feature| (*feature).to_owned())
                .collect(),
            policy.owner(),
            false,
            false,
            false,
            false,
            false,
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
                DependencySourceRouteV1::WorkspacePath,
                &[],
                DependencyOwnerPhaseV1::BaseSchema,
            )
        );
        assert_eq!(
            EVENTUAL_DIRECT_NORMAL_DEPENDENCY_ALLOWLIST_V1,
            [
                DependencyPolicyRowV1::new(
                    "fs-blake3",
                    DependencySourceRouteV1::WorkspacePath,
                    &[],
                    DependencyOwnerPhaseV1::BaseSchema,
                ),
                DependencyPolicyRowV1::new(
                    "asupersync",
                    DependencySourceRouteV1::PinnedConstellationSibling,
                    &[],
                    DependencyOwnerPhaseV1::CapabilityCancellationProcessPersistence,
                ),
                DependencyPolicyRowV1::new(
                    "fsqlite",
                    DependencySourceRouteV1::PinnedConstellationSibling,
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
    }

    #[test]
    fn every_absent_extra_or_order_mutant_refuses() {
        let exact = EVENTUAL_DIRECT_NORMAL_DEPENDENCY_ALLOWLIST_V1.map(route);
        assert!(validate_eventual_direct_dependencies_v1(&exact[..2]).is_err());
        let mut extra = exact.to_vec();
        extra.push(
            PresentedDependencyRouteV1::presented(
                "other",
                DependencySourceRouteV1::WorkspacePath,
                Vec::new(),
                DependencyOwnerPhaseV1::BaseSchema,
                false,
                false,
                false,
                false,
                false,
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
        for flags in [
            (true, false, false, false, false),
            (false, true, false, false, false),
            (false, false, true, false, false),
            (false, false, false, true, false),
            (false, false, false, false, true),
        ] {
            let mutant = PresentedDependencyRouteV1::presented(
                expected.package(),
                expected.source(),
                Vec::new(),
                expected.owner(),
                flags.0,
                flags.1,
                flags.2,
                flags.3,
                flags.4,
            )
            .expect("flag mutant");
            assert!(validate_current_direct_dependencies_v1(&[mutant]).is_err());
        }
        for mutant in [
            PresentedDependencyRouteV1::presented(
                expected.package(),
                DependencySourceRouteV1::AmbientRegistry,
                Vec::new(),
                expected.owner(),
                false,
                false,
                false,
                false,
                false,
            )
            .expect("registry mutant"),
            PresentedDependencyRouteV1::presented(
                expected.package(),
                expected.source(),
                vec!["unexpected".to_owned()],
                expected.owner(),
                false,
                false,
                false,
                false,
                false,
            )
            .expect("feature mutant"),
            PresentedDependencyRouteV1::presented(
                expected.package(),
                expected.source(),
                Vec::new(),
                DependencyOwnerPhaseV1::CapabilityCancellationProcessPersistence,
                false,
                false,
                false,
                false,
                false,
            )
            .expect("owner mutant"),
        ] {
            assert!(validate_current_direct_dependencies_v1(&[mutant]).is_err());
        }
    }

    #[test]
    fn live_phase_one_manifest_has_only_the_owned_normal_row() {
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
        assert!(!manifest.contains("[target."));
        assert!(!manifest.contains("optional = true"));
        assert!(!manifest.contains("package = "));
    }
}
