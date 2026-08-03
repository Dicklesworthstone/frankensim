//! Typed Euler-disc composition declarations and exactly-once energy receipts.
//!
//! This module is deliberately an integration boundary, not a contact, base,
//! gas, or dissipation model. It admits only declarations whose ownership
//! domains cannot overlap ambiguously, and it retains caller-supplied energy
//! terms exactly once per contribution identity. It never infers a force,
//! impulse, power, or energy closure from a port declaration.

use core::cmp::Ordering;
use core::fmt;
use std::collections::BTreeSet;

use fs_couple::{CoordinateBinding, PortKind, PortTimestamp, StableId};

/// Version of this Euler-local composition declaration surface.
pub const EULER_COMPOSITION_PORT_SCHEMA_VERSION: u16 = 1;

/// Maximum number of declared ports in one Euler composition registry.
pub const MAX_EULER_PORT_DECLARATIONS: usize = 64;

/// Typed channel families that may eventually contribute to Euler-disc motion.
///
/// Enum order is the canonical deterministic channel order.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum EulerChannel {
    /// Uniform or spatially declared gravity contribution.
    Gravity,
    /// Normal contact contribution.
    NormalContact,
    /// Tangential, partial-slip, creepage, or microslip contribution.
    TangentialContact,
    /// Rolling resistance, contour, and spin contribution.
    RollingContourSpin,
    /// Instantaneous impact contribution.
    Impact,
    /// Flexible or rigid base contribution.
    Base,
    /// External gas contribution outside the thin film.
    ExternalGas,
    /// Thin-gap gas-film contribution.
    GasFilm,
}

impl EulerChannel {
    /// All channel families in canonical order.
    pub const ALL: [Self; 8] = [
        Self::Gravity,
        Self::NormalContact,
        Self::TangentialContact,
        Self::RollingContourSpin,
        Self::Impact,
        Self::Base,
        Self::ExternalGas,
        Self::GasFilm,
    ];
}

/// Whether a declared channel may contribute in its declared interval.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ChannelActivity {
    /// The contribution is declared active and can own an energy receipt.
    Active,
    /// The model is present but declared inactive for this interval.
    Inactive {
        /// Identity of the model/card making the inactive declaration.
        model_identity: StableId,
    },
    /// The model is not available in this composition.
    Unavailable {
        /// Identity of the requested but unavailable model/card.
        model_identity: StableId,
        /// Stable identity of the missing capability or refusal source.
        reason_identity: StableId,
    },
}

impl ChannelActivity {
    /// Whether this declaration can receive a ledger contribution.
    #[must_use]
    pub const fn is_active(&self) -> bool {
        matches!(self, Self::Active)
    }
}

/// Two named surfaces participating in one canonical interface pair.
///
/// The constructor sorts the two identities, so action/reaction descriptions
/// cannot obtain separate ownership domains merely by reversing their order.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub struct SurfacePair {
    first: StableId,
    second: StableId,
}

impl SurfacePair {
    /// Creates a canonical pair of distinct surface identities.
    pub fn try_new(first: StableId, second: StableId) -> Result<Self, EulerPortError> {
        if first == second {
            return Err(EulerPortError::IdenticalSurfacePair);
        }
        if first < second {
            Ok(Self { first, second })
        } else {
            Ok(Self {
                first: second,
                second: first,
            })
        }
    }

    /// First canonical surface identity.
    #[must_use]
    pub fn first(&self) -> &StableId {
        &self.first
    }

    /// Second canonical surface identity.
    #[must_use]
    pub fn second(&self) -> &StableId {
        &self.second
    }
}

/// A half-open region interval inside one named surface patch.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PatchRegion {
    patch_identity: StableId,
    first: u64,
    end_exclusive: u64,
}

impl PatchRegion {
    /// Creates a non-empty half-open patch region `[first, end_exclusive)`.
    pub fn try_new(
        patch_identity: StableId,
        first: u64,
        end_exclusive: u64,
    ) -> Result<Self, EulerPortError> {
        if first >= end_exclusive {
            return Err(EulerPortError::InvalidPatchRegion {
                first,
                end_exclusive,
            });
        }
        Ok(Self {
            patch_identity,
            first,
            end_exclusive,
        })
    }

    /// Stable patch/region identity.
    #[must_use]
    pub fn patch_identity(&self) -> &StableId {
        &self.patch_identity
    }

    /// Inclusive first logical patch coordinate.
    #[must_use]
    pub const fn first(&self) -> u64 {
        self.first
    }

    /// Exclusive final logical patch coordinate.
    #[must_use]
    pub const fn end_exclusive(&self) -> u64 {
        self.end_exclusive
    }

    fn overlaps(&self, other: &Self) -> bool {
        self.patch_identity == other.patch_identity
            && self.first < other.end_exclusive
            && other.first < self.end_exclusive
    }
}

/// A half-open interval in one named logical clock domain.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PortInterval {
    start: PortTimestamp,
    end: PortTimestamp,
}

impl PortInterval {
    /// Creates a non-empty interval with a single clock identity.
    pub fn try_new(start: PortTimestamp, end: PortTimestamp) -> Result<Self, EulerPortError> {
        if start.clock() != end.clock() || start.tick() >= end.tick() {
            return Err(EulerPortError::InvalidPortInterval);
        }
        Ok(Self { start, end })
    }

    /// Interval start, inclusive.
    #[must_use]
    pub fn start(&self) -> &PortTimestamp {
        &self.start
    }

    /// Interval end, exclusive.
    #[must_use]
    pub fn end(&self) -> &PortTimestamp {
        &self.end
    }

    fn overlaps(&self, other: &Self) -> bool {
        self.start.clock() == other.start.clock()
            && self.start.tick() < other.end.tick()
            && other.start.tick() < self.end.tick()
    }
}

/// A named generalized velocity coordinate and its frame/sign binding.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct GeneralizedVelocityCoordinate {
    identity: StableId,
    binding: CoordinateBinding,
}

impl GeneralizedVelocityCoordinate {
    /// Binds a coordinate identity to basis, frame, and positive orientation.
    #[must_use]
    pub fn new(identity: StableId, binding: CoordinateBinding) -> Self {
        Self { identity, binding }
    }

    /// Stable coordinate identity.
    #[must_use]
    pub fn identity(&self) -> &StableId {
        &self.identity
    }

    /// Basis/frame/sign declaration used by this coordinate.
    #[must_use]
    pub fn binding(&self) -> &CoordinateBinding {
        &self.binding
    }
}

/// The spatial, temporal, and generalized-coordinate ownership domain of one port.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ContributionDomain {
    surface_pair: SurfacePair,
    patch_region: PatchRegion,
    interval: PortInterval,
    coordinate: GeneralizedVelocityCoordinate,
}

impl ContributionDomain {
    /// Creates one complete ownership domain.
    #[must_use]
    pub fn new(
        surface_pair: SurfacePair,
        patch_region: PatchRegion,
        interval: PortInterval,
        coordinate: GeneralizedVelocityCoordinate,
    ) -> Self {
        Self {
            surface_pair,
            patch_region,
            interval,
            coordinate,
        }
    }

    /// Canonical surface pair.
    #[must_use]
    pub fn surface_pair(&self) -> &SurfacePair {
        &self.surface_pair
    }

    /// Declared patch region.
    #[must_use]
    pub fn patch_region(&self) -> &PatchRegion {
        &self.patch_region
    }

    /// Declared temporal interval.
    #[must_use]
    pub fn interval(&self) -> &PortInterval {
        &self.interval
    }

    /// Declared generalized coordinate and frame/sign convention.
    #[must_use]
    pub fn coordinate(&self) -> &GeneralizedVelocityCoordinate {
        &self.coordinate
    }

    fn overlaps(&self, other: &Self) -> bool {
        self.surface_pair == other.surface_pair
            && self.patch_region.overlaps(&other.patch_region)
            && self.interval.overlaps(&other.interval)
            // The coordinate identity denotes the physical generalized degree
            // of freedom. A contradictory frame/sign binding must therefore
            // collide rather than provide a second ownership escape hatch.
            && self.coordinate.identity == other.coordinate.identity
    }
}

/// A structural additive-decomposition receipt.
///
/// This receipt binds one exact domain and the complete contributor set. It is
/// not signed evidence and cannot authenticate an external decomposition;
/// that future authority remains outside this skeleton.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DecompositionReceipt {
    identity: StableId,
    domain: ContributionDomain,
    contributor_port_ids: Vec<StableId>,
}

impl DecompositionReceipt {
    /// Creates a canonical structural decomposition receipt for two or more contributors.
    pub fn try_new(
        identity: StableId,
        domain: ContributionDomain,
        contributor_port_ids: impl IntoIterator<Item = StableId>,
    ) -> Result<Self, EulerPortError> {
        let mut contributor_port_ids: Vec<_> = contributor_port_ids.into_iter().collect();
        if contributor_port_ids.len() < 2 {
            return Err(EulerPortError::IncompleteDecompositionReceipt);
        }
        contributor_port_ids.sort();
        let original_len = contributor_port_ids.len();
        contributor_port_ids.dedup();
        if contributor_port_ids.len() != original_len {
            return Err(EulerPortError::DuplicateDecompositionContributor);
        }
        Ok(Self {
            identity,
            domain,
            contributor_port_ids,
        })
    }

    /// Stable decomposition receipt identity.
    #[must_use]
    pub fn identity(&self) -> &StableId {
        &self.identity
    }

    fn covers(&self, first: &PortDeclaration, second: &PortDeclaration) -> bool {
        self.domain == first.domain
            && self.domain == second.domain
            && self
                .contributor_port_ids
                .binary_search(&first.identity)
                .is_ok()
            && self
                .contributor_port_ids
                .binary_search(&second.identity)
                .is_ok()
    }
}

/// Ownership discipline for one port contribution.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum ContributionOwnership {
    /// This declaration is the sole owner of its domain.
    Exclusive,
    /// Multiple declarations share exactly one domain under one receipt.
    AdditiveWithProof {
        /// Structural decomposition receipt binding all contributors.
        decomposition_receipt: DecompositionReceipt,
    },
}

/// Typed declaration for one possible mechanical contribution channel.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PortDeclaration {
    identity: StableId,
    channel: EulerChannel,
    port_kind: PortKind,
    activity: ChannelActivity,
    law_identity: StableId,
    source_identity: StableId,
    domain: ContributionDomain,
    ownership: ContributionOwnership,
}

impl PortDeclaration {
    /// Creates one complete channel declaration.
    #[allow(clippy::too_many_arguments)]
    #[must_use]
    pub fn new(
        identity: StableId,
        channel: EulerChannel,
        port_kind: PortKind,
        activity: ChannelActivity,
        law_identity: StableId,
        source_identity: StableId,
        domain: ContributionDomain,
        ownership: ContributionOwnership,
    ) -> Self {
        Self {
            identity,
            channel,
            port_kind,
            activity,
            law_identity,
            source_identity,
            domain,
            ownership,
        }
    }

    /// Stable port declaration identity.
    #[must_use]
    pub fn identity(&self) -> &StableId {
        &self.identity
    }

    /// Typed channel family.
    #[must_use]
    pub const fn channel(&self) -> EulerChannel {
        self.channel
    }

    /// Mechanical effort/flow vocabulary used by this contribution.
    #[must_use]
    pub const fn port_kind(&self) -> PortKind {
        self.port_kind
    }

    /// Active, inactive, or unavailable state.
    #[must_use]
    pub fn activity(&self) -> &ChannelActivity {
        &self.activity
    }

    /// Constitutive-law/card/model identity.
    #[must_use]
    pub fn law_identity(&self) -> &StableId {
        &self.law_identity
    }

    /// Provenance/source identity for the declaration.
    #[must_use]
    pub fn source_identity(&self) -> &StableId {
        &self.source_identity
    }

    /// Complete force/moment ownership domain.
    #[must_use]
    pub fn domain(&self) -> &ContributionDomain {
        &self.domain
    }

    /// Ownership mode.
    #[must_use]
    pub fn ownership(&self) -> &ContributionOwnership {
        &self.ownership
    }
}

/// Deterministic, overlap-refusing Euler composition registry.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct EulerPortRegistry {
    identity: StableId,
    declarations: Vec<PortDeclaration>,
}

impl EulerPortRegistry {
    /// Admits a deterministic registry with non-overlapping active ownership.
    pub fn try_new(
        identity: StableId,
        declarations: impl IntoIterator<Item = PortDeclaration>,
    ) -> Result<Self, EulerPortError> {
        let mut declarations: Vec<_> = declarations.into_iter().collect();
        if declarations.len() > MAX_EULER_PORT_DECLARATIONS {
            return Err(EulerPortError::TooManyPortDeclarations {
                maximum: MAX_EULER_PORT_DECLARATIONS,
                actual: declarations.len(),
            });
        }
        declarations.sort_by(canonical_declaration_order);
        for pair in declarations.windows(2) {
            let [first, second] = pair else {
                continue;
            };
            if first.identity == second.identity {
                return Err(EulerPortError::DuplicatePortIdentity {
                    identity: first.identity.clone(),
                });
            }
        }
        for declaration in &declarations {
            let ContributionOwnership::AdditiveWithProof {
                decomposition_receipt,
            } = &declaration.ownership
            else {
                continue;
            };
            if decomposition_receipt.domain != declaration.domain
                || decomposition_receipt
                    .contributor_port_ids
                    .binary_search(&declaration.identity)
                    .is_err()
            {
                return Err(EulerPortError::AdditiveProofMismatch {
                    first: declaration.identity.clone(),
                    second: decomposition_receipt.identity.clone(),
                });
            }
            for contributor_identity in &decomposition_receipt.contributor_port_ids {
                let Some(contributor_index) = declarations
                    .binary_search_by(|candidate| candidate.identity.cmp(contributor_identity))
                    .ok()
                else {
                    return Err(EulerPortError::AdditiveProofMismatch {
                        first: declaration.identity.clone(),
                        second: contributor_identity.clone(),
                    });
                };
                let Some(contributor) = declarations.get(contributor_index) else {
                    return Err(EulerPortError::AdditiveProofMismatch {
                        first: declaration.identity.clone(),
                        second: contributor_identity.clone(),
                    });
                };
                if !contributor.activity.is_active()
                    || contributor.domain != decomposition_receipt.domain
                    || !matches!(
                        &contributor.ownership,
                        ContributionOwnership::AdditiveWithProof {
                            decomposition_receipt: contributor_receipt,
                        } if contributor_receipt == decomposition_receipt
                    )
                {
                    return Err(EulerPortError::AdditiveProofMismatch {
                        first: declaration.identity.clone(),
                        second: contributor_identity.clone(),
                    });
                }
            }
        }
        for (index, left) in declarations.iter().enumerate() {
            if !left.activity.is_active() {
                continue;
            }
            for right in declarations.iter().skip(index + 1) {
                if !right.activity.is_active() || !left.domain.overlaps(&right.domain) {
                    continue;
                }
                match (&left.ownership, &right.ownership) {
                    (
                        ContributionOwnership::AdditiveWithProof {
                            decomposition_receipt: left_receipt,
                        },
                        ContributionOwnership::AdditiveWithProof {
                            decomposition_receipt: right_receipt,
                        },
                    ) if left_receipt == right_receipt && left_receipt.covers(left, right) => {}
                    (ContributionOwnership::AdditiveWithProof { .. }, _)
                    | (_, ContributionOwnership::AdditiveWithProof { .. }) => {
                        return Err(EulerPortError::AdditiveProofMismatch {
                            first: left.identity.clone(),
                            second: right.identity.clone(),
                        });
                    }
                    _ => {
                        return Err(EulerPortError::OverlappingExclusiveOwnership {
                            first: left.identity.clone(),
                            second: right.identity.clone(),
                        });
                    }
                }
            }
        }
        Ok(Self {
            identity,
            declarations,
        })
    }

    /// Caller-declared stable registry/checkpoint identity.
    #[must_use]
    pub fn identity(&self) -> &StableId {
        &self.identity
    }

    /// Canonically ordered port declarations.
    #[must_use]
    pub fn declarations(&self) -> &[PortDeclaration] {
        &self.declarations
    }

    fn declaration(&self, identity: &StableId) -> Option<&PortDeclaration> {
        self.declarations
            .binary_search_by(|candidate| candidate.identity.cmp(identity))
            .ok()
            .and_then(|index| self.declarations.get(index))
    }

    fn unavailable_channels(&self) -> Vec<EulerChannel> {
        let mut channels = BTreeSet::new();
        for declaration in &self.declarations {
            if matches!(declaration.activity, ChannelActivity::Unavailable { .. }) {
                channels.insert(declaration.channel);
            }
        }
        channels.into_iter().collect()
    }
}

fn canonical_declaration_order(left: &PortDeclaration, right: &PortDeclaration) -> Ordering {
    left.identity
        .cmp(&right.identity)
        .then_with(|| left.channel.cmp(&right.channel))
        .then_with(|| {
            left.domain
                .surface_pair
                .first
                .cmp(&right.domain.surface_pair.first)
        })
        .then_with(|| {
            left.domain
                .surface_pair
                .second
                .cmp(&right.domain.surface_pair.second)
        })
        .then_with(|| {
            left.domain
                .patch_region
                .patch_identity
                .cmp(&right.domain.patch_region.patch_identity)
        })
        .then_with(|| {
            left.domain
                .patch_region
                .first
                .cmp(&right.domain.patch_region.first)
        })
        .then_with(|| {
            left.domain
                .patch_region
                .end_exclusive
                .cmp(&right.domain.patch_region.end_exclusive)
        })
        .then_with(|| {
            left.domain
                .interval
                .start
                .clock()
                .cmp(right.domain.interval.start.clock())
        })
        .then_with(|| {
            left.domain
                .interval
                .start
                .tick()
                .cmp(&right.domain.interval.start.tick())
        })
}

/// Signed storage changes and non-negative loss/uncertainty terms for one event.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EnergyTerms {
    /// Signed kinetic-energy change in joules.
    kinetic_j: f64,
    /// Signed potential-energy change in joules.
    potential_j: f64,
    /// Signed recoverable-storage change in joules.
    recoverable_j: f64,
    /// Non-negative dissipated mechanical energy in joules.
    dissipated_j: f64,
    /// Non-negative heat energy in joules.
    heat_j: f64,
    /// Signed numerical balance defect in joules.
    numerical_defect_j: f64,
    /// Non-negative unresolved-energy magnitude in joules.
    unresolved_j: f64,
}

impl EnergyTerms {
    /// Zero energy contribution.
    pub const ZERO: Self = Self {
        kinetic_j: 0.0,
        potential_j: 0.0,
        recoverable_j: 0.0,
        dissipated_j: 0.0,
        heat_j: 0.0,
        numerical_defect_j: 0.0,
        unresolved_j: 0.0,
    };

    /// Validates one caller-supplied accounting contribution.
    #[allow(clippy::too_many_arguments)]
    pub fn try_new(
        kinetic_j: f64,
        potential_j: f64,
        recoverable_j: f64,
        dissipated_j: f64,
        heat_j: f64,
        numerical_defect_j: f64,
        unresolved_j: f64,
    ) -> Result<Self, EulerPortError> {
        let values = [
            kinetic_j,
            potential_j,
            recoverable_j,
            dissipated_j,
            heat_j,
            numerical_defect_j,
            unresolved_j,
        ];
        if values.iter().any(|value| !value.is_finite()) {
            return Err(EulerPortError::NonFiniteEnergyTerms);
        }
        if dissipated_j < 0.0 || heat_j < 0.0 || unresolved_j < 0.0 {
            return Err(EulerPortError::NegativeLossOrUnresolvedEnergy);
        }
        Ok(Self {
            kinetic_j,
            potential_j,
            recoverable_j,
            dissipated_j,
            heat_j,
            numerical_defect_j,
            unresolved_j,
        })
    }

    /// Signed kinetic-energy change in joules.
    #[must_use]
    pub const fn kinetic_j(self) -> f64 {
        self.kinetic_j
    }

    /// Signed potential-energy change in joules.
    #[must_use]
    pub const fn potential_j(self) -> f64 {
        self.potential_j
    }

    /// Signed recoverable-storage change in joules.
    #[must_use]
    pub const fn recoverable_j(self) -> f64 {
        self.recoverable_j
    }

    /// Non-negative dissipated mechanical energy in joules.
    #[must_use]
    pub const fn dissipated_j(self) -> f64 {
        self.dissipated_j
    }

    /// Non-negative heat energy in joules.
    #[must_use]
    pub const fn heat_j(self) -> f64 {
        self.heat_j
    }

    /// Signed numerical balance defect in joules.
    #[must_use]
    pub const fn numerical_defect_j(self) -> f64 {
        self.numerical_defect_j
    }

    /// Non-negative unresolved-energy magnitude in joules.
    #[must_use]
    pub const fn unresolved_j(self) -> f64 {
        self.unresolved_j
    }

    fn checked_add(self, other: Self) -> Result<Self, EulerPortError> {
        Self::try_new(
            finite_energy_sum(self.kinetic_j, other.kinetic_j)?,
            finite_energy_sum(self.potential_j, other.potential_j)?,
            finite_energy_sum(self.recoverable_j, other.recoverable_j)?,
            finite_energy_sum(self.dissipated_j, other.dissipated_j)?,
            finite_energy_sum(self.heat_j, other.heat_j)?,
            finite_energy_sum(self.numerical_defect_j, other.numerical_defect_j)?,
            finite_energy_sum(self.unresolved_j, other.unresolved_j)?,
        )
    }
}

/// One exactly-once energy contribution owned by one active port declaration.
#[derive(Clone, Debug, PartialEq)]
pub struct EnergyContribution {
    identity: StableId,
    port_identity: StableId,
    channel: EulerChannel,
    timestamp: PortTimestamp,
    terms: EnergyTerms,
}

impl EnergyContribution {
    /// Creates one caller-supplied energy contribution.
    #[must_use]
    pub fn new(
        identity: StableId,
        port_identity: StableId,
        channel: EulerChannel,
        timestamp: PortTimestamp,
        terms: EnergyTerms,
    ) -> Self {
        Self {
            identity,
            port_identity,
            channel,
            timestamp,
            terms,
        }
    }

    /// Exactly-once contribution identity.
    #[must_use]
    pub fn identity(&self) -> &StableId {
        &self.identity
    }

    /// Owning port declaration identity.
    #[must_use]
    pub fn port_identity(&self) -> &StableId {
        &self.port_identity
    }

    /// Typed channel asserted by this contribution.
    #[must_use]
    pub const fn channel(&self) -> EulerChannel {
        self.channel
    }

    /// Logical event timestamp.
    #[must_use]
    pub fn timestamp(&self) -> &PortTimestamp {
        &self.timestamp
    }

    /// Accounting terms.
    #[must_use]
    pub const fn terms(&self) -> EnergyTerms {
        self.terms
    }
}

/// Checkpoint bound to one ledger identity, registry identity, and retained prefix.
#[derive(Clone, Debug, PartialEq)]
pub struct EnergyLedgerCheckpoint {
    identity: StableId,
    ledger_identity: StableId,
    registry_identity: StableId,
    entry_count: usize,
    final_entry_identity: Option<StableId>,
    timestamp: PortTimestamp,
    cumulative: EnergyTerms,
}

impl EnergyLedgerCheckpoint {
    /// Stable checkpoint identity.
    #[must_use]
    pub fn identity(&self) -> &StableId {
        &self.identity
    }

    /// Retained entry count in this checkpoint.
    #[must_use]
    pub const fn entry_count(&self) -> usize {
        self.entry_count
    }

    /// Exact cumulative terms at this checkpoint.
    #[must_use]
    pub const fn cumulative(&self) -> EnergyTerms {
        self.cumulative
    }

    /// Logical checkpoint timestamp.
    #[must_use]
    pub fn timestamp(&self) -> &PortTimestamp {
        &self.timestamp
    }
}

/// Explicit non-closure result for this skeleton.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum EnergyClosureDisposition {
    /// One or more requested channel families are unavailable.
    NoClaimUnavailableChannels {
        /// Canonically ordered unavailable channel families.
        channels: Vec<EulerChannel>,
    },
    /// All declared channels are active/inactive, but this skeleton has no
    /// physics or closed-window audit authority.
    NoClaimIntegrationSkeleton,
}

/// Exactly-once cumulative energy accounting for an admitted port registry.
#[derive(Clone, Debug, PartialEq)]
pub struct EulerEnergyLedger {
    identity: StableId,
    registry: EulerPortRegistry,
    contributions: Vec<EnergyContribution>,
    contribution_ids: BTreeSet<StableId>,
    cumulative: EnergyTerms,
}

impl EulerEnergyLedger {
    /// Creates an empty ledger bound to one admitted registry.
    #[must_use]
    pub fn new(identity: StableId, registry: EulerPortRegistry) -> Self {
        Self {
            identity,
            registry,
            contributions: Vec::new(),
            contribution_ids: BTreeSet::new(),
            cumulative: EnergyTerms::ZERO,
        }
    }

    /// Stable ledger identity.
    #[must_use]
    pub fn identity(&self) -> &StableId {
        &self.identity
    }

    /// Registry bound to this ledger.
    #[must_use]
    pub fn registry(&self) -> &EulerPortRegistry {
        &self.registry
    }

    /// Retained contributions in exactly-once insertion order.
    #[must_use]
    pub fn contributions(&self) -> &[EnergyContribution] {
        &self.contributions
    }

    /// Cumulative signed stores and non-negative losses/uncertainty magnitude.
    #[must_use]
    pub const fn cumulative(&self) -> EnergyTerms {
        self.cumulative
    }

    /// Records one contribution once, without partial mutation on refusal.
    pub fn record(&mut self, contribution: EnergyContribution) -> Result<(), EulerPortError> {
        if self.contribution_ids.contains(&contribution.identity) {
            return Err(EulerPortError::DuplicateEnergyContribution {
                identity: contribution.identity,
            });
        }
        let port = self
            .registry
            .declaration(&contribution.port_identity)
            .ok_or_else(|| EulerPortError::UnknownEnergyPort {
                identity: contribution.port_identity.clone(),
            })?;
        if contribution.channel != port.channel {
            return Err(EulerPortError::EnergyChannelMismatch {
                port: port.identity.clone(),
                expected: port.channel,
                actual: contribution.channel,
            });
        }
        if contribution.timestamp.clock() != port.domain.interval.start.clock()
            || contribution.timestamp.tick() < port.domain.interval.start.tick()
            || contribution.timestamp.tick() >= port.domain.interval.end.tick()
        {
            return Err(EulerPortError::EnergyTimestampOutsidePortInterval {
                port: port.identity.clone(),
            });
        }
        match &port.activity {
            ChannelActivity::Active => {}
            ChannelActivity::Inactive { .. } => {
                return Err(EulerPortError::EnergyPortInactive {
                    identity: port.identity.clone(),
                });
            }
            ChannelActivity::Unavailable { .. } => {
                return Err(EulerPortError::EnergyPortUnavailable {
                    identity: port.identity.clone(),
                });
            }
        }
        if let Some(previous) = self.contributions.last()
            && canonical_energy_contribution_order(previous, &contribution) == Ordering::Greater
        {
            return Err(EulerPortError::NonDeterministicEnergyContributionOrder);
        }
        let next_cumulative = self.cumulative.checked_add(contribution.terms)?;
        self.contribution_ids.insert(contribution.identity.clone());
        self.contributions.push(contribution);
        self.cumulative = next_cumulative;
        Ok(())
    }

    /// Creates an identity-bound checkpoint for deterministic rollback.
    #[must_use]
    pub fn checkpoint(
        &self,
        identity: StableId,
        timestamp: PortTimestamp,
    ) -> EnergyLedgerCheckpoint {
        EnergyLedgerCheckpoint {
            identity,
            ledger_identity: self.identity.clone(),
            registry_identity: self.registry.identity.clone(),
            entry_count: self.contributions.len(),
            final_entry_identity: self
                .contributions
                .last()
                .map(|entry| entry.identity.clone()),
            timestamp,
            cumulative: self.cumulative,
        }
    }

    /// Rolls back to one checkpoint from this exact ledger and registry.
    ///
    /// The retained prefix is recomputed before mutation, so an invalid or
    /// stale checkpoint cannot partially alter the visible ledger.
    pub fn rollback(&mut self, checkpoint: &EnergyLedgerCheckpoint) -> Result<(), EulerPortError> {
        if checkpoint.ledger_identity != self.identity
            || checkpoint.registry_identity != self.registry.identity
            || checkpoint.entry_count > self.contributions.len()
        {
            return Err(EulerPortError::CheckpointIdentityMismatch);
        }
        let expected_final = checkpoint.entry_count.checked_sub(1).and_then(|index| {
            self.contributions
                .get(index)
                .map(|contribution| contribution.identity.clone())
        });
        if expected_final != checkpoint.final_entry_identity {
            return Err(EulerPortError::CheckpointIdentityMismatch);
        }
        let retained = self.contributions[..checkpoint.entry_count].to_vec();
        let cumulative = retained
            .iter()
            .try_fold(EnergyTerms::ZERO, |sum, contribution| {
                sum.checked_add(contribution.terms)
            })?;
        if cumulative != checkpoint.cumulative {
            return Err(EulerPortError::CheckpointIdentityMismatch);
        }
        let contribution_ids = retained
            .iter()
            .map(|contribution| contribution.identity.clone())
            .collect();
        self.contributions = retained;
        self.contribution_ids = contribution_ids;
        self.cumulative = cumulative;
        Ok(())
    }

    /// Returns a no-claim disposition; this module never certifies closure.
    #[must_use]
    pub fn closure_disposition(&self) -> EnergyClosureDisposition {
        let unavailable_channels = self.registry.unavailable_channels();
        if unavailable_channels.is_empty() {
            EnergyClosureDisposition::NoClaimIntegrationSkeleton
        } else {
            EnergyClosureDisposition::NoClaimUnavailableChannels {
                channels: unavailable_channels,
            }
        }
    }
}

fn finite_energy_sum(left: f64, right: f64) -> Result<f64, EulerPortError> {
    let sum = left + right;
    if sum.is_finite() {
        Ok(sum)
    } else {
        Err(EulerPortError::NonFiniteCumulativeEnergy)
    }
}

fn canonical_energy_contribution_order(
    left: &EnergyContribution,
    right: &EnergyContribution,
) -> Ordering {
    left.timestamp
        .clock()
        .cmp(right.timestamp.clock())
        .then_with(|| left.timestamp.tick().cmp(&right.timestamp.tick()))
        .then_with(|| left.identity.cmp(&right.identity))
}

/// Structural refusal returned by Euler composition declarations and ledger operations.
#[derive(Clone, Debug, PartialEq)]
pub enum EulerPortError {
    /// A surface pair named the same surface twice.
    IdenticalSurfacePair,
    /// A patch region was empty or backwards.
    InvalidPatchRegion {
        /// Declared inclusive start.
        first: u64,
        /// Declared exclusive end.
        end_exclusive: u64,
    },
    /// A port interval mixed clocks or was empty/backwards.
    InvalidPortInterval,
    /// A decomposition receipt named fewer than two contributors.
    IncompleteDecompositionReceipt,
    /// A decomposition receipt repeated a contributor identity.
    DuplicateDecompositionContributor,
    /// Registry size exceeded its declared resource limit.
    TooManyPortDeclarations {
        /// Maximum admitted declarations.
        maximum: usize,
        /// Attempted declaration count.
        actual: usize,
    },
    /// Two declarations used one port identity.
    DuplicatePortIdentity {
        /// Duplicate identity.
        identity: StableId,
    },
    /// Two active exclusive owners overlap in surface, patch, time, and coordinate.
    OverlappingExclusiveOwnership {
        /// First conflicting declaration identity.
        first: StableId,
        /// Second conflicting declaration identity.
        second: StableId,
    },
    /// An additive overlap lacked one exact structural decomposition proof.
    AdditiveProofMismatch {
        /// First conflicting declaration identity.
        first: StableId,
        /// Second conflicting declaration identity.
        second: StableId,
    },
    /// Any energy term was NaN or infinite.
    NonFiniteEnergyTerms,
    /// A loss or unresolved-energy magnitude was negative.
    NegativeLossOrUnresolvedEnergy,
    /// A cumulative energy term overflowed or became non-finite.
    NonFiniteCumulativeEnergy,
    /// The contribution identity was already retained by this ledger timeline.
    DuplicateEnergyContribution {
        /// Duplicate contribution identity.
        identity: StableId,
    },
    /// A contribution referred to no declaration in this registry.
    UnknownEnergyPort {
        /// Unknown port identity.
        identity: StableId,
    },
    /// A contribution channel differed from its port declaration.
    EnergyChannelMismatch {
        /// Port identity.
        port: StableId,
        /// Declared channel.
        expected: EulerChannel,
        /// Contribution channel.
        actual: EulerChannel,
    },
    /// A contribution timestamp was outside its declared half-open interval.
    EnergyTimestampOutsidePortInterval {
        /// Port identity.
        port: StableId,
    },
    /// An inactive declaration cannot receive energy.
    EnergyPortInactive {
        /// Port identity.
        identity: StableId,
    },
    /// An unavailable declaration cannot receive energy.
    EnergyPortUnavailable {
        /// Port identity.
        identity: StableId,
    },
    /// A checkpoint did not originate from this exact ledger prefix.
    CheckpointIdentityMismatch,
    /// A contribution was not appended in canonical timestamp/identity order.
    NonDeterministicEnergyContributionOrder,
}

impl fmt::Display for EulerPortError {
    fn fmt(&self, formatter: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::IdenticalSurfacePair => {
                formatter.write_str("surface pair must name two surfaces")
            }
            Self::InvalidPatchRegion { .. } => {
                formatter.write_str("patch region must be non-empty")
            }
            Self::InvalidPortInterval => {
                formatter.write_str("port interval must be non-empty and use one clock")
            }
            Self::IncompleteDecompositionReceipt => formatter
                .write_str("additive decomposition receipt needs at least two contributors"),
            Self::DuplicateDecompositionContributor => {
                formatter.write_str("additive decomposition receipt has a duplicate contributor")
            }
            Self::TooManyPortDeclarations { maximum, actual } => {
                write!(
                    formatter,
                    "port declaration count {actual} exceeds {maximum}"
                )
            }
            Self::DuplicatePortIdentity { identity } => {
                write!(formatter, "duplicate port identity {}", identity.as_str())
            }
            Self::OverlappingExclusiveOwnership { first, second } => write!(
                formatter,
                "exclusive port ownership overlaps: {} and {}",
                first.as_str(),
                second.as_str()
            ),
            Self::AdditiveProofMismatch { first, second } => write!(
                formatter,
                "additive decomposition proof mismatch: {} and {}",
                first.as_str(),
                second.as_str()
            ),
            Self::NonFiniteEnergyTerms => formatter.write_str("energy terms must be finite"),
            Self::NegativeLossOrUnresolvedEnergy => {
                formatter.write_str("loss and unresolved-energy magnitudes must be non-negative")
            }
            Self::NonFiniteCumulativeEnergy => {
                formatter.write_str("cumulative energy must be finite")
            }
            Self::DuplicateEnergyContribution { identity } => {
                write!(
                    formatter,
                    "duplicate energy contribution {}",
                    identity.as_str()
                )
            }
            Self::UnknownEnergyPort { identity } => {
                write!(formatter, "unknown energy port {}", identity.as_str())
            }
            Self::EnergyChannelMismatch { port, .. } => {
                write!(
                    formatter,
                    "energy channel differs from port {}",
                    port.as_str()
                )
            }
            Self::EnergyTimestampOutsidePortInterval { port } => {
                write!(
                    formatter,
                    "energy timestamp lies outside port {} interval",
                    port.as_str()
                )
            }
            Self::EnergyPortInactive { identity } => {
                write!(
                    formatter,
                    "inactive port {} cannot receive energy",
                    identity.as_str()
                )
            }
            Self::EnergyPortUnavailable { identity } => {
                write!(
                    formatter,
                    "unavailable port {} cannot receive energy",
                    identity.as_str()
                )
            }
            Self::CheckpointIdentityMismatch => {
                formatter.write_str("checkpoint does not match this ledger identity and prefix")
            }
            Self::NonDeterministicEnergyContributionOrder => formatter
                .write_str("energy contribution is not in canonical timestamp and identity order"),
        }
    }
}

impl std::error::Error for EulerPortError {}
