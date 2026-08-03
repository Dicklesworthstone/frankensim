//! G0/G3 tests for the Euler-local composition-port and ledger boundary.

use fs_couple::{CoordinateBinding, PortKind, PortOrientation, PortTimestamp, StableId};
use fs_euler_disc_e2e::{
    ChannelActivity, ContributionDomain, ContributionOwnership, DecompositionReceipt,
    EnergyClosureDisposition, EnergyContribution, EnergyTerms, EulerChannel, EulerEnergyLedger,
    EulerPortError, EulerPortRegistry, GeneralizedVelocityCoordinate, PatchRegion, PortDeclaration,
    PortInterval, SurfacePair,
};

fn stable(value: &str) -> StableId {
    StableId::new(value).expect("test stable identifier")
}

fn timestamp(tick: u64) -> PortTimestamp {
    PortTimestamp::new(stable("clock-main"), tick)
}

fn domain(
    first_surface: &str,
    second_surface: &str,
    patch: &str,
    patch_first: u64,
    patch_end: u64,
    interval_start: u64,
    interval_end: u64,
    coordinate: &str,
    frame: &str,
    orientation: PortOrientation,
) -> ContributionDomain {
    ContributionDomain::new(
        SurfacePair::try_new(stable(first_surface), stable(second_surface))
            .expect("distinct test surfaces"),
        PatchRegion::try_new(stable(patch), patch_first, patch_end).expect("test patch range"),
        PortInterval::try_new(timestamp(interval_start), timestamp(interval_end))
            .expect("test interval"),
        GeneralizedVelocityCoordinate::new(
            stable(coordinate),
            CoordinateBinding::new(stable("basis-world"), stable(frame), orientation),
        ),
    )
}

fn active_port(
    identity: &str,
    channel: EulerChannel,
    domain: ContributionDomain,
    ownership: ContributionOwnership,
) -> PortDeclaration {
    PortDeclaration::new(
        stable(identity),
        channel,
        match channel {
            EulerChannel::RollingContourSpin => PortKind::RotationalTorqueAngularVelocity,
            _ => PortKind::MechanicalForceVelocity,
        },
        ChannelActivity::Active,
        stable(&format!("law-{identity}")),
        stable(&format!("source-{identity}")),
        domain,
        ownership,
    )
}

fn single_active_registry() -> EulerPortRegistry {
    let domain = domain(
        "disc",
        "ground",
        "contact-patch",
        0,
        1,
        0,
        10,
        "disc-twist",
        "world",
        PortOrientation::AlongFrame,
    );
    EulerPortRegistry::try_new(
        stable("registry-main"),
        [active_port(
            "gravity-port",
            EulerChannel::Gravity,
            domain,
            ContributionOwnership::Exclusive,
        )],
    )
    .expect("one active port is admissible")
}

fn terms(kinetic_j: f64) -> EnergyTerms {
    EnergyTerms::try_new(kinetic_j, 2.0, -0.5, 0.25, 0.125, 0.0, 0.0)
        .expect("finite non-negative constrained energy terms")
}

#[test]
fn all_channels_are_typed_and_registry_order_is_permutation_invariant() {
    let mut forward = Vec::new();
    for (index, channel) in EulerChannel::ALL.into_iter().enumerate() {
        let identity = format!("port-{index}");
        forward.push(active_port(
            &identity,
            channel,
            domain(
                "disc",
                "ground",
                &format!("patch-{index}"),
                0,
                1,
                0,
                10,
                &format!("coordinate-{index}"),
                "world",
                PortOrientation::AlongFrame,
            ),
            ContributionOwnership::Exclusive,
        ));
    }
    let mut reverse = forward.clone();
    reverse.reverse();

    let forward = EulerPortRegistry::try_new(stable("registry-channels"), forward)
        .expect("all typed channels can be declared");
    let reverse = EulerPortRegistry::try_new(stable("registry-channels"), reverse)
        .expect("caller permutation is canonicalized");

    assert_eq!(forward, reverse);
    assert_eq!(
        forward
            .declarations()
            .iter()
            .map(PortDeclaration::channel)
            .collect::<Vec<_>>(),
        vec![
            EulerChannel::Gravity,
            EulerChannel::NormalContact,
            EulerChannel::TangentialContact,
            EulerChannel::RollingContourSpin,
            EulerChannel::Impact,
            EulerChannel::Base,
            EulerChannel::ExternalGas,
            EulerChannel::GasFilm,
        ]
    );
    assert_eq!(
        forward
            .declarations()
            .first()
            .expect("declaration")
            .domain()
            .coordinate()
            .binding()
            .orientation(),
        PortOrientation::AlongFrame
    );
}

#[test]
fn duplicate_and_partially_overlapping_owners_refuse() {
    let shared = domain(
        "disc",
        "ground",
        "patch",
        0,
        10,
        0,
        10,
        "twist",
        "world",
        PortOrientation::AlongFrame,
    );
    let duplicate = EulerPortRegistry::try_new(
        stable("registry-duplicate"),
        [
            active_port(
                "same-port",
                EulerChannel::Gravity,
                shared.clone(),
                ContributionOwnership::Exclusive,
            ),
            active_port(
                "same-port",
                EulerChannel::NormalContact,
                shared.clone(),
                ContributionOwnership::Exclusive,
            ),
        ],
    );
    assert!(matches!(
        duplicate,
        Err(EulerPortError::DuplicatePortIdentity { .. })
    ));

    let partial = EulerPortRegistry::try_new(
        stable("registry-partial"),
        [
            active_port(
                "first-owner",
                EulerChannel::NormalContact,
                shared,
                ContributionOwnership::Exclusive,
            ),
            active_port(
                "second-owner",
                EulerChannel::TangentialContact,
                domain(
                    "disc",
                    "ground",
                    "patch",
                    5,
                    15,
                    5,
                    15,
                    "twist",
                    "world",
                    PortOrientation::AlongFrame,
                ),
                ContributionOwnership::Exclusive,
            ),
        ],
    );
    assert!(matches!(
        partial,
        Err(EulerPortError::OverlappingExclusiveOwnership { .. })
    ));
}

#[test]
fn action_reaction_reversal_and_frame_disagreement_cannot_bypass_ownership() {
    let common = (
        "rim-patch",
        0,
        1,
        0,
        10,
        "normal-velocity",
        "world",
        PortOrientation::AlongFrame,
    );
    let action_reaction = EulerPortRegistry::try_new(
        stable("registry-action-reaction"),
        [
            active_port(
                "action",
                EulerChannel::NormalContact,
                domain(
                    "disc", "base", common.0, common.1, common.2, common.3, common.4, common.5,
                    common.6, common.7,
                ),
                ContributionOwnership::Exclusive,
            ),
            active_port(
                "reaction",
                EulerChannel::NormalContact,
                domain(
                    "base", "disc", common.0, common.1, common.2, common.3, common.4, common.5,
                    common.6, common.7,
                ),
                ContributionOwnership::Exclusive,
            ),
        ],
    );
    assert!(matches!(
        action_reaction,
        Err(EulerPortError::OverlappingExclusiveOwnership { .. })
    ));

    let frame_disagreement = EulerPortRegistry::try_new(
        stable("registry-frame-disagreement"),
        [
            active_port(
                "frame-a",
                EulerChannel::TangentialContact,
                domain(
                    "disc",
                    "base",
                    "rim-patch",
                    0,
                    1,
                    0,
                    10,
                    "tangent-velocity",
                    "world-a",
                    PortOrientation::AlongFrame,
                ),
                ContributionOwnership::Exclusive,
            ),
            active_port(
                "frame-b",
                EulerChannel::TangentialContact,
                domain(
                    "disc",
                    "base",
                    "rim-patch",
                    0,
                    1,
                    0,
                    10,
                    "tangent-velocity",
                    "world-b",
                    PortOrientation::AgainstFrame,
                ),
                ContributionOwnership::Exclusive,
            ),
        ],
    );
    assert!(matches!(
        frame_disagreement,
        Err(EulerPortError::OverlappingExclusiveOwnership { .. })
    ));
}

#[test]
fn additive_overlap_requires_one_exact_receipt() {
    let shared = domain(
        "disc",
        "base",
        "rim-patch",
        0,
        1,
        0,
        10,
        "tangent-velocity",
        "world",
        PortOrientation::AlongFrame,
    );
    let receipt = DecompositionReceipt::try_new(
        stable("receipt-additive"),
        shared.clone(),
        [stable("additive-a"), stable("additive-b")],
    )
    .expect("exact contributor set");
    let admitted = EulerPortRegistry::try_new(
        stable("registry-additive"),
        [
            active_port(
                "additive-a",
                EulerChannel::RollingContourSpin,
                shared.clone(),
                ContributionOwnership::AdditiveWithProof {
                    decomposition_receipt: receipt.clone(),
                },
            ),
            active_port(
                "additive-b",
                EulerChannel::RollingContourSpin,
                shared.clone(),
                ContributionOwnership::AdditiveWithProof {
                    decomposition_receipt: receipt,
                },
            ),
        ],
    );
    assert!(admitted.is_ok());

    let incomplete = DecompositionReceipt::try_new(
        stable("receipt-incomplete"),
        shared,
        [stable("additive-a"), stable("wrong-port")],
    )
    .expect("structurally shaped but semantically wrong receipt");
    let refused = EulerPortRegistry::try_new(
        stable("registry-additive-refused"),
        [
            active_port(
                "additive-a",
                EulerChannel::RollingContourSpin,
                domain(
                    "disc",
                    "base",
                    "rim-patch",
                    0,
                    1,
                    0,
                    10,
                    "tangent-velocity",
                    "world",
                    PortOrientation::AlongFrame,
                ),
                ContributionOwnership::AdditiveWithProof {
                    decomposition_receipt: incomplete.clone(),
                },
            ),
            active_port(
                "additive-b",
                EulerChannel::RollingContourSpin,
                domain(
                    "disc",
                    "base",
                    "rim-patch",
                    0,
                    1,
                    0,
                    10,
                    "tangent-velocity",
                    "world",
                    PortOrientation::AlongFrame,
                ),
                ContributionOwnership::AdditiveWithProof {
                    decomposition_receipt: incomplete,
                },
            ),
        ],
    );
    assert!(matches!(
        refused,
        Err(EulerPortError::AdditiveProofMismatch { .. })
    ));

    let missing_contributor = DecompositionReceipt::try_new(
        stable("receipt-missing-contributor"),
        domain(
            "disc",
            "base",
            "rim-patch",
            0,
            1,
            0,
            10,
            "tangent-velocity",
            "world",
            PortOrientation::AlongFrame,
        ),
        [
            stable("complete-a"),
            stable("complete-b"),
            stable("absent-c"),
        ],
    )
    .expect("a receipt can structurally name three contributors");
    let absent = EulerPortRegistry::try_new(
        stable("registry-additive-absent"),
        [
            active_port(
                "complete-a",
                EulerChannel::RollingContourSpin,
                domain(
                    "disc",
                    "base",
                    "rim-patch",
                    0,
                    1,
                    0,
                    10,
                    "tangent-velocity",
                    "world",
                    PortOrientation::AlongFrame,
                ),
                ContributionOwnership::AdditiveWithProof {
                    decomposition_receipt: missing_contributor.clone(),
                },
            ),
            active_port(
                "complete-b",
                EulerChannel::RollingContourSpin,
                domain(
                    "disc",
                    "base",
                    "rim-patch",
                    0,
                    1,
                    0,
                    10,
                    "tangent-velocity",
                    "world",
                    PortOrientation::AlongFrame,
                ),
                ContributionOwnership::AdditiveWithProof {
                    decomposition_receipt: missing_contributor,
                },
            ),
        ],
    );
    assert!(matches!(
        absent,
        Err(EulerPortError::AdditiveProofMismatch { .. })
    ));
}

#[test]
fn ledger_is_exactly_once_and_unavailable_channels_remain_no_claim() {
    let active = single_active_registry();
    let unavailable = PortDeclaration::new(
        stable("gas-film-port"),
        EulerChannel::GasFilm,
        PortKind::MechanicalForceVelocity,
        ChannelActivity::Unavailable {
            model_identity: stable("gas-film-model"),
            reason_identity: stable("gas-film-capability-missing"),
        },
        stable("gas-film-law"),
        stable("gas-film-source"),
        domain(
            "disc",
            "gas",
            "gas-patch",
            0,
            1,
            0,
            10,
            "gas-velocity",
            "world",
            PortOrientation::AlongFrame,
        ),
        ContributionOwnership::Exclusive,
    );
    let registry = EulerPortRegistry::try_new(
        stable("registry-ledger"),
        active.declarations().iter().cloned().chain([unavailable]),
    )
    .expect("active and unavailable non-overlapping declarations");
    let mut ledger = EulerEnergyLedger::new(stable("ledger-main"), registry);
    let entry = EnergyContribution::new(
        stable("energy-entry-a"),
        stable("gravity-port"),
        EulerChannel::Gravity,
        timestamp(2),
        terms(3.0),
    );
    ledger.record(entry.clone()).expect("first receipt");
    assert_eq!(ledger.cumulative().kinetic_j(), 3.0);
    assert!(matches!(
        ledger.record(entry),
        Err(EulerPortError::DuplicateEnergyContribution { .. })
    ));
    assert_eq!(ledger.contributions().len(), 1);
    assert_eq!(
        ledger.closure_disposition(),
        EnergyClosureDisposition::NoClaimUnavailableChannels {
            channels: vec![EulerChannel::GasFilm],
        }
    );
    assert!(matches!(
        ledger.record(EnergyContribution::new(
            stable("unavailable-entry"),
            stable("gas-film-port"),
            EulerChannel::GasFilm,
            timestamp(2),
            terms(0.0),
        )),
        Err(EulerPortError::EnergyPortUnavailable { .. })
    ));
    assert_eq!(ledger.contributions().len(), 1);
}

#[test]
fn absent_and_inactive_channels_do_not_become_implicit_energy_sources() {
    let active_only =
        EulerEnergyLedger::new(stable("ledger-active-only"), single_active_registry());
    assert_eq!(
        active_only.closure_disposition(),
        EnergyClosureDisposition::NoClaimIntegrationSkeleton
    );

    let inactive = PortDeclaration::new(
        stable("base-port"),
        EulerChannel::Base,
        PortKind::MechanicalForceVelocity,
        ChannelActivity::Inactive {
            model_identity: stable("base-model-disabled"),
        },
        stable("base-law"),
        stable("base-source"),
        domain(
            "disc",
            "base",
            "base-patch",
            0,
            1,
            0,
            10,
            "base-velocity",
            "world",
            PortOrientation::AlongFrame,
        ),
        ContributionOwnership::Exclusive,
    );
    let registry = EulerPortRegistry::try_new(stable("registry-inactive"), [inactive])
        .expect("inactive declaration has no active ownership");
    let mut ledger = EulerEnergyLedger::new(stable("ledger-inactive"), registry);
    assert!(matches!(
        ledger.record(EnergyContribution::new(
            stable("inactive-entry"),
            stable("base-port"),
            EulerChannel::Base,
            timestamp(1),
            terms(0.0),
        )),
        Err(EulerPortError::EnergyPortInactive { .. })
    ));
    assert!(ledger.contributions().is_empty());
}

#[test]
fn ledger_refuses_permuted_receipts_and_checkpoint_rollback_is_prefix_bound() {
    let mut ledger = EulerEnergyLedger::new(stable("ledger-rollback"), single_active_registry());
    ledger
        .record(EnergyContribution::new(
            stable("entry-b"),
            stable("gravity-port"),
            EulerChannel::Gravity,
            timestamp(5),
            terms(1.0),
        ))
        .expect("first canonical contribution");
    let before_refusal = ledger.clone();
    assert!(matches!(
        ledger.record(EnergyContribution::new(
            stable("entry-a"),
            stable("gravity-port"),
            EulerChannel::Gravity,
            timestamp(4),
            terms(2.0),
        )),
        Err(EulerPortError::NonDeterministicEnergyContributionOrder)
    ));
    assert_eq!(ledger, before_refusal);

    let checkpoint = ledger.checkpoint(stable("checkpoint-one"), timestamp(5));
    ledger
        .record(EnergyContribution::new(
            stable("entry-c"),
            stable("gravity-port"),
            EulerChannel::Gravity,
            timestamp(6),
            terms(4.0),
        ))
        .expect("later canonical contribution");
    ledger.rollback(&checkpoint).expect("bound prefix rollback");
    assert_eq!(ledger.contributions().len(), 1);
    assert_eq!(ledger.cumulative().kinetic_j(), 1.0);

    let other = EulerEnergyLedger::new(stable("ledger-other"), single_active_registry());
    assert!(matches!(
        other.clone().rollback(&checkpoint),
        Err(EulerPortError::CheckpointIdentityMismatch)
    ));
}

#[test]
fn non_finite_and_negative_constrained_energy_terms_refuse() {
    assert!(matches!(
        EnergyTerms::try_new(f64::NAN, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0),
        Err(EulerPortError::NonFiniteEnergyTerms)
    ));
    assert!(matches!(
        EnergyTerms::try_new(0.0, 0.0, 0.0, -1.0, 0.0, 0.0, 0.0),
        Err(EulerPortError::NegativeLossOrUnresolvedEnergy)
    ));
}
