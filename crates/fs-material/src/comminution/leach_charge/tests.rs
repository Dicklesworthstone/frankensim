//! Authored arithmetic fixtures, not measured mineralogy, grain sizes, or ore data.

use super::*;
use crate::separation::ParticleClass;
use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, ObservationDataset, PropertyClaim, PropertyValue,
    UncertaintyModel,
};
use fs_qty::{Dims, MassFlowRate};

fn provenance() -> Provenance {
    Provenance {
        source: "authored leach-charge arithmetic fixture; not measured ore".into(),
        license: "project test fixture".into(),
        artifact: Some(hash_domain(
            "org.frankensim.test.leach-charge.v1",
            b"synthetic allocation",
        )),
    }
}
fn point(value: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "reference-context",
            QuantitySpec::dimensional(Dims::NONE),
            value,
        )
        .unwrap()
}
fn pack(property: &str, value: f64, dims: Dims) -> NormalizedPack {
    let provenance = provenance();
    let artifact = provenance.artifact.unwrap();
    let mut claims = ClaimSet::new();
    let observation = claims
        .register_observation(ObservationDataset {
            specimen: "synthetic mineral".into(),
            method: "authored arithmetic fixture".into(),
            artifact,
            caveats: "not a measured molar mass or identified mineral".into(),
            provenance: provenance.clone(),
        })
        .unwrap();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::with_quantity(property, QuantitySpec::dimensional(dims)),
            value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with_quantity(
                "reference-context",
                QuantitySpec::dimensional(Dims::NONE),
                1.0,
                2.0,
            ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![observation],
            provenance,
        })
        .unwrap();
    NormalizedPack::new(
        "synthetic-mineral",
        "authored-fixture-v1",
        artifact,
        "project test fixture",
        claims,
        vec![],
        vec![],
    )
    .unwrap()
}
fn molar_mass(value: f64) -> MolarMassSource {
    MolarMassSource::resolve(
        &pack("molar-mass", value, MolarMass::DIMS),
        &point(1.0),
        SelectionPolicy::SingleClaimOnly,
    )
    .unwrap()
}
fn element(id: &str) -> ElementId {
    ElementId::new(id).unwrap()
}
fn class(id: &str, diameter: f64, dry: f64, copper: f64, tin: f64) -> ParticleClass {
    ParticleClass::new(
        id,
        Length::new(diameter),
        MassFlowRate::new(dry),
        [
            (element("Cu"), MassFlowRate::new(copper)),
            (element("Sn"), MassFlowRate::new(tin)),
        ]
        .into(),
    )
    .unwrap()
}
fn stream() -> ParticleStream {
    ParticleStream::new(
        vec![
            class("fine", 40.0e-6, 1.0, 0.01, 0.02),
            class("coarse", 200.0e-6, 2.0, 0.10, 0.02),
        ],
        MassFlowRate::new(4.0),
    )
    .unwrap()
}
fn allocation(id: &str, class_id: &str, fraction: f64, mm: f64) -> MineralClassAllocation {
    MineralClassAllocation {
        population_id: id.into(),
        class_id: class_id.into(),
        solid: SpeciesId::new("declared-solid").unwrap(),
        mass_fraction: Dimensionless::new(fraction),
        fraction_source: provenance(),
        molar_mass: molar_mass(mm),
        radius: ReactiveRadius::LiberatedSphere {
            provenance: provenance(),
        },
    }
}
fn collect(allocations: Vec<MineralClassAllocation>) -> Result<CollectedLeachCharge> {
    CollectedLeachCharge::collect(&stream(), Time::new(10.0), allocations)
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 8.0 * f64::EPSILON * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn typed_molar_mass_queries_keep_real_receipts_and_refuse_wrong_sources() {
    let selected_pack = pack("molar-mass", 0.2, MolarMass::DIMS);
    let source = MolarMassSource::resolve(
        &selected_pack,
        &point(1.5),
        SelectionPolicy::SingleClaimOnly,
    )
    .unwrap();
    assert_eq!(source.molar_mass(), MolarMass::new(0.2));
    assert_eq!(source.pack_id(), selected_pack.pack_id());
    assert_eq!(source.pack_hash(), selected_pack.content_hash());
    assert_eq!(source.query(), &point(1.5));
    selected_pack
        .claims()
        .verify_receipt(&source.answer().receipt)
        .unwrap();
    assert_eq!(source.answer().evidence.value.value, 0.2);
    assert!(matches!(
        MolarMassSource::resolve(
            &selected_pack,
            &point(3.0),
            SelectionPolicy::SingleClaimOnly
        ),
        Err(LeachChargeError::Material(_))
    ));
    for bad in [
        pack("density", 0.2, MolarMass::DIMS),
        pack("molar-mass", 0.2, Dims::NONE),
    ] {
        assert!(matches!(
            MolarMassSource::resolve(&bad, &point(1.0), SelectionPolicy::SingleClaimOnly),
            Err(LeachChargeError::Material(_))
        ));
    }
    for value in [0.0, -0.2] {
        let bad = pack("molar-mass", value, MolarMass::DIMS);
        assert!(matches!(
            MolarMassSource::resolve(&bad, &point(1.0), SelectionPolicy::SingleClaimOnly),
            Err(LeachChargeError::Invalid(_))
        ));
    }
}

#[test]
fn collection_retains_class_species_radius_sources_and_independent_assay_accounts() {
    let b = allocation("b", "coarse", 0.25, 0.2);
    let mut a = allocation("a", "coarse", 0.125, 0.1);
    a.radius = ReactiveRadius::Sourced {
        radius: Length::new(2.0e-6),
        provenance: provenance(),
    };
    let c = allocation("c", "fine", 0.5, 0.2);
    let original_allocations = vec![b.clone(), c.clone(), a.clone()];
    let charge = collect(original_allocations.clone()).unwrap();
    assert_eq!(charge.original_stream(), &stream());
    assert_eq!(charge.collection_duration(), Time::new(10.0));
    assert_eq!(
        charge
            .populations()
            .iter()
            .map(CollectedMineralPopulation::id)
            .collect::<Vec<_>>(),
        ["a", "b", "c"]
    );
    assert_eq!(charge.populations()[0].allocation(), &a);
    assert_eq!(charge.populations()[1].allocation(), &b);
    assert_eq!(charge.populations()[2].allocation(), &c);
    for population in charge.populations() {
        assert_eq!(
            population.solid(),
            &SpeciesId::new("declared-solid").unwrap()
        );
        close(population.amount().value(), 25.0);
        assert_eq!(population.mass_basis().mass(), population.mass());
    }
    close(charge.populations()[0].radius().value(), 2.0e-6);
    close(charge.populations()[1].radius().value(), 100.0e-6);
    close(charge.populations()[2].radius().value(), 20.0e-6);
    close(charge.dry_mass().value(), 30.0);
    close(charge.water_mass().value(), 40.0);
    close(charge.unallocated_mass().value(), 17.5);
    close(charge.elements()[&element("Cu")].value(), 1.1);
    close(charge.elements()[&element("Sn")].value(), 0.4);
    let coarse = &charge.class_accounts()[0];
    assert_eq!(coarse.class_id(), "coarse");
    close(coarse.dry_mass().value(), 20.0);
    close(coarse.allocated_mass().value(), 7.5);
    close(coarse.unallocated_mass().value(), 12.5);
    close(coarse.elements()[&element("Cu")].value(), 1.0);
    for class in charge.class_accounts() {
        close(
            class.allocated_mass().value() + class.unallocated_mass().value(),
            class.dry_mass().value(),
        );
    }
    // Collection time changes inventory, never the sourced reactive radius.
    let double =
        CollectedLeachCharge::collect(&stream(), Time::new(20.0), original_allocations).unwrap();
    for (first, second) in charge.populations().iter().zip(double.populations()) {
        close(second.amount().value(), 2.0 * first.amount().value());
        assert_eq!(first.radius(), second.radius());
    }
    close(double.elements()[&element("Cu")].value(), 2.2);
    close(double.water_mass().value(), 80.0);
}

#[test]
fn zero_collection_zero_flows_and_unallocated_mass_remain_explicit() {
    let allocations = vec![
        allocation("zero", "coarse", 0.0, 0.2),
        allocation("nonzero", "fine", 0.5, 0.2),
    ];
    let zero =
        CollectedLeachCharge::collect(&stream(), Time::new(0.0), allocations.clone()).unwrap();
    assert_eq!(zero.populations().len(), 2);
    assert!(zero.populations().iter().all(|p| p.amount().value() == 0.0));
    assert_eq!(zero.dry_mass().value(), 0.0);
    assert_eq!(zero.water_mass().value(), 0.0);
    assert!(zero.elements().values().all(|v| v.value() == 0.0));
    let unallocated = collect(vec![]).unwrap();
    assert!(unallocated.populations().is_empty());
    assert_eq!(unallocated.unallocated_mass(), unallocated.dry_mass());
    assert_eq!(
        unallocated.elements(),
        collect(allocations).unwrap().elements()
    );
    let empty_stream = ParticleStream::new(
        vec![class("coarse", 1.0e-4, 0.0, 0.0, 0.0)],
        MassFlowRate::new(0.0),
    )
    .unwrap();
    let empty = CollectedLeachCharge::collect(
        &empty_stream,
        Time::new(10.0),
        vec![allocation("p", "coarse", 1.0, 0.2)],
    )
    .unwrap();
    assert_eq!(empty.populations()[0].amount().value(), 0.0);
    assert_eq!(empty.class_accounts().len(), 1);
    assert_eq!(empty.class_accounts()[0].elements().len(), 2);
}

#[test]
fn identity_class_fraction_and_collection_bounds_are_enforced() {
    let a = allocation("same", "coarse", 0.5, 0.2);
    assert_eq!(
        collect(vec![a.clone(), a]),
        Err(LeachChargeError::DuplicatePopulation("same".into()))
    );
    assert_eq!(
        collect(vec![allocation("p", "absent", 0.5, 0.2)]),
        Err(LeachChargeError::UnknownClass("absent".into()))
    );
    for id in ["", "  ", "invalid\n", &"x".repeat(257)] {
        assert!(matches!(
            collect(vec![allocation(id, "coarse", 0.5, 0.2)]),
            Err(LeachChargeError::Invalid(_))
        ));
    }
    for fraction in [-0.1, 1.1, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            collect(vec![allocation("p", "coarse", fraction, 0.2)]),
            Err(LeachChargeError::Invalid(_))
        ));
    }
    assert_eq!(
        collect(vec![
            allocation("a", "coarse", 0.75, 0.2),
            allocation("b", "coarse", 0.5, 0.2)
        ]),
        Err(LeachChargeError::OverallocatedClass("coarse".into()))
    );
    for duration in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            CollectedLeachCharge::collect(&stream(), Time::new(duration), vec![]),
            Err(LeachChargeError::Invalid(_))
        ));
    }
    let many = vec![allocation("p", "coarse", 0.0, 0.2); MAX_MINERAL_POPULATIONS + 1];
    assert!(matches!(collect(many), Err(LeachChargeError::Invalid(_))));
    let classes = (0..=MAX_COLLECTED_CLASSES)
        .map(|i| class(&format!("class-{i}"), 1.0e-4, 0.0, 0.0, 0.0))
        .collect();
    let many = ParticleStream::new(classes, MassFlowRate::new(0.0)).unwrap();
    assert!(matches!(
        CollectedLeachCharge::collect(&many, Time::new(1.0), vec![]),
        Err(LeachChargeError::Invalid(_))
    ));
}

#[test]
fn source_evidence_and_radius_are_required_even_for_zero_inventory() {
    let mut a = allocation("p", "coarse", 0.0, 0.2);
    a.fraction_source.artifact = None;
    assert!(matches!(
        collect(vec![a]),
        Err(LeachChargeError::Invalid(_))
    ));
    for invalid in [
        Provenance {
            source: " ".into(),
            ..provenance()
        },
        Provenance {
            license: "".into(),
            ..provenance()
        },
        Provenance {
            artifact: None,
            ..provenance()
        },
    ] {
        let mut a = allocation("p", "coarse", 0.0, 0.2);
        a.radius = ReactiveRadius::LiberatedSphere {
            provenance: invalid.clone(),
        };
        assert!(matches!(
            collect(vec![a.clone()]),
            Err(LeachChargeError::Invalid(_))
        ));
        a.radius = ReactiveRadius::Sourced {
            radius: Length::new(1.0e-6),
            provenance: invalid,
        };
        assert!(matches!(
            collect(vec![a]),
            Err(LeachChargeError::Invalid(_))
        ));
    }
    for radius in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let mut a = allocation("p", "coarse", 0.0, 0.2);
        a.radius = ReactiveRadius::Sourced {
            radius: Length::new(radius),
            provenance: provenance(),
        };
        assert!(matches!(
            collect(vec![a]),
            Err(LeachChargeError::Invalid(_))
        ));
    }
}

#[test]
fn exact_fraction_budget_and_numeric_range_preserve_small_accounts_or_refuse() {
    let tiny = f64::from_bits(1);
    // Ordinary 1 + tiny == 1 must not authorize an extra positive allocation.
    assert_eq!(
        collect(vec![
            allocation("a", "coarse", 1.0, 0.2),
            allocation("b", "coarse", tiny, 0.2)
        ]),
        Err(LeachChargeError::OverallocatedClass("coarse".into()))
    );
    let below_half = f64::from_bits(0.5_f64.to_bits() - 1);
    let charge = collect(vec![
        allocation("a", "coarse", 0.5, 0.2),
        allocation("b", "coarse", below_half, 0.2),
    ])
    .unwrap();
    assert_eq!(
        charge.class_accounts()[0].unallocated_mass().value(),
        20.0 * (0.5 - below_half)
    );
    assert!(charge.class_accounts()[0].unallocated_mass().value() > 0.0);
    let one = |dry, water, diameter| {
        ParticleStream::new(
            vec![class("coarse", diameter, dry, 0.0, 0.0)],
            MassFlowRate::new(water),
        )
        .unwrap()
    };
    for (feed, duration, allocation) in [
        (one(f64::MAX, 0.0, 1.0), 2.0, vec![]),
        (one(1.0, f64::MAX, 1.0), 2.0, vec![]),
        (one(tiny, 0.0, 1.0), 0.5, vec![]),
        (
            one(tiny, 0.0, 1.0),
            1.0,
            vec![allocation("p", "coarse", 0.5, 1.0)],
        ),
        (
            one(tiny, 0.0, 1.0),
            1.0,
            vec![allocation("p", "coarse", 1.0, f64::MAX)],
        ),
        (
            one(0.0, 0.0, tiny),
            1.0,
            vec![allocation("p", "coarse", 0.0, 0.2)],
        ),
    ] {
        assert!(matches!(
            CollectedLeachCharge::collect(&feed, Time::new(duration), allocation),
            Err(LeachChargeError::Arithmetic(_))
        ));
    }
    assert!(matches!(
        CollectedLeachCharge::collect(
            &one(1.0, 0.0, 1.0),
            Time::new(1.0),
            vec![allocation("p", "coarse", 1.0, tiny)]
        ),
        Err(LeachChargeError::Chemistry(_))
    ));
    let trace = ParticleStream::new(
        vec![class("coarse", 1.0, 1.0, tiny, 0.0)],
        MassFlowRate::new(0.0),
    )
    .unwrap();
    assert!(matches!(
        CollectedLeachCharge::collect(&trace, Time::new(0.5), vec![]),
        Err(LeachChargeError::Arithmetic(_))
    ));
    let overflow = ParticleStream::new(
        vec![
            class("a", 1.0, f64::MAX / 4.0, 0.0, 0.0),
            class("b", 1.0, f64::MAX / 4.0, 0.0, 0.0),
        ],
        MassFlowRate::new(0.0),
    )
    .unwrap();
    assert!(matches!(
        CollectedLeachCharge::collect(&overflow, Time::new(3.0), vec![]),
        Err(LeachChargeError::Arithmetic(_))
    ));
}
