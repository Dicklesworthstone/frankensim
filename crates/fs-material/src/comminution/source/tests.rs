use super::*;
use fs_blake3::hash_bytes;
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, Provenance, UncertaintyModel};
use fs_qty::{Dimensionless, Dims, Power, QuantitySpec, Temperature};

fn point(sample: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("sample", QuantitySpec::dimensional(Dims::NONE), sample)
        .unwrap()
}

fn domain() -> BondSizeDomain {
    BondSizeDomain::new([Length::new(0.0025); 2], [Length::new(0.0001); 2]).unwrap()
}

fn fixture(
    basis: WorkIndexMassBasis,
    edit: impl FnOnce(&mut ClaimSet, &mut [PropertyClaim; 4]),
) -> NormalizedPack {
    let provenance = Provenance {
        source: "synthetic source fixture; no measured ore".into(),
        license: "project test fixture".into(),
        artifact: Some(hash_bytes(b"synthetic Bond row")),
    };
    let mut set = ClaimSet::new();
    let observation = set
        .register_observation(ObservationDataset {
            specimen: "authored composite 1".into(),
            method: "authored Bond ball-mill test fixture".into(),
            artifact: hash_bytes(b"synthetic Bond row"),
            caveats: "unknown temperature and preparation; no site calibration".into(),
            provenance: provenance.clone(),
        })
        .unwrap();
    let mut claims = [
        (work_index_key(basis), 54_000.0),
        (
            PropertyKey::new("bond-test-closing-screen", Length::DIMS),
            150e-6,
        ),
        (
            PropertyKey::new("bond-test-feed-f80", Length::DIMS),
            2500e-6,
        ),
        (
            PropertyKey::new("bond-test-product-p80", Length::DIMS),
            100e-6,
        ),
    ]
    .map(|(key, value)| PropertyClaim {
        value: PropertyValue::Scalar {
            value,
            dims: key.dims(),
        },
        key,
        validity: ValidityDomain::unconstrained().with_quantity(
            "sample",
            QuantitySpec::dimensional(Dims::NONE),
            1.0,
            1.0,
        ),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![observation],
        provenance: provenance.clone(),
    });
    edit(&mut set, &mut claims);
    for claim in claims {
        set.insert_claim(claim).unwrap();
    }
    NormalizedPack::new(
        "synthetic-bond-test",
        "unit-test",
        hash_bytes(b"synthetic Bond row"),
        "project fixture",
        set,
        vec![],
        vec![],
    )
    .unwrap()
}

fn resolve(pack: &NormalizedPack) -> Result<BondGrindingCard, BondSourceError> {
    resolve_bond_ball_mill(
        pack,
        &point(1.0),
        WorkIndexMassBasis::MetricTonne,
        SelectionPolicy::SingleClaimOnly,
        domain(),
    )
}

fn close(actual: f64, expected: f64) {
    assert!((actual - expected).abs() <= expected.abs().max(1.0) * 4e-15);
}

#[test]
fn original_source_survives_energy_throughput_and_card_clone() {
    let pack = fixture(WorkIndexMassBasis::MetricTonne, |_, _| {});
    let card = resolve(&pack).unwrap();
    let source = card.calibration_source().unwrap();
    assert_eq!(source.pack_id(), pack.pack_id());
    assert_eq!(source.pack_hash(), pack.content_hash());
    assert_eq!(source.query(), &point(1.0));
    assert_eq!(source.observation().specimen, card.test_context().specimen);
    assert_eq!(source.observation().method, card.test_context().method);
    assert_eq!(
        source.observation().provenance,
        card.test_context().provenance
    );
    assert!(source.observation().caveats.contains("unknown temperature"));
    for (claim, answer) in source.claims().iter().zip(source.answers()) {
        assert_eq!(pack.claims().claim(answer.receipt.selected), Some(claim));
        pack.claims().verify_receipt(&answer.receipt).unwrap();
        assert_eq!(
            answer.evidence.value.uncertainty,
            UncertaintyModel::Unstated
        );
    }
    let [f80, p80] = source.test_sizes();
    close(card.work_index().source_value(), 15.0);
    close(card.specific_energy(f80, p80).unwrap().value(), 43_200.0);
    close(
        card.dry_throughput(f80, p80, Power::new(12e6), Dimensionless::new(0.9))
            .unwrap()
            .value(),
        250.0,
    );
    assert_eq!(card.clone().calibration_source(), Some(source));
}

#[test]
fn original_mass_convention_is_selected_not_relabelled() {
    let metric = fixture(WorkIndexMassBasis::MetricTonne, |_, _| {});
    assert!(
        resolve_bond_ball_mill(
            &metric,
            &point(1.0),
            WorkIndexMassBasis::ShortTon,
            SelectionPolicy::SingleClaimOnly,
            domain(),
        )
        .is_err()
    );
    let short = fixture(WorkIndexMassBasis::ShortTon, |_, _| {});
    let card = resolve_bond_ball_mill(
        &short,
        &point(1.0),
        WorkIndexMassBasis::ShortTon,
        SelectionPolicy::SingleClaimOnly,
        domain(),
    )
    .unwrap();
    assert_eq!(
        card.work_index().source_basis(),
        WorkIndexMassBasis::ShortTon
    );
    close(card.work_index().source_value(), 15.0 * 907.18474 / 1000.0);
    assert_eq!(card.work_index().specific_energy().value(), 54_000.0);
}

#[test]
fn reported_sizes_do_not_widen_the_callers_operational_domain() {
    let pack = fixture(WorkIndexMassBasis::MetricTonne, |_, _| {});
    let operational =
        BondSizeDomain::new([Length::new(0.003); 2], [Length::new(0.00012); 2]).unwrap();
    let card = resolve_bond_ball_mill(
        &pack,
        &point(1.0),
        WorkIndexMassBasis::MetricTonne,
        SelectionPolicy::SingleClaimOnly,
        operational,
    )
    .unwrap();
    assert_eq!(card.domain(), operational);
    let [reported_f, reported_p] = card.calibration_source().unwrap().test_sizes();
    assert!(card.specific_energy(reported_f, reported_p).is_err());
    assert!(
        card.specific_energy(Length::new(0.003), Length::new(0.00012))
            .is_ok()
    );
    close(reported_p.value(), 0.0001);
}

#[test]
fn unknown_or_extra_source_conditions_refuse() {
    let pack = fixture(WorkIndexMassBasis::MetricTonne, |_, _| {});
    for query in [
        QueryPoint::new(),
        point(2.0),
        point(1.0)
            .with_quantity(
                "temperature",
                QuantitySpec::dimensional(Temperature::DIMS),
                298.15,
            )
            .unwrap(),
    ] {
        assert!(
            resolve_bond_ball_mill(
                &pack,
                &query,
                WorkIndexMassBasis::MetricTonne,
                SelectionPolicy::SingleClaimOnly,
                domain(),
            )
            .is_err()
        );
    }
}

#[test]
fn mixed_specimens_and_citation_only_tests_refuse() {
    let mixed = fixture(WorkIndexMassBasis::MetricTonne, |set, claims| {
        let original = set.observation(claims[0].observations[0]).unwrap();
        let mut other = original.clone();
        other.specimen = "different composite".into();
        claims[2].observations = vec![set.register_observation(other).unwrap()];
    });
    assert!(matches!(resolve(&mixed), Err(BondSourceError::Invalid(_))));
    let uncited = fixture(WorkIndexMassBasis::MetricTonne, |_, claims| {
        claims[0].observations.clear();
    });
    assert!(matches!(
        resolve(&uncited),
        Err(BondSourceError::Invalid(_))
    ));
}

#[test]
fn wrong_units_nonphysical_rows_and_invented_interpolation_refuse() {
    for (field, value) in [(0, 0.0), (2, 90e-6), (3, 160e-6)] {
        let pack = fixture(WorkIndexMassBasis::MetricTonne, |_, claims| {
            claims[field].value = PropertyValue::Scalar {
                value,
                dims: claims[field].key.dims(),
            };
        });
        assert!(resolve(&pack).is_err());
    }
    let wrong = fixture(WorkIndexMassBasis::MetricTonne, |_, claims| {
        claims[0].key = PropertyKey::new("bond-ball-mill-work-index-metric-tonne", Dims::NONE);
        claims[0].value = PropertyValue::Scalar {
            value: 15.0,
            dims: Dims::NONE,
        };
    });
    assert!(matches!(resolve(&wrong), Err(BondSourceError::Material(_))));
    let interpolated = fixture(WorkIndexMassBasis::MetricTonne, |_, claims| {
        claims[0].value = PropertyValue::Curve {
            abscissa: "sample".into(),
            abscissa_dims: Dims::NONE,
            dims: SpecificGrindingEnergy::DIMS,
            knots: vec![(1.0, 54_000.0), (2.0, 55_000.0)],
        };
        claims[0].interpolation = InterpolationPolicy::LinearInside;
    });
    assert!(matches!(
        resolve(&interpolated),
        Err(BondSourceError::Invalid(_))
    ));
}
