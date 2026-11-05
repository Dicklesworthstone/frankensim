use super::super::kinetics::UniformizationBudget;
use super::*;
use crate::separation::{ParticleClass, ParticleStream};
use fs_blake3::hash_bytes;
use fs_matdb::{ClaimSet, Provenance, UncertaintyModel};
use fs_qty::{MassFlowRate, chemistry::ElementId};
use std::collections::BTreeMap;

fn conditions() -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "test-condition",
            QuantitySpec::dimensional(Dimensionless::DIMS),
            1.0,
        )
        .unwrap()
}
fn fixture(edit: impl FnOnce(&mut ClaimSet, &mut Vec<PropertyClaim>)) -> NormalizedPack {
    let artifact = hash_bytes(b"authored Austin fit fixture; not measured ore");
    let provenance = Provenance {
        source: "authored analytical fixture".into(),
        license: "project fixture".into(),
        artifact: Some(artifact),
    };
    let mut set = ClaimSet::new();
    let observation = set
        .register_observation(ObservationDataset {
            specimen: "synthetic three-class source".into(),
            method: "authored exact chain".into(),
            artifact,
            caveats: "synthetic; no experimental calibration".into(),
            provenance: provenance.clone(),
        })
        .unwrap();
    let validity = ValidityDomain::unconstrained()
        .with_quantity(
            "test-condition",
            QuantitySpec::dimensional(Dimensionless::DIMS),
            1.0,
            1.0,
        )
        .with_quantity(
            CLASS_INDEX_AXIS,
            QuantitySpec::dimensional(Dimensionless::DIMS),
            0.0,
            2.0,
        )
        .with_quantity(
            CONTACT_TIME_AXIS,
            QuantitySpec::dimensional(Time::DIMS),
            0.0,
            2.0,
        );
    let mut claims: Vec<_> = parameter_keys()
        .into_iter()
        .zip([1.0, 1.0, 0.001, 0.0, 1.0, 1.0, 1.0])
        .map(|(key, value)| PropertyClaim {
            value: PropertyValue::Scalar {
                value,
                dims: key.dims(),
            },
            key,
            validity: validity.clone(),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![observation],
            provenance: provenance.clone(),
        })
        .collect();
    claims.push(PropertyClaim {
        key: class_size_key(),
        value: PropertyValue::Curve {
            abscissa: CLASS_INDEX_AXIS.into(),
            abscissa_dims: Dimensionless::DIMS,
            knots: vec![(0.0, 0.004), (1.0, 0.002), (2.0, 0.001)],
            dims: Length::DIMS,
        },
        validity,
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::TabulatedOnly,
        observations: vec![observation],
        provenance,
    });
    edit(&mut set, &mut claims);
    for claim in claims {
        set.insert_claim(claim).unwrap();
    }
    NormalizedPack::new(
        "synthetic-austin",
        "unit-test",
        artifact,
        "project fixture",
        set,
        vec![],
        vec![],
    )
    .unwrap()
}
fn resolve(pack: &NormalizedPack) -> Result<SelectionBreakageModel> {
    resolve_austin_grinding(
        pack,
        &conditions(),
        Time::new(2.0),
        SelectionPolicy::SingleClaimOnly,
    )
}
fn budget() -> UniformizationBudget {
    UniformizationBudget {
        maximum_l1_tail: Dimensionless::new(1e-13),
        max_poisson_terms: 256,
        max_coefficient_products: 100_000,
    }
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-12 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn exact_daughter_indexing_and_sink_drive_independent_batch_solution() {
    // S_0=2/s, S_1=1/s. B_10=1, B_20=1/2, hence b_10=b_20=1/2.
    // Initial top-class mass 1 gives m_0=e^-2t, m_1=e^-t-e^-2t,
    // m_sink=1-e^-t. A mistaken x_i/x_j offset loses half the parent mass.
    let pack = fixture(|_, _| {});
    let model = resolve(&pack).unwrap();
    let daughters = model.daughters();
    close(
        daughters.fraction("size-001", "size-000").unwrap().value(),
        0.5,
    );
    close(
        daughters
            .fraction("fines-sink", "size-000")
            .unwrap()
            .value(),
        0.5,
    );
    close(
        daughters
            .fraction("fines-sink", "size-001")
            .unwrap()
            .value(),
        1.0,
    );
    for (class, rate) in daughters.classes().iter().zip(model.rates()) {
        assert_eq!(
            rate.value(),
            match class.id.as_str() {
                "size-000" => 2.0,
                "size-001" => 1.0,
                "fines-sink" => 0.0,
                _ => unreachable!(),
            }
        );
    }
    let si = ElementId::new("Si").unwrap();
    let feed = ParticleStream::new(
        daughters
            .classes()
            .iter()
            .map(|class| {
                let flow = if class.id == "size-000" { 1.0 } else { 0.0 };
                ParticleClass::new(
                    class.id.clone(),
                    class.diameter,
                    MassFlowRate::new(flow),
                    BTreeMap::from([(si.clone(), MassFlowRate::new(0.4 * flow))]),
                )
            })
            .collect::<core::result::Result<_, _>>()
            .unwrap(),
        MassFlowRate::new(0.0),
    )
    .unwrap();
    for t in [0.0, 0.5, 1.0, 2.0] {
        let pass = model
            .prepare(Time::new(t), &conditions(), budget(), || false)
            .unwrap();
        let output = pass.transform(&feed).unwrap();
        for class in output.classes() {
            let expected = match class.id() {
                "size-000" => (-2.0 * t).exp(),
                "size-001" => (-t).exp() - (-2.0 * t).exp(),
                "fines-sink" => 1.0 - (-t).exp(),
                _ => unreachable!(),
            };
            close(class.dry_solids().value(), expected);
            close(class.elements()[&si].value(), 0.4 * expected);
        }
        close(output.dry_solids().unwrap().value(), 1.0);
        close(output.element_flow(&si).unwrap().value(), 0.4);
        assert_eq!(output.water().value(), 0.0);
        assert_eq!(
            pass.model().austin_calibration_source(),
            model.austin_calibration_source()
        );
    }
    let source = model.austin_calibration_source().unwrap();
    assert_eq!(source.pack_id(), pack.pack_id());
    assert_eq!(source.pack_hash(), pack.content_hash());
    assert_eq!(source.conditions(), &conditions());
    assert_eq!(source.maximum_contact_time(), Time::new(2.0));
    assert_eq!(source.claims().len(), 8);
    assert_eq!(source.answers().len(), 10);
    for envelope in source.answers() {
        for answer in [&envelope.lower, &envelope.upper] {
            pack.claims().verify_receipt(&answer.receipt).unwrap();
            assert_eq!(
                answer.evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
        }
    }
}

#[test]
fn original_fit_grid_and_complete_time_support_cannot_be_mixed_or_widened() {
    let other = fixture(|set, claims| {
        let mut observation = set.observation(claims[0].observations[0]).unwrap().clone();
        observation.specimen = "different experiment".into();
        claims[4].observations = vec![set.register_observation(observation).unwrap()];
    });
    assert!(resolve(&other).is_err());
    let support = fixture(|_, claims| {
        claims[0].validity = claims[0].validity.clone().with_quantity(
            CONTACT_TIME_AXIS,
            QuantitySpec::dimensional(Time::DIMS),
            0.0,
            3.0,
        );
    });
    assert!(resolve(&support).is_err());
    let missing = fixture(|_, claims| {
        claims[6].key = PropertyKey::new("unrelated-coefficient", Dimensionless::DIMS);
    });
    assert!(resolve(&missing).is_err());
    let pack = fixture(|_, _| {});
    assert!(
        resolve_austin_grinding(
            &pack,
            &conditions(),
            Time::new(2.01),
            SelectionPolicy::SingleClaimOnly
        )
        .is_err()
    );
    for t in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            resolve_austin_grinding(
                &pack,
                &conditions(),
                Time::new(t),
                SelectionPolicy::SingleClaimOnly
            )
            .is_err()
        );
    }
    let model = resolve(&pack).unwrap();
    assert!(
        model
            .prepare(Time::new(2.01), &conditions(), budget(), || false)
            .is_err()
    );
    let extra = conditions()
        .with_quantity(
            "temperature",
            QuantitySpec::dimensional(fs_qty::Temperature::DIMS),
            298.15,
        )
        .unwrap();
    assert!(
        resolve_austin_grinding(
            &pack,
            &extra,
            Time::new(1.0),
            SelectionPolicy::SingleClaimOnly
        )
        .is_err()
    );
    assert!(
        model
            .prepare(Time::new(1.0), &extra, budget(), || false)
            .is_err()
    );
    let changed = conditions().with("test-condition", 2.0).unwrap();
    assert!(
        model
            .prepare(Time::new(1.0), &changed, budget(), || false)
            .is_err()
    );
    let untyped = QueryPoint::new().with("test-condition", 1.0).unwrap();
    assert!(
        resolve_austin_grinding(
            &pack,
            &untyped,
            Time::new(1.0),
            SelectionPolicy::SingleClaimOnly
        )
        .is_err()
    );
}

#[test]
fn unsupported_grids_and_nonphysical_coefficients_refuse_without_repair() {
    let interpolated = fixture(|_, claims| {
        claims[7].interpolation = InterpolationPolicy::LinearInside;
    });
    assert!(resolve(&interpolated).is_err());
    let missing_index = fixture(|_, claims| {
        if let PropertyValue::Curve { knots, .. } = &mut claims[7].value {
            knots[1].0 = 1.5;
        }
    });
    assert!(resolve(&missing_index).is_err());
    let equal_sizes = fixture(|_, claims| {
        if let PropertyValue::Curve { knots, .. } = &mut claims[7].value {
            knots[1].1 = knots[0].1;
        }
    });
    assert!(resolve(&equal_sizes).is_err());
    for (index, value) in [
        (0, 0.0),
        (1, 0.0),
        (2, -1.0),
        (3, -1.0),
        (4, 1.01),
        (5, 0.0),
        (6, -1.0),
    ] {
        let pack = fixture(|_, claims| {
            if let PropertyValue::Scalar { value: stored, .. } = &mut claims[index].value {
                *stored = value;
            }
        });
        assert!(resolve(&pack).is_err());
    }
    let vanished = fixture(|_, claims| {
        for index in [5, 6] {
            if let PropertyValue::Scalar { value, .. } = &mut claims[index].value {
                *value = f64::MAX;
            }
        }
    });
    assert!(resolve(&vanished).is_err());
}
