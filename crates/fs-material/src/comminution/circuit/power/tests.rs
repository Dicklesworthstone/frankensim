use super::*;
use crate::comminution::breakage::{
    BreakageColumn, BreakageError, BreakageTestContext, BreakageTransfer, SizeClass,
    SizeRedistributionCard,
};
use crate::comminution::circuit::MILL_DRY_FLOW_AXIS;
use crate::comminution::{BondSizeDomain, BondTestContext, BondWorkIndex};
use crate::separation::{
    CLASSIFIER_LAW_ID, ClassifierCard, PARTICLE_DIAMETER_AXIS, ParticleClass, SOLIDS_FRACTION_AXIS,
};
use fs_blake3::hash_bytes;
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::{Dims, QuantitySpec, chemistry::ElementId};
use std::collections::BTreeMap;

fn source() -> Provenance {
    Provenance {
        source: "authored Bond/circuit arithmetic fixture; not measured ore".into(),
        license: "project fixture".into(),
        artifact: Some(hash_bytes(b"bond-circuit-power-v1")),
    }
}
fn copper() -> ElementId {
    ElementId::new("Cu").unwrap()
}
fn feed() -> ParticleStream {
    ParticleStream::new(
        [("coarse", 200e-6, 10.0, 2.0), ("fine", 100e-6, 2.0, 0.2)]
            .into_iter()
            .map(|(id, d, mass, cu)| {
                ParticleClass::new(
                    id,
                    Length::new(d),
                    MassFlowRate::new(mass),
                    BTreeMap::from([(copper(), MassFlowRate::new(cu))]),
                )
                .unwrap()
            })
            .collect(),
        MassFlowRate::new(12.0),
    )
    .unwrap()
}
fn circuit(bounded_flow: bool) -> ClosedGrindingCircuit {
    let mut validity = ValidityDomain::unconstrained();
    if bounded_flow {
        validity = validity.with_quantity(
            MILL_DRY_FLOW_AXIS,
            QuantitySpec::dimensional(MassFlowRate::DIMS),
            38.0,
            39.0,
        );
    }
    let mill = SizeRedistributionCard::new(
        [("coarse", 200e-6), ("fine", 100e-6)]
            .into_iter()
            .map(|(id, d)| SizeClass {
                id: id.into(),
                diameter: Length::new(d),
            })
            .collect(),
        [
            ("coarse", vec![("fine", 0.5), ("coarse", 0.5)]),
            ("fine", vec![("fine", 1.0)]),
        ]
        .into_iter()
        .map(|(input, transfers)| BreakageColumn {
            input: input.into(),
            transfers: transfers
                .into_iter()
                .map(|(output, f)| BreakageTransfer {
                    output: output.into(),
                    mass_fraction: Dimensionless::new(f),
                })
                .collect(),
        })
        .collect(),
        BreakageTestContext {
            specimen: "authored composite".into(),
            method: "explicit congruent half-coarse breakage per pass".into(),
            provenance: source(),
        },
        validity,
    )
    .unwrap();
    let classifier = ClassifierCard::from_model_card(&ConstitutiveModelCard {
        law: LawId(CLASSIFIER_LAW_ID.into()),
        law_version: 1,
        state_schema_version: 1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        parameters: [
            ("corrected-cut-size", 100e-6, Length::DIMS),
            ("sharpness", 2.0, Dims::NONE),
            ("water-recovery-to-underflow", 0.2, Dims::NONE),
        ]
        .into_iter()
        .map(|(key, value, dims)| (key.into(), LawParameter { value, dims }))
        .collect(),
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                PARTICLE_DIAMETER_AXIS,
                QuantitySpec::dimensional(Length::DIMS),
                70e-6,
                1e-3,
            )
            .with_quantity(
                SOLIDS_FRACTION_AXIS,
                QuantitySpec::dimensional(Dims::NONE),
                0.0,
                1.0,
            ),
        sources: vec![source().artifact.unwrap()],
        provenance: source(),
    })
    .unwrap();
    ClosedGrindingCircuit::new(mill, classifier)
}
fn bond() -> BondGrindingCard {
    BondGrindingCard::new(
        BondWorkIndex::from_kwh_per_metric_tonne(15.0).unwrap(),
        BondTestContext {
            specimen: "authored composite".into(),
            method: "explicit authored Bond fixture".into(),
            closing_screen: Length::new(125e-6),
            provenance: source(),
        },
        BondSizeDomain::new(
            [Length::new(100e-6), Length::new(1e-3)],
            [Length::new(70e-6), Length::new(500e-6)],
        )
        .unwrap(),
    )
    .unwrap()
}
fn settings() -> CircuitSettings {
    CircuitSettings {
        absolute_balance_tolerance: MassFlowRate::new(1e-11),
        relative_balance_tolerance: Dimensionless::new(1e-12),
        maximum_circulating_load: Dimensionless::new(100.0),
    }
}
fn prototype(bounded_flow: bool, split: Partition) -> Result<BondCircuitPrototype> {
    BondCircuitPrototype::new(
        &circuit(bounded_flow),
        &bond(),
        &feed(),
        &QueryPoint::new(),
        &QueryPoint::new(),
        settings(),
        split,
    )
}
fn return_all() -> Partition {
    Partition::from_product(Dimensionless::new(1.0)).unwrap()
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= expected.abs().max(1e-300) * 3e-12,
        "{actual} != {expected}",
    );
}

#[test]
fn independent_bond_power_uses_fresh_tonnage_and_scales_every_account() {
    let prototype = prototype(false, return_all()).unwrap();
    let preserved = prototype.clone();
    // Wi=15 kWh/t, F80=200 um, P80=100 um; the independent Bond equation
    // is 10*15*(1/sqrt(100)-1/sqrt(200))*3600 J/kg.
    let energy = 150.0 * (0.1 - 1.0 / 200.0_f64.sqrt()) * 3600.0;
    let point = prototype
        .at_net_pinion_power(Power::new(energy * 24.0))
        .unwrap();
    close(
        point.solution.fresh_feed.dry_solids().unwrap().value(),
        24.0,
    );
    close(point.specific_energy.value(), energy);
    close(point.solution.fresh_feed.water().value(), 24.0);
    close(
        point
            .solution
            .product
            .element_flow(&copper())
            .unwrap()
            .value(),
        4.4,
    );
    // Independently: coarse mill feed=10/(1-.95*.5), fine feed=
    // (2+.6*.5*coarse)/.4; doubling gives total230/3, not24 kg/s.
    close(
        point.solution.mill_feed.dry_solids().unwrap().value(),
        230.0 / 3.0,
    );
    close(point.solution.circulating_load.value(), 79.0 / 36.0);
    for (scaled, reference) in point
        .solution
        .fresh_feed
        .classes()
        .iter()
        .zip(prototype.reference().fresh_feed.classes())
    {
        assert_eq!(scaled.id(), reference.id());
        close(
            scaled.dry_solids().value(),
            2.0 * reference.dry_solids().value(),
        );
        close(
            scaled.elements()[&copper()].value(),
            2.0 * reference.elements()[&copper()].value(),
        );
    }
    assert_eq!(point.fresh_f80, Length::new(200e-6));
    assert_eq!(point.external_p80, Length::new(100e-6));
    assert_eq!(point.prototype, preserved);
    assert_eq!(prototype, preserved);
    assert_eq!(point.prototype.bond(), &bond());
    assert_eq!(point.prototype.circuit(), &circuit(false));
    let inverse = prototype
        .at_fresh_throughput(MassFlowRate::new(24.0))
        .unwrap();
    close(inverse.net_pinion_power.value(), energy * 24.0);
    close(
        inverse.solution.fresh_feed.dry_solids().unwrap().value(),
        24.0,
    );
}

#[test]
fn purge_changes_external_p80_even_when_overflow_alone_is_fine() {
    let split = Partition::from_product(Dimensionless::new(0.5)).unwrap();
    let admitted = circuit(false)
        .solve_with_underflow_recycle(
            &feed(),
            &QueryPoint::new(),
            &QueryPoint::new(),
            settings(),
            split,
        )
        .unwrap();
    let fraction = Dimensionless::new(0.8);
    assert_eq!(
        admitted.product.discrete_passing_size(fraction).unwrap(),
        Length::new(100e-6)
    );
    let external = ParticleStream::mix(&[&admitted.product, &admitted.purge]).unwrap();
    assert_eq!(
        external.discrete_passing_size(fraction).unwrap(),
        Length::new(200e-6)
    );
    // Including coarse purge makes external F80=P80, so the original Bond owner
    // refuses zero reduction instead of charging the overflow-only product.
    assert!(matches!(
        prototype(false, split),
        Err(CircuitPowerError::Bond(ComminutionError::Invalid(_)))
    ));
}

#[test]
fn actual_scaled_flow_must_reenter_original_source_domain() {
    let prototype = prototype(true, return_all()).unwrap();
    let original = prototype.clone();
    prototype
        .at_fresh_throughput(MassFlowRate::new(12.0))
        .unwrap();
    assert!(matches!(
        prototype.at_fresh_throughput(MassFlowRate::new(24.0)),
        Err(CircuitPowerError::Circuit(CircuitError::Breakage(
            BreakageError::OutsideValidity
        )))
    ));
    let power = prototype.specific_energy.value() * 24.0;
    assert!(matches!(
        prototype.at_net_pinion_power(Power::new(power)),
        Err(CircuitPowerError::Circuit(CircuitError::Breakage(
            BreakageError::OutsideValidity
        )))
    ));
    assert_eq!(prototype, original);
}

#[test]
fn independent_tiny_purge_and_original_evidence_survive_scaling() {
    let split =
        Partition::from_fractions(Dimensionless::new(1.0), Dimensionless::new(1e-20)).unwrap();
    let prototype = prototype(false, split).unwrap();
    let point = prototype
        .at_fresh_throughput(MassFlowRate::new(24.0))
        .unwrap();
    assert_eq!(point.solution.underflow_recycle, split);
    assert!(point.solution.purge.dry_solids().unwrap().value() > 0.0);
    assert!(point.solution.purge.water().value() > 0.0);
    assert!(
        point
            .solution
            .purge
            .element_flow(&copper())
            .unwrap()
            .value()
            > 0.0
    );
    assert_eq!(
        point.prototype.conditions(),
        (&QueryPoint::new(), &QueryPoint::new())
    );
    assert_eq!(point.prototype.settings(), settings());
    assert!(point.solution.closure.maximum_absolute_error.value() < 1e-11);
}

#[test]
fn invalid_or_unrepresentable_requests_return_no_operating_point() {
    let prototype = prototype(false, return_all()).unwrap();
    for value in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(
            prototype
                .at_fresh_throughput(MassFlowRate::new(value))
                .is_err()
        );
        assert!(prototype.at_net_pinion_power(Power::new(value)).is_err());
    }
    assert!(
        prototype
            .at_fresh_throughput(MassFlowRate::new(f64::MAX))
            .is_err()
    );
    assert!(
        prototype
            .at_fresh_throughput(MassFlowRate::new(f64::from_bits(1)))
            .is_err()
    );
    assert!(
        prototype
            .at_net_pinion_power(Power::new(f64::from_bits(1)))
            .is_err()
    );
}
