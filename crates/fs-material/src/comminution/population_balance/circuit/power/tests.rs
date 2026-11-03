use super::*;
use crate::comminution::{
    BondSizeDomain, BondTestContext, BondWorkIndex,
    breakage::{BreakageTestContext, BreakageTransfer, SizeClass},
    population_balance::{
        ContinuousPopulationBalanceMill, MEAN_RESIDENCE_TIME_AXIS, MILL_DRY_FLOW_AXIS,
        PopulationBalanceError, SelectionBreakageClass,
    },
};
use crate::separation::{
    CLASSIFIER_LAW_ID, ClassifierCard, PARTICLE_DIAMETER_AXIS, ParticleClass, SOLIDS_FRACTION_AXIS,
};
use fs_blake3::hash_bytes;
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::{Dims, Frequency, MassFlowRate, Pressure, QuantitySpec, chemistry::ElementId};
use std::collections::BTreeMap;

fn source() -> Provenance {
    Provenance {
        source: "authored three-size continuous Bond fixture; not measured ore".into(),
        license: "project fixture".into(),
        artifact: Some(hash_bytes(b"continuous-circuit-bond-power-fixture-v1")),
    }
}
fn fixture(maximum_flow: f64, maximum_p80: f64) -> (ContinuousBondCircuit, ParticleStream) {
    let layout = [("coarse", 400e-6), ("middle", 200e-6), ("fine", 100e-6)];
    let mill = ContinuousPopulationBalanceMill::new(
        layout
            .iter()
            .map(|&(id, diameter)| SizeClass {
                id: id.into(),
                diameter: Length::new(diameter),
            })
            .collect(),
        [
            ("coarse", 0.1, Some("middle")),
            ("middle", 0.05, Some("fine")),
            ("fine", 0.0, None),
        ]
        .into_iter()
        .map(|(input, rate, daughter)| SelectionBreakageClass {
            input: input.into(),
            selection_rate: Frequency::new(rate),
            daughters: daughter
                .into_iter()
                .map(|output| BreakageTransfer {
                    output: output.into(),
                    mass_fraction: Dimensionless::new(1.0),
                })
                .collect(),
        })
        .collect(),
        BreakageTestContext {
            specimen: "synthetic assayed three-size feed".into(),
            method: "declared continuous selection and daughter coefficients".into(),
            provenance: source(),
        },
        ValidityDomain::unconstrained()
            .with_quantity(
                MEAN_RESIDENCE_TIME_AXIS,
                QuantitySpec::dimensional(Time::DIMS),
                0.1,
                100.0,
            )
            .with_quantity(
                MILL_DRY_FLOW_AXIS,
                QuantitySpec::dimensional(MassFlowRate::DIMS),
                1.0,
                maximum_flow,
            )
            .with_quantity(
                SOLIDS_FRACTION_AXIS,
                QuantitySpec::dimensional(Dims::NONE),
                0.0,
                1.0,
            ),
    )
    .unwrap();
    let classifier = ClassifierCard::from_model_card(&ConstitutiveModelCard {
        law: LawId(CLASSIFIER_LAW_ID.into()),
        law_version: 1,
        state_schema_version: 1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        parameters: [
            ("corrected-cut-size", 200e-6, Length::DIMS),
            ("sharpness", 2.0, Dims::NONE),
            ("water-recovery-to-underflow", 0.2, Dims::NONE),
        ]
        .into_iter()
        .map(|(name, value, dims)| (name.into(), LawParameter { value, dims }))
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
            )
            .with_quantity(
                "pressure",
                QuantitySpec::dimensional(Pressure::DIMS),
                100_000.0,
                100_000.0,
            ),
        sources: vec![source().artifact.unwrap()],
        provenance: source(),
    })
    .unwrap();
    let bond = BondGrindingCard::new(
        BondWorkIndex::from_kwh_per_metric_tonne(15.0).unwrap(),
        BondTestContext {
            specimen: "synthetic assayed three-size feed".into(),
            method: "independently authored Bond work-index fixture".into(),
            closing_screen: Length::new(125e-6),
            provenance: source(),
        },
        BondSizeDomain::new(
            [Length::new(200e-6), Length::new(1e-3)],
            [Length::new(70e-6), Length::new(maximum_p80)],
        )
        .unwrap(),
    )
    .unwrap();
    let fresh = ParticleStream::new(
        layout
            .iter()
            .map(|&(id, diameter)| {
                let dry = if id == "coarse" { 10.0 } else { 0.0 };
                ParticleClass::new(
                    id,
                    Length::new(diameter),
                    MassFlowRate::new(dry),
                    BTreeMap::from([(ElementId::new("Cu").unwrap(), MassFlowRate::new(0.2 * dry))]),
                )
                .unwrap()
            })
            .collect(),
        MassFlowRate::new(12.0),
    )
    .unwrap();
    (
        ContinuousBondCircuit::new(ContinuousGrindingCircuit::new(mill, classifier), bond),
        fresh,
    )
}
fn classifier_conditions() -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "pressure",
            QuantitySpec::dimensional(Pressure::DIMS),
            100_000.0,
        )
        .unwrap()
}
fn circuit_settings() -> CircuitSettings {
    CircuitSettings {
        absolute_balance_tolerance: MassFlowRate::new(1e-10),
        relative_balance_tolerance: Dimensionless::new(1e-12),
        maximum_circulating_load: Dimensionless::new(100.0),
    }
}
fn holdup_settings() -> HoldupSolveSettings {
    HoldupSolveSettings {
        minimum_residence: Time::new(1.0),
        maximum_residence: Time::new(100.0),
        absolute_mass_tolerance: Mass::new(1e-9),
        relative_mass_tolerance: Dimensionless::new(1e-12),
        maximum_iterations: 80,
    }
}
fn returned(fraction: f64) -> Partition {
    Partition::from_product(Dimensionless::new(fraction)).unwrap()
}
fn residence(
    model: &ContinuousBondCircuit,
    fresh: &ParticleStream,
    tau: f64,
    recycle: f64,
) -> Result<ContinuousBondCircuitPoint> {
    model.at_residence_time(
        fresh,
        Time::new(tau),
        &QueryPoint::new(),
        &classifier_conditions(),
        circuit_settings(),
        returned(recycle),
    )
}
fn holdup(
    model: &ContinuousBondCircuit,
    fresh: &ParticleStream,
    target: f64,
) -> Result<ContinuousBondCircuitPoint> {
    model.at_dry_holdup(
        fresh,
        Mass::new(target),
        &QueryPoint::new(),
        &classifier_conditions(),
        circuit_settings(),
        returned(1.0),
        holdup_settings(),
    )
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-9 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}
fn energy(p80_micrometres: f64) -> f64 {
    10.0 * 15.0 * (1.0 / p80_micrometres.sqrt() - 1.0 / 400.0_f64.sqrt()) * 3600.0
}
// Independent three-class equations, with full recycle and escape fractions
// (.05, .4, .8*2^(-.25)); selections are (.1, .05, 0) per second.
fn discharge(tau: f64, fresh: f64) -> [f64; 3] {
    let coarse = fresh / (0.05 + 0.1 * tau);
    let middle = 0.1 * tau * coarse / (0.4 + 0.05 * tau);
    let fine = 0.05 * tau * middle / (0.8 * 2.0_f64.powf(-0.25));
    [coarse, middle, fine]
}
fn target_holdup() -> f64 {
    40.0 * discharge(40.0, 10.0).iter().sum::<f64>()
}

#[test]
fn residence_energy_uses_original_state_and_fresh_dry_basis() {
    let (model, fresh) = fixture(1000.0, 500e-6);
    let saved = (model.clone(), fresh.clone());
    let point = residence(&model, &fresh, 10.0, 1.0).unwrap();
    let original = model
        .circuit()
        .at_residence_time(
            &fresh,
            Time::new(10.0),
            &QueryPoint::new(),
            &classifier_conditions(),
            circuit_settings(),
            returned(1.0),
        )
        .unwrap();
    assert_eq!(point.solution(), &original);
    assert_eq!(point.model(), &model);
    assert_eq!(point.model().bond(), model.bond());
    assert!(point.holdup_solution().is_none());
    assert_eq!(point.fresh_f80, Length::new(400e-6));
    assert_eq!(point.external_p80, Length::new(200e-6));
    close(point.specific_energy.value(), energy(200.0));
    close(point.net_pinion_power.value(), 10.0 * energy(200.0));
    let circulating = discharge(10.0, 10.0).iter().sum::<f64>();
    close(
        point
            .solution()
            .circuit
            .mill_feed
            .dry_solids()
            .unwrap()
            .value(),
        circulating,
    );
    close(
        point.solution().mill.conditions.axes()[MILL_DRY_FLOW_AXIS],
        circulating,
    );
    assert!(circulating > 2.0 * fresh.dry_solids().unwrap().value());
    assert!(point.net_pinion_power.value() < circulating * energy(200.0));
    assert_eq!((model, fresh), saved);
}

#[test]
fn fixed_holdup_changes_psd_p80_and_energy_with_throughput() {
    let (model, fresh) = fixture(1000.0, 500e-6);
    let doubled = fresh.scaled(Dimensionless::new(2.0)).unwrap();
    let target = target_holdup();
    let low = holdup(&model, &fresh, target).unwrap();
    let high = holdup(&model, &doubled, target).unwrap();
    for (point, feed, tau, p80) in [
        (&low, &fresh, 40.0, 100.0),
        (&high, &doubled, 15.775944375965341, 200.0),
    ] {
        let original = model
            .circuit()
            .at_dry_holdup(
                feed,
                Mass::new(target),
                &QueryPoint::new(),
                &classifier_conditions(),
                circuit_settings(),
                returned(1.0),
                holdup_settings(),
            )
            .unwrap();
        assert_eq!(point.holdup_solution(), Some(&original));
        assert_eq!(point.solution(), &original.solution);
        close(point.solution().mill.mean_residence_time.value(), tau);
        close(point.solution().mill.dry_solids_holdup.value(), target);
        close(
            point
                .solution()
                .circuit
                .mill_discharge
                .dry_solids()
                .unwrap()
                .value(),
            discharge(tau, feed.dry_solids().unwrap().value())
                .iter()
                .sum(),
        );
        close(point.external_p80.value(), p80 * 1e-6);
        close(point.specific_energy.value(), energy(p80));
        close(
            point.net_pinion_power.value(),
            feed.dry_solids().unwrap().value() * energy(p80),
        );
        assert!(original.mass_residual.value().abs() <= original.allowed_mass_error.value());
    }
    assert_eq!(low.external_p80, Length::new(100e-6));
    assert_eq!(high.external_p80, Length::new(200e-6));
    // Quantile coarsening can outweigh a throughput increase. This is not a
    // smooth monotone power/capacity inverse or a scaled invariant-PSD prototype.
    assert!(high.net_pinion_power < low.net_pinion_power);
    assert_ne!(high.specific_energy, low.specific_energy);
}

#[test]
fn purge_participates_in_external_p80_and_zero_reduction_refuses() {
    let (model, fresh) = fixture(1000.0, 500e-6);
    let point = residence(&model, &fresh, 40.0, 0.75).unwrap();
    let flows = &point.solution().circuit;
    assert_eq!(
        flows
            .product
            .discrete_passing_size(Dimensionless::new(0.8))
            .unwrap(),
        Length::new(100e-6)
    );
    assert!(flows.purge.dry_solids().unwrap().value() > 0.0);
    assert_eq!(point.external_p80, Length::new(200e-6));
    close(point.net_pinion_power.value(), 10.0 * energy(200.0));
    close(
        flows.product.dry_solids().unwrap().value() + flows.purge.dry_solids().unwrap().value(),
        10.0,
    );

    let admitted = model
        .circuit()
        .at_residence_time(
            &fresh,
            Time::new(10.0),
            &QueryPoint::new(),
            &classifier_conditions(),
            circuit_settings(),
            returned(0.75),
        )
        .unwrap();
    let external =
        ParticleStream::mix(&[&admitted.circuit.product, &admitted.circuit.purge]).unwrap();
    assert_eq!(
        external
            .discrete_passing_size(Dimensionless::new(0.8))
            .unwrap(),
        Length::new(400e-6)
    );
    let original = model
        .bond()
        .specific_energy(Length::new(400e-6), Length::new(400e-6))
        .unwrap_err();
    assert_eq!(
        residence(&model, &fresh, 10.0, 0.75).unwrap_err(),
        ContinuousBondPowerError::Bond(original)
    );
}

#[test]
fn new_actual_p80_reenters_original_bond_support() {
    let (model, fresh) = fixture(1000.0, 150e-6);
    let doubled = fresh.scaled(Dimensionless::new(2.0)).unwrap();
    assert_eq!(
        holdup(&model, &fresh, target_holdup())
            .unwrap()
            .external_p80,
        Length::new(100e-6)
    );
    let expected = model
        .bond()
        .specific_energy(Length::new(400e-6), Length::new(200e-6))
        .unwrap_err();
    assert!(matches!(
        expected,
        ComminutionError::OutsideValidity { axis: "P80", .. }
    ));
    assert_eq!(
        holdup(&model, &doubled, target_holdup()).unwrap_err(),
        ContinuousBondPowerError::Bond(expected)
    );
}

#[test]
fn each_mode_preserves_original_circulating_flow_source_refusals() {
    // Fresh flow is supported. At tau=1, circulating flow is outside support.
    let (model, fresh) = fixture(50.0, 500e-6);
    let saved = (model.clone(), fresh.clone());
    let original = ContinuousCircuitError::Population(PopulationBalanceError::OutsideValidity);
    assert_eq!(
        residence(&model, &fresh, 1.0, 1.0).unwrap_err(),
        ContinuousBondPowerError::Circuit(original.clone())
    );
    assert_eq!(
        holdup(&model, &fresh, target_holdup()).unwrap_err(),
        ContinuousBondPowerError::Holdup(HoldupError::Circuit {
            residence: Time::new(1.0),
            error: original
        })
    );
    assert_eq!((model, fresh), saved);
}

#[test]
fn holdup_bracket_and_classifier_context_refusals_are_not_replaced() {
    let (model, fresh) = fixture(1000.0, 500e-6);
    let original = model
        .circuit()
        .at_dry_holdup(
            &fresh,
            Mass::new(1.0),
            &QueryPoint::new(),
            &classifier_conditions(),
            circuit_settings(),
            returned(1.0),
            holdup_settings(),
        )
        .unwrap_err();
    assert!(matches!(original, HoldupError::Unbracketed { .. }));
    assert_eq!(
        holdup(&model, &fresh, 1.0).unwrap_err(),
        ContinuousBondPowerError::Holdup(original)
    );
    let original = model
        .circuit()
        .at_residence_time(
            &fresh,
            Time::new(10.0),
            &QueryPoint::new(),
            &QueryPoint::new(),
            circuit_settings(),
            returned(1.0),
        )
        .unwrap_err();
    let refused = model
        .at_residence_time(
            &fresh,
            Time::new(10.0),
            &QueryPoint::new(),
            &QueryPoint::new(),
            circuit_settings(),
            returned(1.0),
        )
        .unwrap_err();
    assert_eq!(refused, ContinuousBondPowerError::Circuit(original));
}
