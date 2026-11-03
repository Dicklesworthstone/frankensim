use super::*;
use crate::comminution::{
    breakage::{BreakageTestContext, BreakageTransfer, SizeClass},
    circuit::CircuitError,
    population_balance::{
        ContinuousPopulationBalanceMill, MEAN_RESIDENCE_TIME_AXIS, MILL_DRY_FLOW_AXIS,
        PopulationBalanceError, SelectionBreakageClass,
    },
};
use crate::separation::{
    CLASSIFIER_LAW_ID, ClassifierCard, PARTICLE_DIAMETER_AXIS, ParticleClass, SOLIDS_FRACTION_AXIS,
    SeparationError,
};
use fs_blake3::hash_bytes;
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::{Dims, Frequency, Length, MassFlowRate, Pressure, QuantitySpec, chemistry::ElementId};
use std::collections::BTreeMap;

fn copper() -> ElementId {
    ElementId::new("Cu").unwrap()
}
fn fixture(inert: bool, maximum_flow: f64) -> (ContinuousGrindingCircuit, ParticleStream) {
    let source = Provenance {
        source: "authored hold-up arithmetic fixture, not ore testwork".into(),
        license: "project fixture".into(),
        artifact: Some(hash_bytes(b"continuous-circuit-holdup-fixture-v1")),
    };
    let layout = if inert {
        vec![("coarse", 200e-6)]
    } else {
        vec![("coarse", 200e-6), ("fine", 100e-6)]
    };
    let mill = ContinuousPopulationBalanceMill::new(
        layout
            .iter()
            .map(|&(id, diameter)| SizeClass {
                id: id.into(),
                diameter: Length::new(diameter),
            })
            .collect(),
        layout
            .iter()
            .map(|&(id, _)| {
                let breaks = !inert && id == "coarse";
                SelectionBreakageClass {
                    input: id.into(),
                    selection_rate: Frequency::new(if breaks { 0.1 } else { 0.0 }),
                    daughters: if breaks {
                        vec![BreakageTransfer {
                            output: "fine".into(),
                            mass_fraction: Dimensionless::new(1.0),
                        }]
                    } else {
                        vec![]
                    },
                }
            })
            .collect(),
        BreakageTestContext {
            specimen: "synthetic assayed feed".into(),
            method: "declared continuous selection/daughter coefficients".into(),
            provenance: source.clone(),
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
            ("corrected-cut-size", 100e-6, Length::DIMS),
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
                1e-6,
                1.0,
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
        sources: vec![source.artifact.unwrap()],
        provenance: source,
    })
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
                    BTreeMap::from([(copper(), MassFlowRate::new(0.2 * dry))]),
                )
                .unwrap()
            })
            .collect(),
        MassFlowRate::new(12.0),
    )
    .unwrap();
    (ContinuousGrindingCircuit::new(mill, classifier), fresh)
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
fn settings() -> HoldupSolveSettings {
    HoldupSolveSettings {
        minimum_residence: Time::new(1.0),
        maximum_residence: Time::new(20.0),
        absolute_mass_tolerance: Mass::new(1e-9),
        relative_mass_tolerance: Dimensionless::new(1e-12),
        maximum_iterations: 80,
    }
}
fn recycle(returned: f64) -> Partition {
    Partition::from_product(Dimensionless::new(returned)).unwrap()
}
fn solve(
    model: &ContinuousGrindingCircuit,
    fresh: &ParticleStream,
    target: f64,
    returned: f64,
    numerics: HoldupSolveSettings,
) -> Result<HoldupSolution, HoldupError> {
    model.at_dry_holdup(
        fresh,
        Mass::new(target),
        &QueryPoint::new(),
        &classifier_conditions(),
        circuit_settings(),
        recycle(returned),
        numerics,
    )
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-9 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn coupled_holdup_and_doubled_feed_match_independent_two_bin_solution() {
    let (model, fresh) = fixture(false, 1000.0);
    let saved = fresh.clone();
    let target = 1000.0 / 3.0;
    let doubled = ParticleStream::mix(&[&fresh, &fresh]).unwrap();
    // H/f = tau*(.4 + .1*tau)/[.4*(.05 + .1*tau)]. At f=20 and
    // H=1000/3, 3*tau^2 - 8*tau - 10 = 0, independent of bisection.
    for (feed, expected_time) in [(&fresh, 10.0), (&doubled, (4.0 + 46.0_f64.sqrt()) / 3.0)] {
        let result = solve(&model, feed, target, 1.0, settings()).unwrap();
        let original = &result.solution;
        close(original.mill.mean_residence_time.value(), expected_time);
        close(original.mill.dry_solids_holdup.value(), target);
        close(
            original.circuit.mill_feed.dry_solids().unwrap().value(),
            target / expected_time,
        );
        close(
            original.circuit.product.dry_solids().unwrap().value(),
            feed.dry_solids().unwrap().value(),
        );
        let fine = original
            .circuit
            .product
            .classes()
            .iter()
            .find(|class| class.id() == "fine")
            .unwrap();
        close(
            fine.dry_solids().value(),
            feed.dry_solids().unwrap().value() * 0.1 * expected_time / (0.05 + 0.1 * expected_time),
        );
        close(
            original
                .circuit
                .product
                .element_flow(&copper())
                .unwrap()
                .value(),
            feed.element_flow(&copper()).unwrap().value(),
        );
        close(
            original.circuit.product.water().value(),
            feed.water().value(),
        );
        close(
            original.mill.water_holdup.value(),
            expected_time * feed.water().value() / 0.8,
        );
        assert_eq!(original.model(), &model);
        assert_eq!(original.mill.model(), model.mill());
        assert_eq!(result.target_dry_holdup, Mass::new(target));
        assert_eq!(result.settings, settings());
        assert_eq!(
            result.mass_residual.value(),
            original.mill.dry_solids_holdup.value() - target
        );
        assert!(result.mass_residual.value().abs() <= result.allowed_mass_error.value());
        assert!(result.iterations > 0 && result.iterations <= settings().maximum_iterations);
        // The returned original source query is the solved circulating flow, not fresh flow.
        close(
            original.mill.conditions.axes()[MILL_DRY_FLOW_AXIS],
            original.circuit.mill_feed.dry_solids().unwrap().value(),
        );
        close(
            original.mill.conditions.axes()[MEAN_RESIDENCE_TIME_AXIS],
            expected_time,
        );
        assert_eq!(
            *original,
            model
                .at_residence_time(
                    feed,
                    original.mill.mean_residence_time,
                    &QueryPoint::new(),
                    &classifier_conditions(),
                    circuit_settings(),
                    recycle(1.0)
                )
                .unwrap()
        );
    }
    assert_eq!(fresh, saved);
}

#[test]
fn purge_and_no_recycle_close_dry_inventory_separately_from_water() {
    let (model, fresh) = fixture(false, 1000.0);
    // Return .75 gives coarse/fine escape .2875/.55. For H=1000/3,
    // f=10, the positive root is (77 + sqrt(13519))/12.
    let target = 1000.0 / 3.0;
    let purged = solve(&model, &fresh, target, 0.75, settings()).unwrap();
    close(
        purged.solution.mill.mean_residence_time.value(),
        (77.0 + 13519.0_f64.sqrt()) / 12.0,
    );
    close(purged.solution.mill.dry_solids_holdup.value(), target);
    let flows = &purged.solution.circuit;
    assert!(flows.purge.dry_solids().unwrap().value() > 0.0);
    close(
        flows.product.dry_solids().unwrap().value() + flows.purge.dry_solids().unwrap().value(),
        10.0,
    );
    close(
        flows.product.element_flow(&copper()).unwrap().value()
            + flows.purge.element_flow(&copper()).unwrap().value(),
        2.0,
    );
    close(
        flows.product.water().value() + flows.purge.water().value(),
        12.0,
    );

    let once = solve(&model, &fresh, 100.0, 0.0, settings()).unwrap();
    // Without recycle H/fresh determines residence, regardless of the classifier split.
    close(once.solution.mill.mean_residence_time.value(), 100.0 / 10.0);
    close(once.solution.mill.dry_solids_holdup.value(), 100.0);
    close(once.solution.mill.water_holdup.value(), 120.0);
    assert_eq!(once.solution.circuit.mill_feed, fresh);
    assert_eq!(
        once.solution.circuit.recycle.dry_solids().unwrap().value(),
        0.0
    );
}

#[test]
fn positive_unbracketed_targets_are_not_extrapolated_or_clipped() {
    let (model, fresh) = fixture(false, 1000.0);
    for target in [10.0, 1000.0] {
        let HoldupError::Unbracketed {
            target: reported,
            bracket,
        } = solve(&model, &fresh, target, 1.0, settings()).unwrap_err()
        else {
            panic!("target must refuse")
        };
        assert_eq!(reported.value(), target);
        assert_eq!(bracket.lower_residence, Time::new(1.0));
        assert_eq!(bracket.upper_residence, Time::new(20.0));
        close(bracket.lower_holdup.value(), 10.0 * 0.5 / (0.4 * 0.15));
        close(bracket.upper_holdup.value(), 200.0 * 2.4 / (0.4 * 2.05));
    }
}

#[test]
fn source_refusals_keep_the_original_error_and_leave_inputs_untouched() {
    // Fresh 10 kg/s is supported, and the desired tau=10 solution is supported;
    // the declared lower bracket at tau=1 circulates >50 kg/s and must refuse.
    let (limited, fresh) = fixture(false, 50.0);
    let saved = fresh.clone();
    assert!(
        limited
            .at_residence_time(
                &fresh,
                Time::new(10.0),
                &QueryPoint::new(),
                &classifier_conditions(),
                circuit_settings(),
                recycle(1.0)
            )
            .is_ok()
    );
    assert_eq!(
        solve(&limited, &fresh, 1000.0 / 3.0, 1.0, settings()).unwrap_err(),
        HoldupError::Circuit {
            residence: Time::new(1.0),
            error: ContinuousCircuitError::Population(PopulationBalanceError::OutsideValidity)
        }
    );
    assert_eq!(fresh, saved);

    let (model, fresh) = fixture(true, 1000.0);
    // Even when the lower endpoint exactly meets the target, the declared upper
    // endpoint must be admitted before any result can succeed.
    let outside = HoldupSolveSettings {
        maximum_residence: Time::new(101.0),
        ..settings()
    };
    assert_eq!(
        solve(&model, &fresh, 10.0, 0.0, outside).unwrap_err(),
        HoldupError::Circuit {
            residence: Time::new(101.0),
            error: ContinuousCircuitError::Population(PopulationBalanceError::OutsideValidity)
        }
    );
    assert_eq!(
        model
            .at_dry_holdup(
                &fresh,
                Mass::new(100.0),
                &QueryPoint::new(),
                &QueryPoint::new(),
                circuit_settings(),
                recycle(0.0),
                settings()
            )
            .unwrap_err(),
        HoldupError::Circuit {
            residence: Time::new(1.0),
            error: ContinuousCircuitError::Circuit(CircuitError::Separation(
                SeparationError::OutsideValidity
            ))
        }
    );
}

#[test]
fn finite_iteration_budget_returns_no_partial_accepted_result() {
    let (model, fresh) = fixture(false, 1000.0);
    let saved = fresh.clone();
    let bounded = HoldupSolveSettings {
        maximum_iterations: 1,
        ..settings()
    };
    let HoldupError::IterationExhausted {
        target,
        bracket,
        iterations,
    } = solve(&model, &fresh, 1000.0 / 3.0, 1.0, bounded).unwrap_err()
    else {
        panic!("budget must refuse")
    };
    assert_eq!(iterations, 1);
    assert_eq!(target, Mass::new(1000.0 / 3.0));
    assert_eq!(bracket.upper_residence, Time::new(10.5));
    assert!(bracket.lower_holdup < target && target < bracket.upper_holdup);
    assert_eq!(fresh, saved);
}

#[test]
fn invalid_inputs_and_overflowed_combined_mass_tolerance_refuse() {
    let (model, fresh) = fixture(false, 1000.0);
    for target in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(matches!(
            solve(&model, &fresh, target, 1.0, settings()),
            Err(HoldupError::Invalid(_))
        ));
    }
    for invalid in [
        HoldupSolveSettings {
            minimum_residence: Time::new(0.0),
            ..settings()
        },
        HoldupSolveSettings {
            minimum_residence: Time::new(20.0),
            ..settings()
        },
        HoldupSolveSettings {
            maximum_residence: Time::new(f64::INFINITY),
            ..settings()
        },
        HoldupSolveSettings {
            absolute_mass_tolerance: Mass::new(0.0),
            ..settings()
        },
        HoldupSolveSettings {
            relative_mass_tolerance: Dimensionless::new(0.0),
            ..settings()
        },
        HoldupSolveSettings {
            relative_mass_tolerance: Dimensionless::new(1.0),
            ..settings()
        },
        HoldupSolveSettings {
            maximum_iterations: 0,
            ..settings()
        },
    ] {
        assert!(matches!(
            solve(&model, &fresh, 1000.0 / 3.0, 1.0, invalid),
            Err(HoldupError::Invalid(_))
        ));
    }
    let overflow = HoldupSolveSettings {
        absolute_mass_tolerance: Mass::new(f64::MAX),
        relative_mass_tolerance: Dimensionless::new(0.5),
        ..settings()
    };
    assert_eq!(
        solve(&model, &fresh, f64::MAX, 1.0, overflow).unwrap_err(),
        HoldupError::Invalid("combined dry-mass tolerance overflows")
    );
}

#[test]
fn adjacent_float_residences_refuse_instead_of_claiming_mass_accuracy() {
    let (model, fresh) = fixture(true, 1000.0);
    // Adjacent residence values; 10*tau spans two representable mass steps.
    // The target mass between them has no representable residence in this bracket.
    let lower = f64::from_bits(1.0_f64.to_bits() + 2);
    let upper = f64::from_bits(1.0_f64.to_bits() + 3);
    let target = f64::from_bits(10.0_f64.to_bits() + 3);
    let numerics = HoldupSolveSettings {
        minimum_residence: Time::new(lower),
        maximum_residence: Time::new(upper),
        absolute_mass_tolerance: Mass::new(1e-30),
        relative_mass_tolerance: Dimensionless::new(1e-30),
        ..settings()
    };
    let HoldupError::NumericalResolution {
        reason,
        bracket,
        iterations,
        ..
    } = solve(&model, &fresh, target, 0.0, numerics).unwrap_err()
    else {
        panic!("resolution must refuse")
    };
    assert_eq!(reason, "no representable interior residence remains");
    assert_eq!(iterations, 0);
    assert_eq!(bracket.lower_holdup.value(), 10.0 * lower);
    assert_eq!(bracket.upper_holdup.value(), 10.0 * upper);
    assert!(bracket.lower_holdup.value() < target && target < bracket.upper_holdup.value());
}
