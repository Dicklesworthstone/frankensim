//! Independent analytical fixtures, not fitted ore or plant performance.

#[path = "sump/tests.rs"]
mod sump;

use super::*;
use crate::comminution::{
    breakage::{BreakageTestContext, BreakageTransfer, SizeClass},
    circuit::CircuitSettings,
    population_balance::{
        ClassHoldup, MEAN_RESIDENCE_TIME_AXIS, MILL_DRY_FLOW_AXIS, SelectionBreakageClass,
    },
};
use crate::separation::{
    CLASSIFIER_LAW_ID, ClassifierCard, PARTICLE_DIAMETER_AXIS, ParticleClass, SOLIDS_FRACTION_AXIS,
};
use fs_blake3::hash_bytes;
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::{Dims, Frequency, Length, Mass, MassFlowRate, Pressure, QuantitySpec};
use std::collections::BTreeMap;

fn element(s: &str) -> ElementId {
    ElementId::new(s).unwrap()
}
fn provenance() -> Provenance {
    Provenance {
        source: "constructed dynamic-circuit analytical fixture; not testwork".into(),
        license: "project fixture".into(),
        artifact: Some(hash_bytes(b"dynamic-closed-grinding-v1")),
    }
}
fn model(max_flow: f64, solids: (f64, f64), water_return: f64) -> ContinuousGrindingCircuit {
    let mill = ContinuousPopulationBalanceMill::new(
        vec![
            SizeClass {
                id: "coarse".into(),
                diameter: Length::new(200e-6),
            },
            SizeClass {
                id: "fine".into(),
                diameter: Length::new(100e-6),
            },
        ],
        vec![
            SelectionBreakageClass {
                input: "coarse".into(),
                selection_rate: Frequency::new(0.1),
                daughters: vec![BreakageTransfer {
                    output: "fine".into(),
                    mass_fraction: Dimensionless::new(1.0),
                }],
            },
            SelectionBreakageClass {
                input: "fine".into(),
                selection_rate: Frequency::new(0.0),
                daughters: vec![],
            },
        ],
        BreakageTestContext {
            specimen: "synthetic two-size feed".into(),
            method: "authored rate and daughters".into(),
            provenance: provenance(),
        },
        ValidityDomain::unconstrained()
            .with_quantity(
                MEAN_RESIDENCE_TIME_AXIS,
                QuantitySpec::dimensional(Time::DIMS),
                1.0,
                100.0,
            )
            .with_quantity(
                MILL_DRY_FLOW_AXIS,
                QuantitySpec::dimensional(MassFlowRate::DIMS),
                1.0,
                max_flow,
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
            ("sharpness", 1.0, Dims::NONE),
            ("water-recovery-to-underflow", water_return, Dims::NONE),
        ]
        .into_iter()
        .map(|(name, value, dims)| (name.into(), LawParameter { value, dims }))
        .collect(),
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                PARTICLE_DIAMETER_AXIS,
                QuantitySpec::dimensional(Length::DIMS),
                100e-6,
                200e-6,
            )
            .with_quantity(
                SOLIDS_FRACTION_AXIS,
                QuantitySpec::dimensional(Dims::NONE),
                solids.0,
                solids.1,
            )
            .with_quantity(
                "pressure",
                QuantitySpec::dimensional(Pressure::DIMS),
                100_000.0,
                100_000.0,
            ),
        sources: vec![provenance().artifact.unwrap()],
        provenance: provenance(),
    })
    .unwrap();
    ContinuousGrindingCircuit::new(mill, classifier)
}
fn standard_model() -> ContinuousGrindingCircuit {
    model(1000.0, (0.0, 1.0), 0.2)
}
fn conditions() -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "pressure",
            QuantitySpec::dimensional(Pressure::DIMS),
            100_000.0,
        )
        .unwrap()
}
fn fresh(copper: f64) -> ParticleStream {
    ParticleStream::new(
        [("coarse", 200e-6, 10.0, copper), ("fine", 100e-6, 0.0, 0.0)]
            .into_iter()
            .map(|(id, diameter, dry, cu)| {
                ParticleClass::new(
                    id,
                    Length::new(diameter),
                    MassFlowRate::new(dry),
                    BTreeMap::from([
                        (element("Cu"), MassFlowRate::new(cu)),
                        (element("Fe"), MassFlowRate::new(2.0 * cu)),
                    ]),
                )
                .unwrap()
            })
            .collect(),
        MassFlowRate::new(5.0),
    )
    .unwrap()
}
fn initial() -> MillInventory {
    MillInventory::new(
        [("coarse", 200e-6, 50.0, 5.0), ("fine", 100e-6, 40.0, 0.8)]
            .into_iter()
            .map(|(id, diameter, dry, cu)| ClassHoldup {
                id: id.into(),
                diameter: Length::new(diameter),
                dry_solids: Mass::new(dry),
                elements: BTreeMap::from([
                    (element("Cu"), Mass::new(cu)),
                    (element("Fe"), Mass::new(2.0 * cu)),
                ]),
            })
            .collect(),
        Mass::new(30.0),
    )
    .unwrap()
}
fn settings(step: f64) -> TransientCircuitSettings {
    TransientCircuitSettings {
        time: TransientMillSettings {
            maximum_step: Time::new(step),
            maximum_steps: 100_000,
            absolute_balance_tolerance: Mass::new(1e-9),
            relative_balance_tolerance: Dimensionless::new(1e-11),
        },
        maximum_circulating_load: Dimensionless::new(10.0),
    }
}
fn recycle(fraction: f64) -> Partition {
    Partition::from_product(Dimensionless::new(fraction)).unwrap()
}
fn advance(
    model: &ContinuousGrindingCircuit,
    initial: &MillInventory,
    fresh: &ParticleStream,
    t: f64,
    dt: f64,
    returned: f64,
) -> TransientCircuitStep {
    model
        .advance_inventory(
            initial,
            fresh,
            Time::new(10.0),
            Time::new(t),
            &QueryPoint::new(),
            &conditions(),
            recycle(returned),
            settings(dt),
            || false,
        )
        .unwrap()
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-10 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}
fn mass(inventory: &MillInventory, class: &str) -> f64 {
    inventory
        .classes()
        .iter()
        .find(|c| c.id == class)
        .unwrap()
        .dry_solids
        .value()
}
fn metal(inventory: &MillInventory, class: &str) -> f64 {
    inventory
        .classes()
        .iter()
        .find(|c| c.id == class)
        .unwrap()
        .elements[&element("Cu")]
        .value()
}
fn closes(step: &TransientCircuitStep) {
    for balance in [&step.mill_balance, &step.boundary_balance] {
        close(balance.dry_solids.value(), 0.0);
        close(balance.water.value(), 0.0);
        for residual in balance.elements.values() {
            close(residual.value(), 0.0);
        }
    }
    for channel in [None, Some(element("Cu")), Some(element("Fe"))] {
        let total = |i: &MillInventory| {
            channel.as_ref().map_or_else(
                || i.dry_solids().unwrap().value(),
                |e| i.element_mass(e).unwrap().value(),
            )
        };
        close(
            total(&step.integrated_mill_discharge),
            total(&step.integrated_product) + total(&step.integrated_underflow),
        );
        close(
            total(&step.integrated_underflow),
            total(&step.integrated_recycle) + total(&step.integrated_purge),
        );
        close(
            total(step.initial()) + total(&step.integrated_fresh),
            total(&step.inventory)
                + total(&step.integrated_product)
                + total(&step.integrated_purge),
        );
    }
    close(
        step.initial().water().value() + step.integrated_fresh.water().value(),
        step.inventory.water().value()
            + step.integrated_product.water().value()
            + step.integrated_purge.water().value(),
    );
}

#[test]
fn independent_two_size_implicit_answer_keeps_every_material_and_internal_account() {
    let model = standard_model();
    let step = advance(&model, &initial(), &fresh(1.0), 1.0, 1.0, 0.75);
    // At d50 and 2*d50 with sharpness1 and water split0.2, underflow is
    // 0.6 and0.8. With return0.75: escape_f=.55, escape_c=.4, escape_w=.85.
    // The independent triangular one-second equations give these exact rationals.
    close(mass(&step.inventory, "coarse"), 1000.0 / 19.0);
    close(mass(&step.inventory, "fine"), 172000.0 / 4009.0);
    close(step.inventory.water().value(), 1000.0 / 31.0);
    close(metal(&step.inventory, "coarse"), 100.0 / 19.0);
    close(metal(&step.inventory, "fine"), 5040.0 / 4009.0);
    close(
        step.initial_point.mill_feed.dry_solids().unwrap().value(),
        14.8,
    );
    assert!(step.final_point.mill_conditions.axes()[MILL_DRY_FLOW_AXIS] > 15.0);
    assert_eq!(step.model(), &model);
    assert_eq!(step.initial(), &initial());
    assert_eq!(step.fresh_feed(), &fresh(1.0));
    assert_eq!(step.substeps, 1);
    closes(&step);
}

#[test]
fn time_refinement_converges_to_independent_continuous_ode_solution() {
    let model = standard_model();
    let t = 20.0;
    let (a, b) = (0.14_f64, 0.055_f64);
    let coarse_steady = 10.0 / a;
    let coarse_exact = coarse_steady + (50.0 - coarse_steady) * (-a * t).exp();
    let fine_exact = 40.0 * (-b * t).exp()
        + (0.1 * coarse_steady / b) * (1.0 - (-b * t).exp())
        + 0.1 * (50.0 - coarse_steady) * ((-a * t).exp() - (-b * t).exp()) / (b - a);
    let water_exact = 5.0 / 0.085 + (30.0 - 5.0 / 0.085) * (-0.085_f64 * t).exp();
    let coarse = advance(&model, &initial(), &fresh(1.0), t, 1.0, 0.75);
    let fine = advance(&model, &initial(), &fresh(1.0), t, 0.5, 0.75);
    for (a, b, expected) in [
        (
            mass(&coarse.inventory, "coarse"),
            mass(&fine.inventory, "coarse"),
            coarse_exact,
        ),
        (
            mass(&coarse.inventory, "fine"),
            mass(&fine.inventory, "fine"),
            fine_exact,
        ),
        (
            coarse.inventory.water().value(),
            fine.inventory.water().value(),
            water_exact,
        ),
    ] {
        let ratio = (a - expected).abs() / (b - expected).abs();
        assert!(
            (1.8..2.2).contains(&ratio),
            "first-order error ratio {ratio}"
        );
    }
    closes(&coarse);
    closes(&fine);
}

#[test]
fn zero_return_reuses_original_open_mill_with_identical_final_inventory() {
    let model = standard_model();
    let start = initial();
    let feed = fresh(1.0);
    let closed = advance(&model, &start, &feed, 20.0, 0.5, 0.0);
    let open = model
        .mill()
        .advance_inventory(
            &start,
            &feed,
            Time::new(10.0),
            Time::new(20.0),
            &QueryPoint::new(),
            settings(0.5).time,
            || false,
        )
        .unwrap();
    assert_eq!(closed.inventory, open.inventory);
    assert_eq!(closed.final_point.mill_discharge, open.outlet);
    assert_eq!(closed.integrated_recycle.dry_solids().unwrap().value(), 0.0);
    close(
        closed
            .integrated_mill_discharge
            .dry_solids()
            .unwrap()
            .value(),
        open.integrated_discharge.dry_solids().unwrap().value(),
    );
    closes(&closed);
}

fn steady_inventory(model: &ContinuousGrindingCircuit) -> MillInventory {
    let steady = model
        .at_residence_time(
            &fresh(1.0),
            Time::new(10.0),
            &QueryPoint::new(),
            &conditions(),
            CircuitSettings {
                absolute_balance_tolerance: MassFlowRate::new(1e-10),
                relative_balance_tolerance: Dimensionless::new(1e-12),
                maximum_circulating_load: Dimensionless::new(10.0),
            },
            recycle(0.75),
        )
        .unwrap();
    MillInventory::new(steady.mill.class_holdup, steady.mill.water_holdup).unwrap()
}

#[test]
fn original_steady_circuit_is_a_fixed_point_and_a_grade_step_has_retained_inventory_delay() {
    let model = standard_model();
    let start = steady_inventory(&model);
    close(mass(&start, "coarse"), 500.0 / 7.0);
    close(mass(&start, "fine"), 10000.0 / 77.0);
    let steady = advance(&model, &start, &fresh(1.0), 50.0, 0.5, 0.75);
    close(mass(&steady.inventory, "coarse"), mass(&start, "coarse"));
    close(mass(&steady.inventory, "fine"), mass(&start, "fine"));
    close(steady.inventory.water().value(), start.water().value());
    closes(&steady);
    let pulse = advance(&model, &start, &fresh(2.0), 1.0, 1.0, 0.75);
    // Only the source Cu/Fe grade changes. Exact independent one-step Cu masses:
    close(metal(&pulse.inventory, "coarse"), 3200.0 / 399.0);
    close(metal(&pulse.inventory, "fine"), 12104000.0 / 926079.0);
    close(
        pulse.inventory.dry_solids().unwrap().value(),
        start.dry_solids().unwrap().value(),
    );
    let external_cu = pulse
        .integrated_product
        .element_mass(&element("Cu"))
        .unwrap()
        .value()
        + pulse
            .integrated_purge
            .element_mass(&element("Cu"))
            .unwrap()
            .value();
    assert!(external_cu > 1.0 && external_cu < 2.0);
    assert!(
        pulse.inventory.element_mass(&element("Cu")).unwrap()
            > start.element_mass(&element("Cu")).unwrap()
    );
    closes(&pulse);
    let later = advance(&model, &pulse.inventory, &fresh(2.0), 50.0, 0.5, 0.75);
    assert!(
        later
            .final_point
            .product
            .element_grade(&element("Cu"))
            .unwrap()
            .unwrap()
            > pulse
                .final_point
                .product
                .element_grade(&element("Cu"))
                .unwrap()
                .unwrap()
    );
    closes(&later);
}

#[test]
fn every_actual_circulating_and_classifier_point_must_remain_in_support() {
    let start = initial();
    let feed = fresh(1.0);
    let evaluate = |model: ContinuousGrindingCircuit, cfg| {
        model.advance_inventory(
            &start,
            &feed,
            Time::new(10.0),
            Time::new(1.0),
            &QueryPoint::new(),
            &conditions(),
            recycle(0.75),
            cfg,
            || false,
        )
    };
    // Initial circulating feed14.8 is in support; the first endpoint exceeds14.9.
    assert!(matches!(
        evaluate(model(14.9, (0.0, 1.0), 0.2), settings(1.0)),
        Err(TransientCircuitError::Mill(TransientMillError::Population(
            PopulationBalanceError::OutsideValidity
        )))
    ));
    // Initial classifier solids fraction.75 is supported; endpoint11873/15882 is not.
    assert!(matches!(
        evaluate(model(1000.0, (0.749, 0.751), 0.2), settings(1.0)),
        Err(TransientCircuitError::Circuit(CircuitError::Separation(
            SeparationError::OutsideValidity
        )))
    ));
    let mut limit = settings(1.0);
    limit.maximum_circulating_load = Dimensionless::new(0.49);
    assert!(matches!(
        evaluate(standard_model(), limit),
        Err(TransientCircuitError::Circuit(
            CircuitError::CirculatingLoadExceeded { .. }
        ))
    ));
    assert_eq!(start, initial());
    assert_eq!(feed, fresh(1.0));
}

#[test]
fn a_trapped_loop_can_accumulate_finite_inventory_without_claiming_steady_flow() {
    let model = model(1000.0, (0.0, 1.0), 1.0);
    let result = advance(&model, &initial(), &fresh(1.0), 1.0, 1.0, 1.0);
    close(result.inventory.dry_solids().unwrap().value(), 100.0);
    close(result.inventory.water().value(), 35.0);
    assert_eq!(result.integrated_product.dry_solids().unwrap().value(), 0.0);
    assert_eq!(result.integrated_purge.dry_solids().unwrap().value(), 0.0);
    closes(&result);
}

#[test]
fn bounded_work_cancellation_and_unsupported_empty_start_never_publish_partial_state() {
    let model = standard_model();
    let start = initial();
    let feed = fresh(1.0);
    let mut cfg = settings(1.0);
    cfg.time.maximum_steps = 1;
    assert!(matches!(
        model.advance_inventory(
            &start,
            &feed,
            Time::new(10.0),
            Time::new(2.0),
            &QueryPoint::new(),
            &conditions(),
            recycle(0.75),
            cfg,
            || false
        ),
        Err(TransientCircuitError::Mill(TransientMillError::Budget))
    ));
    let mut polls = 0;
    assert!(matches!(
        model.advance_inventory(
            &start,
            &feed,
            Time::new(10.0),
            Time::new(2.0),
            &QueryPoint::new(),
            &conditions(),
            recycle(0.75),
            settings(1.0),
            || {
                polls += 1;
                polls == 5
            }
        ),
        Err(TransientCircuitError::Mill(TransientMillError::Cancelled))
    ));
    let empty = MillInventory::empty_like(&feed).unwrap();
    assert!(
        model
            .advance_inventory(
                &empty,
                &feed,
                Time::new(10.0),
                Time::new(1.0),
                &QueryPoint::new(),
                &conditions(),
                recycle(0.75),
                settings(1.0),
                || false
            )
            .is_err()
    );
    assert_eq!(start, initial());
    assert_eq!(feed, fresh(1.0));
}
