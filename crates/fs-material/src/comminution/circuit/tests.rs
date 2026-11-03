//! G0 independent closed-loop balances; G1 applicability/traps; G3 scaling.
use super::*;
use crate::comminution::breakage::{
    BreakageColumn, BreakageTestContext, BreakageTransfer, SizeClass,
};
use crate::separation::{CLASSIFIER_LAW_ID, PARTICLE_DIAMETER_AXIS};
use fs_blake3::hash_bytes;
use fs_evidence::ValidityDomain;
use fs_matdb::{ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance};
use fs_qty::{Dims, Length, Pressure, Temperature};

fn source() -> Provenance {
    Provenance {
        source: "authored closed-circuit arithmetic fixture, not measured ore".into(),
        license: "project fixture".into(),
        artifact: Some(hash_bytes(b"closed-circuit-test-v1")),
    }
}
fn mill(
    classes: &[(&str, f64)],
    columns: &[(&str, Vec<(&str, f64)>)],
    actual_box: bool,
) -> SizeRedistributionCard {
    let mut validity = ValidityDomain::unconstrained().with_quantity(
        "temperature",
        QuantitySpec::dimensional(Temperature::DIMS),
        298.15,
        298.15,
    );
    if actual_box {
        validity = validity
            .with_quantity(
                SOLIDS_FRACTION_AXIS,
                QuantitySpec::dimensional(Dims::NONE),
                0.7,
                0.75,
            )
            .with_quantity(
                MILL_DRY_FLOW_AXIS,
                QuantitySpec::dimensional(MassFlowRate::DIMS),
                38.0,
                39.0,
            );
    }
    SizeRedistributionCard::new(
        classes
            .iter()
            .map(|&(id, d)| SizeClass {
                id: id.into(),
                diameter: Length::new(d),
            })
            .collect(),
        columns
            .iter()
            .map(|(id, entries)| BreakageColumn {
                input: (*id).into(),
                transfers: entries
                    .iter()
                    .map(|&(output, f)| BreakageTransfer {
                        output: output.into(),
                        mass_fraction: Dimensionless::new(f),
                    })
                    .collect(),
            })
            .collect(),
        BreakageTestContext {
            specimen: "authored ore".into(),
            method: "explicit congruent per-pass mass transfers".into(),
            provenance: source(),
        },
        validity,
    )
    .unwrap()
}
fn classifier_source(bypass: f64) -> ConstitutiveModelCard {
    ConstitutiveModelCard {
        law: LawId(CLASSIFIER_LAW_ID.into()),
        law_version: 1,
        state_schema_version: 1,
        initial_state: InitialStatePolicy::ZeroInternalState,
        parameters: [
            ("corrected-cut-size", 100e-6, Length::DIMS),
            ("sharpness", 2.0, Dims::NONE),
            ("water-recovery-to-underflow", bypass, Dims::NONE),
        ]
        .into_iter()
        .map(|(name, value, dims)| (name.into(), LawParameter { value, dims }))
        .collect(),
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                PARTICLE_DIAMETER_AXIS,
                QuantitySpec::dimensional(Length::DIMS),
                1e-6,
                0.1,
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
                100000.0,
                100000.0,
            ),
        sources: vec![source().artifact.unwrap()],
        provenance: source(),
    }
}
fn classifier(bypass: f64) -> ClassifierCard {
    ClassifierCard::from_model_card(&classifier_source(bypass)).unwrap()
}
fn copper() -> ElementId {
    ElementId::new("Cu").unwrap()
}
fn stream(rows: &[(&str, f64, f64, f64)], water: f64) -> ParticleStream {
    ParticleStream::new(
        rows.iter()
            .map(|&(id, d, m, cu)| {
                ParticleClass::new(
                    id,
                    Length::new(d),
                    MassFlowRate::new(m),
                    BTreeMap::from([(copper(), MassFlowRate::new(cu))]),
                )
                .unwrap()
            })
            .collect(),
        MassFlowRate::new(water),
    )
    .unwrap()
}
fn conditions() -> (QueryPoint, QueryPoint) {
    (
        QueryPoint::new()
            .with_quantity(
                "temperature",
                QuantitySpec::dimensional(Temperature::DIMS),
                298.15,
            )
            .unwrap(),
        QueryPoint::new()
            .with_quantity(
                "pressure",
                QuantitySpec::dimensional(Pressure::DIMS),
                100000.0,
            )
            .unwrap(),
    )
}
fn settings() -> CircuitSettings {
    CircuitSettings {
        absolute_balance_tolerance: MassFlowRate::new(1e-11),
        relative_balance_tolerance: Dimensionless::new(1e-12),
        maximum_circulating_load: Dimensionless::new(100.0),
    }
}
fn circuit(actual_box: bool) -> ClosedGrindingCircuit {
    ClosedGrindingCircuit::new(
        mill(
            &[("fine", 100e-6), ("coarse", 200e-6)],
            &[
                ("fine", vec![("fine", 1.0)]),
                ("coarse", vec![("fine", 0.5), ("coarse", 0.5)]),
            ],
            actual_box,
        ),
        classifier(0.2),
    )
}
fn feed() -> ParticleStream {
    stream(
        &[("coarse", 200e-6, 10.0, 2.0), ("fine", 100e-6, 2.0, 0.2)],
        12.0,
    )
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-12 * expected.abs().max(1e-300),
        "{actual} != {expected}"
    );
}

#[test]
fn analytic_recycle_closes_each_class_copper_and_water_at_actual_conditions() {
    let fresh = feed();
    let original = fresh.clone();
    let (mc, cc) = conditions();
    let model = circuit(true);
    let result = model.solve(&fresh, &mc, &cc, settings()).unwrap();
    // rcoarse=.95, rfine=.6; coarse half-retained each pass.
    // xcoarse=10/(1-.95*.5)=400/21; xfine=(2+.6*.5*xcoarse)/.4=135/7.
    close(
        result.mill_feed.classes()[0].dry_solids().value(),
        400.0 / 21.0,
    );
    close(
        result.mill_feed.classes()[1].dry_solids().value(),
        135.0 / 7.0,
    );
    close(
        result.mill_feed.classes()[0].elements()[&copper()].value(),
        80.0 / 21.0,
    );
    close(
        result.mill_feed.classes()[1].elements()[&copper()].value(),
        47.0 / 14.0,
    );
    close(result.mill_feed.dry_solids().unwrap().value(), 115.0 / 3.0);
    close(result.circulating_load.value(), 79.0 / 36.0);
    close(result.mill_feed.water().value(), 15.0);
    close(result.recycle.water().value(), 3.0);
    close(result.product.water().value(), 12.0);
    close(result.product.dry_solids().unwrap().value(), 12.0);
    close(result.product.element_flow(&copper()).unwrap().value(), 2.2);
    close(
        result.product.classes()[0].dry_solids().value(),
        10.0 / 21.0,
    );
    close(
        result.product.classes()[1].dry_solids().value(),
        242.0 / 21.0,
    );
    close(
        result.mill_feed.solids_mass_fraction().unwrap().value(),
        23.0 / 32.0,
    );
    // Fresh 50% solids is OUTSIDE the mill's 70–75% support. The admitted
    // circuit correctly uses its solved 71.875% feed, not a guessed start state.
    assert_eq!(fresh.solids_mass_fraction().unwrap().value(), 0.5);
    assert!(result.closure.maximum_absolute_error.value() < 1e-12);
    assert_eq!(fresh, original);
    assert_eq!(
        result
            .product
            .discrete_passing_size(Dimensionless::new(0.8))
            .unwrap(),
        Length::new(100e-6)
    );
    assert_eq!(
        fresh
            .discrete_passing_size(Dimensionless::new(0.8))
            .unwrap(),
        Length::new(200e-6)
    );
}

#[test]
fn direct_solution_agrees_with_independent_geometric_passes() {
    let (mc, cc) = conditions();
    let result = circuit(false).solve(&feed(), &mc, &cc, settings()).unwrap();
    let (mut coarse, mut fine) = (0.0, 0.0);
    for _ in 0..100 {
        let next_coarse = 10.0 + 0.95 * 0.5 * coarse;
        let next_fine = 2.0 + 0.6 * (fine + 0.5 * coarse);
        coarse = next_coarse;
        fine = next_fine;
    }
    close(result.mill_feed.classes()[0].dry_solids().value(), coarse);
    close(result.mill_feed.classes()[1].dry_solids().value(), fine);
}

#[test]
fn small_product_tail_supports_large_circulation_without_iterative_convergence() {
    let model = ClosedGrindingCircuit::new(
        mill(
            &[("coarse", 0.001)],
            &[("coarse", vec![("coarse", 1.0)])],
            false,
        ),
        classifier(0.0),
    );
    let (mc, cc) = conditions();
    let fresh = stream(&[("coarse", 0.001, 1.0, 0.2)], 1.0);
    let mut allowance = settings();
    allowance.maximum_circulating_load = Dimensionless::new(1e40);
    let result = model.solve(&fresh, &mc, &cc, allowance).unwrap();
    // d/d50=10, m=2 => overflow=2^-100, evaluated independently of underflow.
    close(
        result.mill_feed.dry_solids().unwrap().value(),
        2.0_f64.powi(100),
    );
    close(result.product.dry_solids().unwrap().value(), 1.0);
    close(result.product.element_flow(&copper()).unwrap().value(), 0.2);
    assert!(result.circulating_load.value() > 1e30);
    assert!(matches!(
        model.solve(&fresh, &mc, &cc, settings()),
        Err(CircuitError::CirculatingLoadExceeded { .. })
    ));
    assert!(matches!(
        model.solve(
            &fresh.scaled(Dimensionless::new(1e290)).unwrap(),
            &mc,
            &cc,
            allowance
        ),
        Err(CircuitError::Arithmetic(_))
    ));
}

#[test]
fn fed_trap_refuses_but_unfed_closed_size_stays_empty() {
    let model = ClosedGrindingCircuit::new(
        mill(
            &[("fine", 100e-6), ("trapped", 0.1)],
            &[
                ("fine", vec![("fine", 1.0)]),
                ("trapped", vec![("trapped", 1.0)]),
            ],
            false,
        ),
        classifier(0.0),
    );
    let (mc, cc) = conditions();
    let fresh = stream(
        &[("fine", 100e-6, 1.0, 0.1), ("trapped", 0.1, 0.0, 0.0)],
        1.0,
    );
    let result = model.solve(&fresh, &mc, &cc, settings()).unwrap();
    assert_eq!(result.mill_feed.classes()[1].dry_solids().value(), 0.0);
    close(result.product.dry_solids().unwrap().value(), 1.0);
    let trapped = stream(
        &[("fine", 100e-6, 1.0, 0.1), ("trapped", 0.1, 1.0, 0.2)],
        1.0,
    );
    let unchanged = trapped.clone();
    assert!(
        matches!(model.solve(&trapped,&mc,&cc,settings()),Err(CircuitError::TrappedClass{class}) if class=="trapped")
    );
    assert_eq!(trapped, unchanged);
}

#[test]
fn actual_source_domains_and_owned_coordinates_cannot_be_bypassed() {
    let (mc, cc) = conditions();
    let fresh = feed();
    let model = circuit(true);
    assert!(
        model
            .solve(&fresh, &QueryPoint::new(), &cc, settings())
            .is_err()
    );
    assert!(
        model
            .solve(&fresh, &mc, &QueryPoint::new(), settings())
            .is_err()
    );
    let spoof = mc
        .clone()
        .with_quantity(
            MILL_DRY_FLOW_AXIS,
            QuantitySpec::dimensional(MassFlowRate::DIMS),
            38.5,
        )
        .unwrap();
    assert!(matches!(
        model.solve(&fresh, &spoof, &cc, settings()),
        Err(CircuitError::Invalid(_))
    ));
    let spoof = cc.clone().with(SOLIDS_FRACTION_AXIS, 0.5).unwrap();
    assert!(model.solve(&fresh, &mc, &spoof, settings()).is_err());
    // Double all extensive flows preserves composition but leaves the declared
    // mill dry-throughput window; the solver must evaluate this actual state.
    assert!(matches!(
        model.solve(
            &fresh.scaled(Dimensionless::new(2.0)).unwrap(),
            &mc,
            &cc,
            settings()
        ),
        Err(CircuitError::Breakage(BreakageError::OutsideValidity))
    ));
    let mut source = classifier_source(0.2);
    source.validity = source.validity.with_quantity(
        SOLIDS_FRACTION_AXIS,
        QuantitySpec::dimensional(Dims::NONE),
        0.0,
        0.6,
    );
    let narrow = ClosedGrindingCircuit::new(
        model.mill().clone(),
        ClassifierCard::from_model_card(&source).unwrap(),
    );
    assert!(matches!(
        narrow.solve(&fresh, &mc, &cc, settings()),
        Err(CircuitError::Separation(SeparationError::OutsideValidity))
    ));
}

#[test]
fn tiny_coefficient_does_not_erase_representable_physical_fragment_flow() {
    let model = ClosedGrindingCircuit::new(
        mill(
            &[("coarse", 10e-6), ("fine", 1e-6)],
            &[
                ("coarse", vec![("coarse", 1.0), ("fine", f64::from_bits(1))]),
                ("fine", vec![("fine", 1.0)]),
            ],
            false,
        ),
        classifier(0.0),
    );
    let (mc, cc) = conditions();
    let fresh = stream(
        &[("coarse", 10e-6, 1e300, 2e299), ("fine", 1e-6, 0.0, 0.0)],
        1e300,
    );
    let result = model.solve(&fresh, &mc, &cc, settings()).unwrap();
    let transferred = result.mill_feed.classes()[0].dry_solids().value() * f64::from_bits(1);
    assert!(transferred > 0.0);
    close(
        result.product.classes()[1].dry_solids().value(),
        transferred,
    );
    assert!(result.product.classes()[1].elements()[&copper()].value() > 0.0);
    close(
        result.product.dry_solids().unwrap().value(),
        fresh.dry_solids().unwrap().value(),
    );
}

#[test]
fn zero_recycle_and_scaling_preserve_physical_accounts() {
    let (mc, cc) = conditions();
    let model = circuit(false);
    let mut source = classifier_source(0.0);
    source
        .parameters
        .get_mut("corrected-cut-size")
        .unwrap()
        .value = 0.01;
    source.parameters.get_mut("sharpness").unwrap().value = 1000.0;
    let no_recycle = ClosedGrindingCircuit::new(
        model.mill().clone(),
        ClassifierCard::from_model_card(&source).unwrap(),
    )
    .solve(&feed(), &mc, &cc, settings())
    .unwrap();
    assert_eq!(no_recycle.circulating_load.value(), 0.0);
    assert_eq!(no_recycle.mill_feed, feed());
    close(no_recycle.product.dry_solids().unwrap().value(), 12.0);
    let base = model.solve(&feed(), &mc, &cc, settings()).unwrap();
    let scaled = model
        .solve(
            &feed().scaled(Dimensionless::new(3.0)).unwrap(),
            &mc,
            &cc,
            settings(),
        )
        .unwrap();
    close(
        scaled.circulating_load.value(),
        base.circulating_load.value(),
    );
    for (a, b) in [
        (&scaled.mill_feed, &base.mill_feed),
        (&scaled.recycle, &base.recycle),
        (&scaled.product, &base.product),
    ] {
        close(
            a.dry_solids().unwrap().value(),
            3.0 * b.dry_solids().unwrap().value(),
        );
        close(
            a.element_flow(&copper()).unwrap().value(),
            3.0 * b.element_flow(&copper()).unwrap().value(),
        );
        close(a.water().value(), 3.0 * b.water().value());
    }
    let mut invalid = settings();
    invalid.relative_balance_tolerance = Dimensionless::new(f64::NAN);
    assert!(matches!(
        model.solve(&feed(), &mc, &cc, invalid),
        Err(CircuitError::Invalid(_))
    ));
    assert!(
        model
            .solve(
                &feed().scaled(Dimensionless::new(0.0)).unwrap(),
                &mc,
                &cc,
                settings()
            )
            .is_err()
    );
    let wrong = stream(
        &[("coarse", 300e-6, 10.0, 2.0), ("fine", 100e-6, 2.0, 0.2)],
        12.0,
    );
    assert!(matches!(
        model.solve(&wrong, &mc, &cc, settings()),
        Err(CircuitError::Invalid(_))
    ));
}
