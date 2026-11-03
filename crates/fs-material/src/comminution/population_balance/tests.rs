use super::*;
use fs_matdb::Provenance;

fn element(name: &str) -> ElementId {
    ElementId::new(name).unwrap()
}
fn context() -> BreakageTestContext {
    BreakageTestContext {
        specimen: "synthetic three-size ore, not measured testwork".into(),
        method: "authored selection frequencies and daughter mass fractions".into(),
        provenance: Provenance {
            source: "independent analytical balances".into(),
            license: "project fixture".into(),
            artifact: None,
        },
    }
}
fn validity() -> ValidityDomain {
    ValidityDomain::unconstrained()
        .with_quantity(
            MEAN_RESIDENCE_TIME_AXIS,
            QuantitySpec::dimensional(Time::DIMS),
            0.1,
            1000.0,
        )
        .with_quantity(
            MILL_DRY_FLOW_AXIS,
            QuantitySpec::dimensional(MassFlowRate::DIMS),
            1.0,
            100.0,
        )
        .with_quantity(
            SOLIDS_FRACTION_AXIS,
            QuantitySpec::dimensional(Dimensionless::DIMS),
            0.5,
            1.0,
        )
        .with_quantity(
            "rotation-frequency",
            QuantitySpec::dimensional(Frequency::DIMS),
            1.0,
            3.0,
        )
}
fn conditions() -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "rotation-frequency",
            QuantitySpec::dimensional(Frequency::DIMS),
            2.0,
        )
        .unwrap()
}
fn layout() -> Vec<SizeClass> {
    [("medium", 150e-6), ("fine", 75e-6), ("coarse", 300e-6)]
        .into_iter()
        .map(|(id, size)| SizeClass {
            id: id.into(),
            diameter: Length::new(size),
        })
        .collect()
}
fn column(input: &str, rate: f64, daughters: &[(&str, f64)]) -> SelectionBreakageClass {
    SelectionBreakageClass {
        input: input.into(),
        selection_rate: Frequency::new(rate),
        daughters: daughters
            .iter()
            .map(|&(output, fraction)| BreakageTransfer {
                output: output.into(),
                mass_fraction: Dimensionless::new(fraction),
            })
            .collect(),
    }
}
fn kinetics() -> Vec<SelectionBreakageClass> {
    vec![
        column("medium", 0.2, &[("fine", 1.0)]),
        column("fine", 0.0, &[]),
        column("coarse", 0.1, &[("medium", 0.6), ("fine", 0.4)]),
    ]
}
fn model() -> ContinuousPopulationBalanceMill {
    ContinuousPopulationBalanceMill::new(layout(), kinetics(), context(), validity()).unwrap()
}
fn particle(id: &str, size: f64, dry: f64, cu: f64, fe: f64) -> ParticleClass {
    ParticleClass::new(
        id,
        Length::new(size),
        MassFlowRate::new(dry),
        BTreeMap::from([
            (element("Cu"), MassFlowRate::new(cu)),
            (element("Fe"), MassFlowRate::new(fe)),
        ]),
    )
    .unwrap()
}
fn feed() -> ParticleStream {
    ParticleStream::new(
        vec![
            particle("medium", 150e-6, 10.0, 3.0, 1.0),
            particle("coarse", 300e-6, 20.0, 2.0, 4.0),
            particle("fine", 75e-6, 5.0, 0.5, 1.0),
        ],
        MassFlowRate::new(7.0),
    )
    .unwrap()
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 32.0 * f64::EPSILON * expected.abs().max(1.0),
        "actual={actual:.17e}, expected={expected:.17e}"
    );
}
fn class<'a>(stream: &'a ParticleStream, id: &str) -> &'a ParticleClass {
    stream.classes().iter().find(|c| c.id() == id).unwrap()
}

#[test]
fn two_class_analytical_cstr_response_is_not_a_batch_survival_curve() {
    let model = ContinuousPopulationBalanceMill::new(
        vec![
            SizeClass {
                id: "coarse".into(),
                diameter: Length::new(300e-6),
            },
            SizeClass {
                id: "fine".into(),
                diameter: Length::new(75e-6),
            },
        ],
        vec![
            column("coarse", 0.1, &[("fine", 1.0)]),
            column("fine", 0.0, &[]),
        ],
        context(),
        validity(),
    )
    .unwrap();
    let feed = ParticleStream::new(
        vec![
            particle("coarse", 300e-6, 8.0, 1.6, 0.8),
            particle("fine", 75e-6, 2.0, 0.2, 0.4),
        ],
        MassFlowRate::new(3.0),
    )
    .unwrap();
    let solution = model
        .at_residence_time(&feed, Time::new(10.0), &conditions())
        .unwrap();
    close(class(&solution.outlet, "coarse").dry_solids().value(), 4.0);
    close(class(&solution.outlet, "fine").dry_solids().value(), 6.0);
    close(
        class(&solution.outlet, "coarse").elements()[&element("Cu")].value(),
        0.8,
    );
    close(
        class(&solution.outlet, "fine").elements()[&element("Cu")].value(),
        1.0,
    );
    close(solution.dry_solids_holdup.value(), 100.0);
    close(solution.water_holdup.value(), 30.0);
    assert_eq!(solution.outlet.water(), feed.water());
    assert_eq!(solution.model(), &model);
    assert_eq!(solution.model().test_context(), &context());
}

#[test]
fn recursive_multibin_solution_conserves_every_account_and_hold_up() {
    let input = feed();
    let original = input.clone();
    let model = model();
    let solution = model
        .at_residence_time(&input, Time::new(10.0), &conditions())
        .unwrap();
    for (id, dry, cu, fe) in [
        ("coarse", 10.0, 1.0, 2.0),
        ("medium", 16.0 / 3.0, 1.2, 11.0 / 15.0),
        ("fine", 59.0 / 3.0, 3.3, 49.0 / 15.0),
    ] {
        let out = class(&solution.outlet, id);
        close(out.dry_solids().value(), dry);
        close(out.elements()[&element("Cu")].value(), cu);
        close(out.elements()[&element("Fe")].value(), fe);
        let hold = solution.class_holdup.iter().find(|h| h.id == id).unwrap();
        close(hold.dry_solids.value(), 10.0 * dry);
        close(hold.elements[&element("Cu")].value(), 10.0 * cu);
        close(hold.elements[&element("Fe")].value(), 10.0 * fe);
    }
    close(
        solution.outlet.dry_solids().unwrap().value(),
        input.dry_solids().unwrap().value(),
    );
    for name in ["Cu", "Fe"] {
        close(
            solution
                .outlet
                .element_flow(&element(name))
                .unwrap()
                .value(),
            input.element_flow(&element(name)).unwrap().value(),
        );
    }
    assert!(solution.outlet.element_flow(&element("Au")).is_err());
    close(solution.dry_solids_holdup.value(), 350.0);
    close(solution.water_holdup.value(), 70.0);
    assert_eq!(solution.outlet.water(), input.water());
    assert_eq!(input, original);
}

#[test]
fn fixed_hold_up_and_increased_throughput_change_product_size() {
    let model = model();
    let input = feed();
    let doubled = ParticleStream::mix(&[&input, &input]).unwrap();
    let slow = model
        .at_dry_holdup(&input, Mass::new(350.0), &conditions())
        .unwrap();
    let fast = model
        .at_dry_holdup(&doubled, Mass::new(350.0), &conditions())
        .unwrap();
    close(slow.mean_residence_time.value(), 10.0);
    close(fast.mean_residence_time.value(), 5.0);
    close(fast.dry_solids_holdup.value(), 350.0);
    close(
        class(&fast.outlet, "coarse").dry_solids().value(),
        40.0 / 1.5,
    );
    assert!(
        class(&fast.outlet, "fine").dry_solids().value() / 70.0
            < class(&slow.outlet, "fine").dry_solids().value() / 35.0
    );
    assert_eq!(fast.conditions.axes()[MILL_DRY_FLOW_AXIS], 70.0);
    assert_eq!(fast.conditions.axes()[MEAN_RESIDENCE_TIME_AXIS], 5.0);
}

#[test]
fn inactive_classes_and_input_order_retain_the_same_declared_material() {
    let inert = layout().iter().map(|c| column(&c.id, 0.0, &[])).collect();
    let identity =
        ContinuousPopulationBalanceMill::new(layout(), inert, context(), validity()).unwrap();
    assert_eq!(
        identity
            .at_residence_time(&feed(), Time::new(2.0), &conditions())
            .unwrap()
            .outlet,
        feed()
    );
    let mut classes = layout();
    classes.reverse();
    let mut rates = kinetics();
    rates.reverse();
    for c in &mut rates {
        c.daughters.reverse();
    }
    let permuted =
        ContinuousPopulationBalanceMill::new(classes, rates, context(), validity()).unwrap();
    assert_eq!(permuted, model());
}

#[test]
fn source_support_uses_actual_conditions_and_refuses_forged_coordinates() {
    let model = model();
    let input = feed();
    assert_eq!(
        model.at_residence_time(&input, Time::new(1001.0), &conditions()),
        Err(PopulationBalanceError::OutsideValidity)
    );
    assert_eq!(
        model.at_residence_time(&input, Time::new(1.0), &QueryPoint::new()),
        Err(PopulationBalanceError::OutsideValidity)
    );
    let wrong_units = QueryPoint::new()
        .with_quantity(
            "rotation-frequency",
            QuantitySpec::dimensional(Time::DIMS),
            2.0,
        )
        .unwrap();
    assert_eq!(
        model.at_residence_time(&input, Time::new(1.0), &wrong_units),
        Err(PopulationBalanceError::OutsideValidity)
    );
    let too_much = ParticleStream::mix(&[&input, &input, &input]).unwrap();
    assert_eq!(
        model.at_dry_holdup(&too_much, Mass::new(350.0), &conditions()),
        Err(PopulationBalanceError::OutsideValidity)
    );
    for axis in [
        MEAN_RESIDENCE_TIME_AXIS,
        MILL_DRY_FLOW_AXIS,
        SOLIDS_FRACTION_AXIS,
    ] {
        let forged = conditions()
            .with_quantity(axis, QuantitySpec::dimensional(Time::DIMS), 1.0)
            .unwrap();
        assert!(matches!(
            model.at_residence_time(&input, Time::new(1.0), &forged),
            Err(PopulationBalanceError::Invalid(_))
        ));
    }
    let wrong_domain = validity().with_quantity(
        MEAN_RESIDENCE_TIME_AXIS,
        QuantitySpec::dimensional(Frequency::DIMS),
        0.1,
        1000.0,
    );
    assert!(
        ContinuousPopulationBalanceMill::new(layout(), kinetics(), context(), wrong_domain)
            .is_err()
    );
    assert!(
        ContinuousPopulationBalanceMill::new(
            layout(),
            kinetics(),
            context(),
            ValidityDomain::unconstrained()
        )
        .is_err()
    );
}

#[test]
fn invalid_kinetics_layout_and_nonphysical_operating_inputs_refuse() {
    let malformed = [
        column("coarse", -0.1, &[("fine", 1.0)]),
        column("coarse", f64::NAN, &[("fine", 1.0)]),
        column("coarse", 0.1, &[("fine", 0.8)]),
        column("coarse", 0.1, &[("coarse", 1.0)]),
        column("coarse", 0.1, &[("unknown", 1.0)]),
        column("coarse", 0.0, &[("fine", 1.0)]),
    ];
    for bad in malformed {
        let mut rates = kinetics();
        rates[2] = bad;
        assert!(
            ContinuousPopulationBalanceMill::new(layout(), rates, context(), validity()).is_err()
        );
    }
    let mut fine_breaks = kinetics();
    fine_breaks[1] = column("fine", 0.1, &[("medium", 1.0)]);
    assert!(
        ContinuousPopulationBalanceMill::new(layout(), fine_breaks, context(), validity()).is_err()
    );
    let model = model();
    let mut wrong_sizes = layout();
    wrong_sizes[0].diameter = Length::new(200e-6);
    let wrong_model =
        ContinuousPopulationBalanceMill::new(wrong_sizes, kinetics(), context(), validity())
            .unwrap();
    assert!(
        wrong_model
            .at_residence_time(&feed(), Time::new(10.0), &conditions())
            .is_err()
    );
    for value in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            model
                .at_residence_time(&feed(), Time::new(value), &conditions())
                .is_err()
        );
        assert!(
            model
                .at_dry_holdup(&feed(), Mass::new(value), &conditions())
                .is_err()
        );
    }
}
