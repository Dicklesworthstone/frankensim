use super::*;
use fs_qty::{Frequency, QuantitySpec, chemistry::ElementId};

fn context() -> BreakageTestContext {
    BreakageTestContext {
        specimen: "synthetic three-size ore fixture".into(),
        method: "authored congruent single-pass matrix, not a fitted plant law".into(),
        provenance: Provenance {
            source: "independent hand mass balance".into(),
            license: "project fixture".into(),
            artifact: None,
        },
    }
}
fn validity() -> ValidityDomain {
    ValidityDomain::unconstrained().with_quantity(
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
    vec![
        SizeClass {
            id: "medium".into(),
            diameter: Length::new(150e-6),
        },
        SizeClass {
            id: "fine".into(),
            diameter: Length::new(75e-6),
        },
        SizeClass {
            id: "coarse".into(),
            diameter: Length::new(300e-6),
        },
    ]
}
fn column(input: &str, transfers: &[(&str, f64)]) -> BreakageColumn {
    BreakageColumn {
        input: input.into(),
        transfers: transfers
            .iter()
            .map(|&(output, mass_fraction)| BreakageTransfer {
                output: output.into(),
                mass_fraction: Dimensionless::new(mass_fraction),
            })
            .collect(),
    }
}
fn columns() -> Vec<BreakageColumn> {
    vec![
        column("medium", &[("medium", 0.4), ("fine", 0.6)]),
        column(
            "coarse",
            &[("medium", 0.5), ("coarse", 0.25), ("fine", 0.25)],
        ),
        column("fine", &[("fine", 1.0)]),
    ]
}
fn card() -> SizeRedistributionCard {
    SizeRedistributionCard::new(layout(), columns(), context(), validity()).unwrap()
}
fn element(name: &str) -> ElementId {
    ElementId::new(name).unwrap()
}
fn particle(id: &str, diameter: f64, mass: f64, cu: f64, fe: f64) -> ParticleClass {
    ParticleClass::new(
        id,
        Length::new(diameter),
        MassFlowRate::new(mass),
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
        (actual - expected).abs() <= 8.0 * f64::EPSILON * expected.abs(),
        "actual={actual:e}, expected={expected:e}"
    );
}

#[test]
fn hand_balance_conserves_solids_each_element_and_exact_carrier_water() {
    let model = card();
    let input = feed();
    let original = input.clone();
    let output = model.transform(&input, &conditions()).unwrap();
    for (id, dry, cu, fe) in [
        ("coarse", 5.0, 0.5, 1.0),
        ("fine", 16.0, 2.8, 2.6),
        ("medium", 14.0, 2.2, 2.4),
    ] {
        let class = output
            .classes()
            .iter()
            .find(|class| class.id() == id)
            .unwrap();
        close(class.dry_solids().value(), dry);
        close(class.elements()[&element("Cu")].value(), cu);
        close(class.elements()[&element("Fe")].value(), fe);
        assert_eq!(
            class.diameter(),
            input
                .classes()
                .iter()
                .find(|c| c.id() == id)
                .unwrap()
                .diameter()
        );
    }
    close(output.dry_solids().unwrap().value(), 35.0);
    for name in ["Cu", "Fe"] {
        close(
            output.element_flow(&element(name)).unwrap().value(),
            input.element_flow(&element(name)).unwrap().value(),
        );
    }
    assert_eq!(
        output.water().value().to_bits(),
        input.water().value().to_bits()
    );
    assert_eq!(input, original);
    assert_eq!(model.test_context(), &context());
    assert_eq!(model.validity(), &validity());
    assert_eq!(
        model
            .classes()
            .iter()
            .map(|c| c.id.as_str())
            .collect::<Vec<_>>(),
        ["coarse", "fine", "medium"]
    );
    assert_eq!(model.fraction("coarse", "fine").unwrap().value(), 0.0);
    assert_eq!(model.retained_fraction("coarse").unwrap().value(), 0.25);
    assert_eq!(model.fragmented_fraction("coarse").unwrap().value(), 0.75);
    assert!(model.fraction("unknown", "coarse").is_err());
    assert!(model.retained_fraction("unknown").is_err());
}

#[test]
fn identity_and_input_order_are_preserved_without_invented_assays() {
    let classes = layout();
    let identity = classes
        .iter()
        .map(|c| column(&c.id, &[(&c.id, 1.0)]))
        .collect();
    let model = SizeRedistributionCard::new(classes, identity, context(), validity()).unwrap();
    assert_eq!(model.transform(&feed(), &conditions()).unwrap(), feed());
    assert_eq!(model.fragmented_fraction("medium").unwrap().value(), 0.0);
    let mut shuffled_classes = layout();
    shuffled_classes.reverse();
    let mut shuffled_columns = columns();
    shuffled_columns.reverse();
    for column in &mut shuffled_columns {
        column.transfers.reverse();
    }
    let shuffled =
        SizeRedistributionCard::new(shuffled_classes, shuffled_columns, context(), validity())
            .unwrap();
    assert_eq!(shuffled, card());
    assert!(
        shuffled
            .transform(&feed(), &conditions())
            .unwrap()
            .element_flow(&element("Au"))
            .is_err()
    );
}

#[test]
fn tiny_fragmented_share_and_many_small_inputs_survive_compensated_accumulation() {
    let tiny = SizeRedistributionCard::new(
        layout(),
        vec![
            column("coarse", &[("coarse", 1.0), ("fine", 1e-30)]),
            column("fine", &[("fine", 1.0)]),
            column("medium", &[("medium", 1.0)]),
        ],
        context(),
        validity(),
    )
    .unwrap();
    assert_eq!(tiny.retained_fraction("coarse").unwrap().value(), 1.0);
    assert_eq!(tiny.fragmented_fraction("coarse").unwrap().value(), 1e-30);
    let input = ParticleStream::new(
        vec![
            particle("coarse", 300e-6, 1.0, 0.25, 0.125),
            particle("medium", 150e-6, 0.0, 0.0, 0.0),
            particle("fine", 75e-6, 0.0, 0.0, 0.0),
        ],
        MassFlowRate::new(1.0),
    )
    .unwrap();
    let output = tiny.transform(&input, &conditions()).unwrap();
    let fine = output.classes().iter().find(|c| c.id() == "fine").unwrap();
    close(fine.dry_solids().value(), 1e-30);
    close(fine.elements()[&element("Cu")].value(), 0.25e-30);

    let mut classes = vec![SizeClass {
        id: "z-fine".into(),
        diameter: Length::new(75e-6),
    }];
    let mut matrix = vec![column("z-fine", &[("z-fine", 1.0)])];
    let mut feed_classes = vec![particle("z-fine", 75e-6, 0.0, 0.0, 0.0)];
    for i in 0..127 {
        let id = format!("c{i:03}");
        let mass = if i == 0 { 1.0 } else { 1e-16 };
        classes.push(SizeClass {
            id: id.clone(),
            diameter: Length::new(200e-6),
        });
        matrix.push(column(&id, &[("z-fine", 1.0)]));
        feed_classes.push(particle(&id, 200e-6, mass, 0.25 * mass, 0.125 * mass));
    }
    let model = SizeRedistributionCard::new(classes, matrix, context(), validity()).unwrap();
    let input = ParticleStream::new(feed_classes, MassFlowRate::new(1.0)).unwrap();
    let output = model.transform(&input, &conditions()).unwrap();
    assert_eq!(
        output.classes().last().unwrap().dry_solids().value(),
        1.0 + 126.0e-16
    );
    assert_eq!(output.dry_solids().unwrap(), input.dry_solids().unwrap());
    assert_eq!(
        output.element_flow(&element("Cu")).unwrap(),
        input.element_flow(&element("Cu")).unwrap()
    );
}

#[test]
fn malformed_columns_and_nonbreaking_relabels_are_refused() {
    for bad in [-0.1, 1.1, f64::NAN, f64::INFINITY, 0.5] {
        let mut matrix = columns();
        matrix[2] = column("fine", &[("fine", bad)]);
        assert!(SizeRedistributionCard::new(layout(), matrix, context(), validity()).is_err());
    }
    for bad in [
        column("fine", &[("fine", 0.5), ("fine", 0.5)]),
        column("fine", &[("unknown", 1.0)]),
        column("fine", &[("coarse", 1.0)]),
        column("coarse", &[("coarse", 1.0)]),
    ] {
        let mut matrix = columns();
        matrix[2] = bad;
        assert!(SizeRedistributionCard::new(layout(), matrix, context(), validity()).is_err());
    }
    let mut equal_sizes = layout();
    equal_sizes[0].diameter = Length::new(300e-6); // coarse -> medium becomes a relabel.
    assert!(SizeRedistributionCard::new(equal_sizes, columns(), context(), validity()).is_err());
    let mut duplicates = layout();
    duplicates[0].id = "coarse".into();
    assert!(SizeRedistributionCard::new(duplicates, columns(), context(), validity()).is_err());
    assert!(
        SizeRedistributionCard::new(layout(), columns()[..2].to_vec(), context(), validity())
            .is_err()
    );
    let too_many = vec![
        SizeClass {
            id: "x".into(),
            diameter: Length::new(1.0)
        };
        MAX_SIZE_CLASSES + 1
    ];
    assert!(SizeRedistributionCard::new(too_many, vec![], context(), validity()).is_err());
}

#[test]
fn source_conditions_and_feed_layout_cannot_be_bypassed() {
    let model = card();
    assert_eq!(
        model.transform(&feed(), &QueryPoint::new()),
        Err(BreakageError::OutsideValidity)
    );
    let untyped = QueryPoint::new().with("rotation-frequency", 2.0).unwrap();
    assert_eq!(
        model.transform(&feed(), &untyped),
        Err(BreakageError::OutsideValidity)
    );
    let outside = conditions()
        .with("rotation-frequency", 3.0_f64.next_up())
        .unwrap();
    assert_eq!(
        model.transform(&feed(), &outside),
        Err(BreakageError::OutsideValidity)
    );
    let mut actual = feed().classes().to_vec();
    actual[0] = particle("coarse", 300e-6_f64.next_up(), 20.0, 2.0, 4.0);
    let changed_size = ParticleStream::new(actual, MassFlowRate::new(7.0)).unwrap();
    assert!(model.transform(&changed_size, &conditions()).is_err());
    let missing =
        ParticleStream::new(feed().classes()[1..].to_vec(), MassFlowRate::new(7.0)).unwrap();
    assert!(model.transform(&missing, &conditions()).is_err());
    let mut absent = context();
    absent.provenance.source.clear();
    assert!(SizeRedistributionCard::new(layout(), columns(), absent, validity()).is_err());
    let untyped_domain = ValidityDomain::unconstrained().with("rotation-frequency", 1.0, 3.0);
    assert!(SizeRedistributionCard::new(layout(), columns(), context(), untyped_domain).is_err());
}

#[test]
fn positive_mass_or_trace_element_underflow_refuses_atomically() {
    let model = card();
    for (mass, cu) in [(f64::from_bits(1), 0.0), (1.0, f64::from_bits(1))] {
        let input = ParticleStream::new(
            vec![
                particle("coarse", 300e-6, mass, cu, 0.0),
                particle("medium", 150e-6, 0.0, 0.0, 0.0),
                particle("fine", 75e-6, 0.0, 0.0, 0.0),
            ],
            MassFlowRate::new(1.0),
        )
        .unwrap();
        let original = input.clone();
        assert!(matches!(
            model.transform(&input, &conditions()),
            Err(BreakageError::Arithmetic(_))
        ));
        assert_eq!(input, original);
    }
}
