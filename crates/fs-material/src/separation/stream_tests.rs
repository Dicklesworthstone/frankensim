use super::*;

fn copper() -> ElementId {
    ElementId::new("Cu").unwrap()
}
fn iron() -> ElementId {
    ElementId::new("Fe").unwrap()
}

fn class(id: &str, diameter: f64, dry: f64, copper_flow: f64) -> ParticleClass {
    ParticleClass::new(
        id,
        Length::new(diameter),
        MassFlowRate::new(dry),
        BTreeMap::from([
            (copper(), MassFlowRate::new(copper_flow)),
            (iron(), MassFlowRate::new(0.0)),
        ]),
    )
    .unwrap()
}

fn feed() -> ParticleStream {
    ParticleStream::new(
        vec![
            class("fine", 100e-6, 20.0, 2.0),
            class("coarse", 1_000e-6, 30.0, 9.0),
        ],
        MassFlowRate::new(40.0),
    )
    .unwrap()
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 8.0 * f64::EPSILON * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn mixing_conserves_explicit_classes_elements_water_and_split_recombination() {
    let first = feed();
    let second = ParticleStream::new(
        vec![
            class("coarse", 1_000e-6, 7.0, 1.0),
            class("fine", 100e-6, 3.0, 0.5),
        ],
        MassFlowRate::new(5.0),
    )
    .unwrap();
    let combined = ParticleStream::mix(&[&first, &second]).unwrap();
    assert_eq!(
        combined
            .classes()
            .iter()
            .map(ParticleClass::id)
            .collect::<Vec<_>>(),
        ["coarse", "fine"]
    );
    assert_eq!(combined.classes()[0].dry_solids().value(), 37.0);
    assert_eq!(combined.classes()[1].dry_solids().value(), 23.0);
    assert_eq!(combined.water().value(), 45.0);
    assert_eq!(combined.element_flow(&copper()).unwrap().value(), 12.5);
    assert_eq!(combined.element_flow(&iron()).unwrap().value(), 0.0);
    assert!(
        combined
            .element_flow(&ElementId::new("Au").unwrap())
            .is_err()
    );
    assert_eq!(combined, ParticleStream::mix(&[&second, &first]).unwrap());
    let split = first
        .split(
            &[
                (
                    "fine".into(),
                    Partition::from_product(Dimensionless::new(0.3)).unwrap(),
                ),
                (
                    "coarse".into(),
                    Partition::from_product(Dimensionless::new(0.8)).unwrap(),
                ),
            ],
            Partition::from_product(Dimensionless::new(0.1)).unwrap(),
        )
        .unwrap();
    let recombined = ParticleStream::mix(&[&split.reject, &split.product]).unwrap();
    close(recombined.dry_solids().unwrap().value(), 50.0);
    close(recombined.element_flow(&copper()).unwrap().value(), 11.0);
    close(recombined.water().value(), 40.0);
    assert_eq!(first, feed());
}

#[test]
fn mixing_rejects_layout_size_and_measured_zero_versus_missing_mismatches() {
    let first = feed();
    assert!(ParticleStream::mix(&[]).is_err());
    let short =
        ParticleStream::new(vec![first.classes()[0].clone()], MassFlowRate::new(0.0)).unwrap();
    assert!(ParticleStream::mix(&[&first, &short]).is_err());
    let renamed = ParticleStream::new(
        vec![
            class("other", 100e-6, 20.0, 2.0),
            first.classes()[0].clone(),
        ],
        MassFlowRate::new(0.0),
    )
    .unwrap();
    assert!(ParticleStream::mix(&[&first, &renamed]).is_err());
    let resized = ParticleStream::new(
        vec![class("fine", 200e-6, 20.0, 2.0), first.classes()[0].clone()],
        MassFlowRate::new(0.0),
    )
    .unwrap();
    assert!(ParticleStream::mix(&[&first, &resized]).is_err());
    let missing_iron = ParticleStream::new(
        first
            .classes()
            .iter()
            .map(|class| {
                ParticleClass::new(
                    class.id(),
                    class.diameter(),
                    MassFlowRate::new(0.0),
                    BTreeMap::from([(copper(), MassFlowRate::new(0.0))]),
                )
                .unwrap()
            })
            .collect(),
        MassFlowRate::new(0.0),
    )
    .unwrap();
    assert!(ParticleStream::mix(&[&first, &missing_iron]).is_err());
}

#[test]
fn scaling_preserves_layout_and_zero_assays_without_normalizing_composition() {
    let original = feed();
    let doubled = original.scaled(Dimensionless::new(2.0)).unwrap();
    assert_eq!(doubled.dry_solids().unwrap().value(), 100.0);
    assert_eq!(doubled.water().value(), 80.0);
    assert_eq!(doubled.element_flow(&copper()).unwrap().value(), 22.0);
    close(
        doubled.element_grade(&copper()).unwrap().unwrap().value(),
        0.22,
    );
    let zero = original.scaled(Dimensionless::new(0.0)).unwrap();
    assert_eq!(zero.classes().len(), original.classes().len());
    for (empty, source) in zero.classes().iter().zip(original.classes()) {
        assert_eq!(empty.id(), source.id());
        assert_eq!(empty.diameter(), source.diameter());
        assert_eq!(
            empty.elements().keys().collect::<Vec<_>>(),
            source.elements().keys().collect::<Vec<_>>()
        );
        assert!(empty.elements().values().all(|flow| flow.value() == 0.0));
    }
    assert_eq!(zero.element_grade(&copper()).unwrap(), None);
    assert!(zero.element_flow(&ElementId::new("Au").unwrap()).is_err());
    assert_eq!(original, ParticleStream::mix(&[&zero, &original]).unwrap());
}

#[test]
fn compensated_mixing_and_nonzero_arithmetic_refusals() {
    let one = |mass| {
        ParticleStream::new(vec![class("only", 1e-3, mass, 0.0)], MassFlowRate::new(0.0)).unwrap()
    };
    let huge = one(1e16);
    let small = one(1.0);
    let recovered = ParticleStream::mix(&[&huge, &small, &small]).unwrap();
    assert_eq!(recovered.dry_solids().unwrap().value(), 1e16 + 2.0);
    assert_eq!(
        recovered,
        ParticleStream::mix(&[&small, &huge, &small]).unwrap()
    );
    let largest = one(f64::MAX);
    assert!(ParticleStream::mix(&[&largest, &largest]).is_err());
    assert!(largest.scaled(Dimensionless::new(2.0)).is_err());
    let tiniest = one(f64::from_bits(1));
    assert!(tiniest.scaled(Dimensionless::new(0.5)).is_err());
    let trace = ParticleStream::new(
        vec![class("only", 1e-3, 1.0, f64::from_bits(1))],
        MassFlowRate::new(1.0),
    )
    .unwrap();
    assert!(trace.scaled(Dimensionless::new(0.5)).is_err());
    let water_trace = ParticleStream::new(
        vec![class("only", 1e-3, 1.0, 0.0)],
        MassFlowRate::new(f64::from_bits(1)),
    )
    .unwrap();
    assert!(water_trace.scaled(Dimensionless::new(0.5)).is_err());
    for factor in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(small.scaled(Dimensionless::new(factor)).is_err());
    }
}

#[test]
fn discrete_quantiles_group_equal_sizes_ignore_water_and_do_not_interpolate() {
    let stream = ParticleStream::new(
        vec![
            class("a-coarse", 1e-3, 20.0, 0.0),
            class("z-fine", 100e-6, 20.0, 0.0),
            class("b-mid", 200e-6, 30.0, 0.0),
            class("c-mid", 200e-6, 30.0, 0.0),
            class("zero-small", 10e-6, 0.0, 0.0),
            class("zero-large", 1.0, 0.0, 0.0),
        ],
        MassFlowRate::new(1e9),
    )
    .unwrap();
    for (fraction, size) in [
        (0.01, 100e-6),
        (0.2, 100e-6),
        (0.2_f64.next_up(), 200e-6),
        (0.5, 200e-6),
        (0.8, 200e-6),
        (0.8_f64.next_up(), 1e-3),
        (1.0, 1e-3),
    ] {
        assert_eq!(
            stream
                .discrete_passing_size(Dimensionless::new(fraction))
                .unwrap(),
            Length::new(size)
        );
    }
    assert_eq!(
        stream
            .scaled(Dimensionless::new(3.0))
            .unwrap()
            .discrete_passing_size(Dimensionless::new(0.8))
            .unwrap(),
        Length::new(200e-6)
    );
}

#[test]
fn discrete_quantile_boundaries_preserve_tiny_positive_tails() {
    let stream = ParticleStream::new(
        vec![
            class("fine", 100e-6, 1.0, 0.0),
            class("coarse", 1e-3, f64::from_bits(1), 0.0),
        ],
        MassFlowRate::new(0.0),
    )
    .unwrap();
    assert_eq!(
        stream
            .discrete_passing_size(Dimensionless::new(1.0))
            .unwrap(),
        Length::new(1e-3)
    );
    let tiny = ParticleStream::new(
        vec![
            class("zero", 1e-6, 0.0, 0.0),
            class("positive", 1e-3, f64::from_bits(1), 0.0),
        ],
        MassFlowRate::new(0.0),
    )
    .unwrap();
    assert_eq!(
        tiny.discrete_passing_size(Dimensionless::new(f64::from_bits(1)))
            .unwrap(),
        Length::new(1e-3)
    );
    assert!(
        stream
            .scaled(Dimensionless::new(0.0))
            .unwrap()
            .discrete_passing_size(Dimensionless::new(0.8))
            .is_err()
    );
    for fraction in [0.0, -0.1, 1.01, f64::NAN, f64::INFINITY] {
        assert!(
            stream
                .discrete_passing_size(Dimensionless::new(fraction))
                .is_err()
        );
    }
}
