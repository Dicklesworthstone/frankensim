use super::*;
use crate::comminution::{BondSizeDomain, BondTestContext, BondWorkIndex};
use fs_matdb::Provenance;

fn feed(name: &str, mass: f64, wi: f64, f80: f64, p80: [f64; 2]) -> GrindingFeed {
    GrindingFeed {
        card: BondGrindingCard::new(
            BondWorkIndex::from_kwh_per_metric_tonne(wi).unwrap(),
            BondTestContext {
                specimen: name.into(),
                method: "authored arithmetic fixture, not an ore test".into(),
                closing_screen: Length::new(150e-6),
                provenance: Provenance {
                    source: "authored additive-energy scenario".into(),
                    license: "project fixture".into(),
                    artifact: None,
                },
            },
            BondSizeDomain::new(
                [Length::new(200e-6), Length::new(0.01)],
                p80.map(Length::new),
            )
            .unwrap(),
        )
        .unwrap(),
        dry_mass: Mass::new(mass),
        feed_p80: Length::new(f80),
    }
}

fn ordinary_feeds() -> Vec<GrindingFeed> {
    vec![
        feed("A", 25.0, 10.0, 2_500e-6, [70e-6, 150e-6]),
        feed("B", 75.0, 20.0, 400e-6, [90e-6, 200e-6]),
    ]
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 4e-15 * expected.abs(),
        "{actual} != {expected}"
    );
}

#[test]
fn blends_specific_energy_before_converting_power_to_throughput() {
    let inputs = ordinary_feeds();
    let blend = AdditiveGrindingBlend::new(inputs.clone()).unwrap();
    assert_eq!(blend.feeds(), inputs);
    let product = Length::new(100e-6);
    // Independently hand calculated in micrometres:
    // A: 10*10*(1/10-1/50) = 8 kWh/t.
    // B: 10*20*(1/10-1/20) = 10 kWh/t.
    // 25% A + 75% B = 9.5 kWh/t = 34200 J/kg.
    close(blend.specific_energy(product).unwrap().value(), 34_200.0);
    let throughput = blend
        .dry_throughput(product, Power::new(1e6), Dimensionless::new(0.9))
        .unwrap()
        .value();
    close(throughput, 900_000.0 / 34_200.0);
    let average_throughputs = 0.25 * 900_000.0 / 28_800.0 + 0.75 * 900_000.0 / 36_000.0;
    assert!((throughput - average_throughputs).abs() > 0.2);
    // A mean Wi and mean F80 would give a different, unsupported energy.
    let naive = feed("naive", 100.0, 17.5, 925e-6, [70e-6, 200e-6]);
    let naive_energy = naive
        .card
        .specific_energy(naive.feed_p80, product)
        .unwrap()
        .value();
    assert!((naive_energy - 34_200.0).abs() > 1_000.0);
    assert_eq!(
        blend
            .dry_throughput(product, Power::new(1e6), Dimensionless::new(0.0))
            .unwrap()
            .value(),
        0.0
    );
}

#[test]
fn common_domain_retains_each_source_limit_and_actual_feed_size() {
    let blend = AdditiveGrindingBlend::new(ordinary_feeds()).unwrap();
    assert_eq!(
        blend.product_bounds(),
        [Length::new(90e-6), Length::new(150e-6)]
    );
    assert!(blend.specific_energy(Length::new(80e-6)).is_err());
    assert!(blend.specific_energy(Length::new(175e-6)).is_err());
    assert!(blend.specific_energy(Length::new(f64::NAN)).is_err());
    assert!(
        AdditiveGrindingBlend::new(vec![
            feed("A", 1.0, 10.0, 2_500e-6, [70e-6, 100e-6]),
            feed("B", 1.0, 20.0, 400e-6, [150e-6, 200e-6]),
        ])
        .is_err()
    );
    let mut invalid = ordinary_feeds();
    invalid[0].feed_p80 = Length::new(0.02);
    assert!(AdditiveGrindingBlend::new(invalid).is_err());
    let one =
        AdditiveGrindingBlend::new(vec![feed("A", 1.0, 10.0, 200e-6, [70e-6, 300e-6])]).unwrap();
    assert!(one.product_bounds()[1] < Length::new(200e-6));
    assert!(one.specific_energy(Length::new(200e-6)).is_err());
}

#[test]
fn mass_scaling_feed_order_and_split_lots_preserve_physical_energy() {
    let product = Length::new(100e-6);
    let original = ordinary_feeds();
    let expected = AdditiveGrindingBlend::new(original.clone())
        .unwrap()
        .specific_energy(product)
        .unwrap()
        .value();
    let mut reversed = original.clone();
    reversed.reverse();
    close(
        AdditiveGrindingBlend::new(reversed)
            .unwrap()
            .specific_energy(product)
            .unwrap()
            .value(),
        expected,
    );
    let mut scaled = original.clone();
    for input in &mut scaled {
        input.dry_mass = Mass::new(input.dry_mass.value() * 2e306);
    }
    // Raw m*E and sum(m) overflow; the specific-energy mean remains finite.
    close(
        AdditiveGrindingBlend::new(scaled)
            .unwrap()
            .specific_energy(product)
            .unwrap()
            .value(),
        expected,
    );
    let mut split = original;
    let mut extra = split[1].clone();
    split[1].dry_mass = Mass::new(50.0);
    extra.dry_mass = Mass::new(25.0);
    split.push(extra);
    close(
        AdditiveGrindingBlend::new(split)
            .unwrap()
            .specific_energy(product)
            .unwrap()
            .value(),
        expected,
    );
}

#[test]
fn absent_invalid_or_unrepresentable_mass_contributions_refuse() {
    assert!(AdditiveGrindingBlend::new(Vec::new()).is_err());
    for mass in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        let mut inputs = ordinary_feeds();
        inputs[0].dry_mass = Mass::new(mass);
        assert!(AdditiveGrindingBlend::new(inputs).is_err());
    }
    let mut inputs = ordinary_feeds();
    inputs[0].dry_mass = Mass::new(f64::from_bits(1));
    inputs[1].dry_mass = Mass::new(f64::MAX);
    assert!(AdditiveGrindingBlend::new(inputs).is_err());
}
