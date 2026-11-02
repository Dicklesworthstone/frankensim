use super::*;
use crate::comminution::{BondSizeDomain, BondTestContext, BondWorkIndex};
use fs_matdb::Provenance;

fn card() -> BondGrindingCard {
    BondGrindingCard::new(
        BondWorkIndex::from_kwh_per_metric_tonne(15.0).unwrap(),
        BondTestContext {
            specimen: "authored inverse-Bond fixture".into(),
            method: "authored Bond ball-mill test, prepared feed".into(),
            closing_screen: Length::new(150e-6),
            provenance: Provenance {
                source: "synthetic known answer, not site calibration".into(),
                license: "project fixture".into(),
                artifact: None,
            },
        },
        BondSizeDomain::new(
            [Length::new(100e-6), Length::new(0.1)],
            [Length::new(70e-6), Length::new(2_500e-6)],
        )
        .unwrap(),
    )
    .unwrap()
}

#[test]
fn known_energy_and_power_throughput_recover_one_hundred_micrometres() {
    let card = card();
    let original = card.clone();
    // Wi=15 kWh/t, F80=2500 um, P80=100 um gives 12 kWh/t = 43200 J/kg.
    let direct = card
        .product_size_for_specific_energy(
            Length::new(0.0025),
            SpecificGrindingEnergy::new(43_200.0),
        )
        .unwrap();
    assert!((direct.value() - 100e-6).abs() < 3e-19);
    // 12 MW * 0.9 / 250 kg/s = the same specific energy.
    let throughput = card
        .product_size_for_throughput(
            Length::new(0.0025),
            Power::new(12e6),
            Dimensionless::new(0.9),
            MassFlowRate::new(250.0),
        )
        .unwrap();
    assert!((throughput.value() - direct.value()).abs() < 3e-19);
    assert_eq!(card, original);
}

#[test]
fn forward_inverse_round_trip_includes_declared_endpoints_and_tiny_reduction() {
    let card = card();
    for (feed, product) in [
        (0.0025, 70e-6),
        (0.01, 2_500e-6),
        (0.0025, 212e-6),
        (0.001, 600e-6),
        (0.0005, 0.0005_f64.next_down()),
    ] {
        let energy = card
            .specific_energy(Length::new(feed), Length::new(product))
            .unwrap();
        let inverse = card
            .product_size_for_specific_energy(Length::new(feed), energy)
            .unwrap();
        assert!(
            (inverse.value() - product).abs() <= 4e-15 * product,
            "feed={feed}, product={product}, inverted={}",
            inverse.value()
        );
        assert!(inverse.value() < feed);
        if product == 70e-6 || product == 2_500e-6 {
            assert_eq!(inverse.value(), product);
        }
    }
}

#[test]
fn outside_size_energy_window_and_zero_reduction_refuse_without_clamping() {
    let card = card();
    let feed = Length::new(0.01);
    let high = card
        .specific_energy(feed, Length::new(70e-6))
        .unwrap()
        .value();
    let low = card
        .specific_energy(feed, Length::new(2_500e-6))
        .unwrap()
        .value();
    assert!(
        card.product_size_for_specific_energy(feed, SpecificGrindingEnergy::new(high.next_up()))
            .is_err()
    );
    assert!(
        card.product_size_for_specific_energy(feed, SpecificGrindingEnergy::new(low.next_down()))
            .is_err()
    );
    assert!(
        card.product_size_for_specific_energy(
            Length::new(0.2),
            SpecificGrindingEnergy::new(43_200.0)
        )
        .is_err()
    );
    for value in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            card.product_size_for_specific_energy(feed, SpecificGrindingEnergy::new(value))
                .is_err()
        );
    }
    assert!(
        card.product_size_for_specific_energy(
            Length::new(0.0005),
            SpecificGrindingEnergy::new(f64::MIN_POSITIVE)
        )
        .is_err()
    );
}

#[test]
fn throughput_inverse_tracks_power_and_flow_and_rejects_downtime() {
    let card = card();
    let solve = |power, utilization, flow| {
        card.product_size_for_throughput(
            Length::new(0.0025),
            Power::new(power),
            Dimensionless::new(utilization),
            MassFlowRate::new(flow),
        )
    };
    let reference = solve(12e6, 0.9, 250.0).unwrap().value();
    assert!(solve(6e6, 0.9, 250.0).unwrap().value() > reference);
    assert!(solve(12e6, 0.9, 500.0).unwrap().value() > reference);
    assert!(solve(12e6, 0.9, 10.0).is_err()); // Would require a sub-70-um grind.
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(solve(bad, 0.9, 250.0).is_err());
        assert!(solve(12e6, 0.9, bad).is_err());
        assert!(solve(12e6, bad, 250.0).is_err());
    }
    assert!(solve(12e6, 1.01, 250.0).is_err());
}

#[test]
fn scaled_power_energy_preserves_cancelling_extreme_scales() {
    let tiny = f64::from_bits(1);
    // Both naive operation orders lose this representable quotient:
    // power*utilization underflows, while power/flow overflows.
    assert_eq!(positive_product_ratio(1e-15, tiny, tiny).unwrap(), 1e-15);
    // Also retain a useful grinding energy when power/flow overflows.
    assert_eq!(
        positive_product_ratio(43_200.0, tiny, tiny).unwrap(),
        43_200.0
    );
    let card = card();
    let size = card
        .product_size_for_throughput(
            Length::new(0.0025),
            Power::new(43_200.0),
            Dimensionless::new(tiny),
            MassFlowRate::new(tiny),
        )
        .unwrap();
    assert!((size.value() - 100e-6).abs() < 3e-19);
    assert!(positive_product_ratio(f64::MAX, 1.0, tiny).is_err());
    assert!(positive_product_ratio(tiny, tiny, 1.0).is_err());
}
