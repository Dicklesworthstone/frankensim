//! G0 known answers and refusal boundaries; G3 unit and monotonicity checks.
//! Authored values verify the formula, not a site's plant calibration.
use super::*;

fn context() -> BondTestContext {
    BondTestContext {
        specimen: "synthetic-composite-A".into(),
        method: "authored Bond ball mill fixture; prepared feed".into(),
        closing_screen: Length::new(150.0e-6),
        provenance: Provenance {
            source: "authored test fixture, not a measured ore".into(),
            license: "project fixture".into(),
            artifact: None,
        },
    }
}

fn domain() -> BondSizeDomain {
    BondSizeDomain::new(
        [Length::new(100.0e-6), Length::new(0.1)],
        [Length::new(70.0e-6), Length::new(2_500.0e-6)],
    )
    .unwrap()
}

fn card(index: BondWorkIndex) -> BondGrindingCard {
    BondGrindingCard::new(index, context(), domain()).unwrap()
}

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs(),
        "{actual} != {expected}"
    );
}

#[test]
fn bond_known_answer_and_power_to_dry_tonnage() {
    // Independently hand-calculated: 10*15*(1/sqrt(100)-1/sqrt(2500)) = 12 kWh/t.
    let law = card(BondWorkIndex::from_kwh_per_metric_tonne(15.0).unwrap());
    let feed = Length::new(0.0025);
    let product = Length::new(0.0001);
    close(
        law.specific_energy(feed, product).unwrap().value(),
        43_200.0,
        2e-15,
    );
    // 12,000 kW * 0.9 / 12 kWh/t = 900 dry t/h = 250 kg/s.
    close(
        law.dry_throughput(feed, product, Power::new(12.0e6), Dimensionless::new(0.9))
            .unwrap()
            .value(),
        250.0,
        2e-15,
    );
    assert_eq!(
        law.dry_throughput(feed, product, Power::new(12.0e6), Dimensionless::new(0.0))
            .unwrap()
            .value(),
        0.0
    );
}

#[test]
fn published_gmg_generic_circuit_energy_to_printed_precision() {
    // GMG01-MP-2021 §4.2.1, Table 1: Wi=16.1 kWh/t, F80=2500 µm,
    // intermediate=1000 µm, P80=212 µm. Printed stage energies 1.87 and 5.97.
    // Context is an authored use of those published inputs, not a real lab report.
    let law = card(BondWorkIndex::from_kwh_per_metric_tonne(16.1).unwrap());
    let first = law
        .specific_energy(Length::new(0.0025), Length::new(0.001))
        .unwrap()
        .value()
        / 3_600.0;
    let second = law
        .specific_energy(Length::new(0.001), Length::new(0.000212))
        .unwrap()
        .value()
        / 3_600.0;
    assert!((first - 1.87).abs() < 0.005);
    assert!((second - 5.97).abs() < 0.005);
    assert!((first + second - 7.84).abs() < 0.01);
    let whole = law
        .specific_energy(Length::new(0.0025), Length::new(0.000212))
        .unwrap()
        .value()
        / 3_600.0;
    close(first + second, whole, 2e-15);
}

#[test]
fn exact_short_ton_conversion_preserves_same_physical_work_index() {
    let metric = BondWorkIndex::from_kwh_per_metric_tonne(15.0).unwrap();
    let short = BondWorkIndex::from_kwh_per_short_ton(15.0 * 907.18474 / 1_000.0).unwrap();
    close(
        short.specific_energy().value(),
        metric.specific_energy().value(),
        2e-15,
    );
    assert_eq!(short.source_basis(), WorkIndexMassBasis::ShortTon);
    assert_eq!(metric.source_basis(), WorkIndexMassBasis::MetricTonne);
    let wrong_same_literal = BondWorkIndex::from_kwh_per_short_ton(15.0).unwrap();
    close(
        wrong_same_literal.specific_energy().value() / metric.specific_energy().value(),
        1_000.0 / 907.18474,
        2e-15,
    );
    let f = Length::new(0.0025);
    let p = Length::new(0.0001);
    close(
        card(short).specific_energy(f, p).unwrap().value(),
        card(metric).specific_energy(f, p).unwrap().value(),
        2e-15,
    );
}

#[test]
fn energy_monotonicity_and_near_equal_sizes() {
    let law = card(BondWorkIndex::from_kwh_per_metric_tonne(15.0).unwrap());
    let at = |f, p| {
        law.specific_energy(Length::new(f), Length::new(p))
            .unwrap()
            .value()
    };
    assert!(at(0.004, 0.0001) > at(0.0025, 0.0001));
    assert!(at(0.0025, 0.0001) > at(0.0025, 0.0002));
    let product = 0.0001_f64;
    let feed = f64::from_bits(product.to_bits() + 1);
    assert!(at(feed, product) > 0.0);
    let result = at(0.0025, 0.0001).to_bits();
    assert_eq!(result, at(0.0025, 0.0001).to_bits());
}

#[test]
fn size_validity_and_nonphysical_values_refuse() {
    let law = card(BondWorkIndex::from_kwh_per_metric_tonne(15.0).unwrap());
    assert!(matches!(
        law.specific_energy(Length::new(0.2), Length::new(0.0001)),
        Err(ComminutionError::OutsideValidity { axis: "F80", .. })
    ));
    assert!(matches!(
        law.specific_energy(Length::new(0.0025), Length::new(69.0e-6)),
        Err(ComminutionError::OutsideValidity { axis: "P80", .. })
    ));
    for (f, p) in [
        (0.001, 0.001),
        (0.0001, 0.001),
        (0.0, 0.0001),
        (-1.0, 0.0001),
        (f64::NAN, 0.0001),
        (0.0025, f64::INFINITY),
    ] {
        assert!(law.specific_energy(Length::new(f), Length::new(p)).is_err());
    }
    assert!(
        BondSizeDomain::new(
            [Length::new(0.0025); 2],
            [Length::new(60.0e-6), Length::new(100.0e-6)]
        )
        .is_err()
    );
    assert!(
        BondSizeDomain::new(
            [Length::new(0.0025); 2],
            [Length::new(100.0e-6), Length::new(70.0e-6)]
        )
        .is_err()
    );
    assert!(BondSizeDomain::new([Length::new(0.0001); 2], [Length::new(0.001); 2]).is_err());
    for value in [0.0, -1.0, f64::NAN, f64::INFINITY, f64::MAX] {
        assert!(BondWorkIndex::from_kwh_per_metric_tonne(value).is_err());
        assert!(BondWorkIndex::from_kwh_per_short_ton(value).is_err());
    }
}

#[test]
fn source_and_closing_screen_are_required_and_retained() {
    let wi = BondWorkIndex::from_kwh_per_metric_tonne(15.0).unwrap();
    let law = card(wi);
    assert_eq!(law.test_context(), &context());
    for field in 0..5 {
        let mut test = context();
        match field {
            0 => test.specimen = " ".into(),
            1 => test.method.clear(),
            2 => test.provenance.source.clear(),
            3 => test.provenance.license.clear(),
            _ => test.closing_screen = Length::new(f64::NAN),
        }
        assert!(BondGrindingCard::new(wi, test, domain()).is_err());
    }
}

#[test]
fn invalid_power_energy_and_utilization_refuse() {
    for value in [-1.0, 0.0, f64::NAN, f64::INFINITY] {
        assert!(
            throughput_from_net_power(
                SpecificGrindingEnergy::new(value),
                Power::new(1_000.0),
                Dimensionless::new(1.0)
            )
            .is_err()
        );
        assert!(
            throughput_from_net_power(
                SpecificGrindingEnergy::new(1_000.0),
                Power::new(value),
                Dimensionless::new(1.0)
            )
            .is_err()
        );
    }
    for value in [-0.1, 1.01, f64::NAN, f64::INFINITY] {
        assert!(
            throughput_from_net_power(
                SpecificGrindingEnergy::new(1_000.0),
                Power::new(1_000.0),
                Dimensionless::new(value)
            )
            .is_err()
        );
    }
    assert!(
        throughput_from_net_power(
            SpecificGrindingEnergy::new(f64::MIN_POSITIVE),
            Power::new(f64::MAX),
            Dimensionless::new(1.0)
        )
        .is_err()
    );
}

#[test]
fn critical_speed_balances_centripetal_acceleration_at_top() {
    let d = 4.1;
    let ball = 0.1;
    let gravity = 9.80665;
    let n = ball_mill_critical_speed(
        Length::new(d),
        Length::new(ball),
        Acceleration::new(gravity),
    )
    .unwrap()
    .value();
    let omega = n * core::f64::consts::TAU;
    close(omega * omega * (d - ball) / 2.0, gravity, 2e-15);
    // Independent constant: 60/(2*pi)*sqrt(2*9.80665/4).
    assert!((n * 60.0 - 21.145).abs() < 0.001);
    let larger = ball_mill_critical_speed(
        Length::new(8.1),
        Length::new(ball),
        Acceleration::new(gravity),
    )
    .unwrap()
    .value();
    close(n / larger, 2.0_f64.sqrt(), 2e-15);
}

#[test]
fn critical_speed_requires_finite_physical_clearance_and_gravity() {
    for (d, b, g) in [
        (0.0, 0.1, 9.8),
        (4.0, 0.0, 9.8),
        (4.0, 4.0, 9.8),
        (4.0, 5.0, 9.8),
        (4.0, 0.1, 0.0),
        (f64::NAN, 0.1, 9.8),
        (4.0, 0.1, f64::INFINITY),
    ] {
        assert!(
            ball_mill_critical_speed(Length::new(d), Length::new(b), Acceleration::new(g)).is_err()
        );
    }
    // Finite output despite a naive 2*g overflow: g/(D-d) is approximately 1.
    let wide = ball_mill_critical_speed(
        Length::new(f64::MAX),
        Length::new(1.0),
        Acceleration::new(f64::MAX),
    )
    .unwrap();
    close(
        wide.value(),
        core::f64::consts::SQRT_2 / core::f64::consts::TAU,
        2e-15,
    );
}
