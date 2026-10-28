use super::*;

fn nitrogen_source() -> Nasa9StandardStateModelV1 {
    let text = include_str!(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/../../data/matdb/seed-v1/nitrogen-nasa9-mcbride2002/nasa9.tsv"
    ));
    let mut regions = Vec::new();
    for (name, lo, hi) in [("low", 200.0, 1000.0), ("high", 1000.0, 6000.0)] {
        let mut coefficients = [0.0; 9];
        let mut count = 0;
        for line in text.lines() {
            let fields: Vec<_> = line.split('\t').collect();
            if fields.first() == Some(&"coefficient") && fields[2] == name {
                let i: usize = fields[3].strip_prefix('a').unwrap().parse().unwrap();
                coefficients[i] = fields[4].parse().unwrap();
                count += 1;
            }
        }
        assert_eq!(count, 9);
        regions.push(Nasa9RegionV1::from_card(card(lo, hi, 100000.0, coefficients)).unwrap());
    }
    Nasa9StandardStateModelV1::new(
        SpeciesId::new("N2").unwrap(),
        MolarMass::new(0.0280134),
        convention(100000.0),
        regions,
    )
    .unwrap()
}

fn continuous(source: &Nasa9StandardStateModelV1, anchor: f64) -> Nasa9StandardStateModelV1 {
    source
        .with_caloric_continuity(
            Temperature::new(anchor),
            MolarEnergyQuantityV1::new(1e-3),
            MolarThermalQuantityV1::new(1e-5),
        )
        .unwrap()
}

#[test]
fn g0_caloric_reconciliation_preserves_source_cp_and_the_declared_anchor() {
    let raw = nitrogen_source();
    assert!(raw.caloric_reference().is_none());
    let corrected = continuous(&raw, 298.15);
    assert_eq!(raw.regions(), corrected.regions());
    let policy = corrected.caloric_reference().unwrap();
    assert_eq!(policy.anchor_region(), 0);
    assert_eq!(policy.anchor_temperature().value(), 298.15);
    assert_eq!(policy.regions()[0].enthalpy_offset().value(), 0.0);
    assert_eq!(policy.regions()[0].entropy_offset().value(), 0.0);
    assert!(policy.regions()[1].enthalpy_offset().value().abs() > 1e-7);
    assert!(policy.regions()[1].entropy_offset().value().abs() > 1e-9);
    for t in [200.0, 298.15, 900.0, 1000.0, 1500.0, 6000.0] {
        let original = raw.evaluate_properties(Temperature::new(t)).unwrap();
        let full = corrected.evaluate(Temperature::new(t)).unwrap();
        let p = full.properties();
        assert_eq!(p.cp().value().to_bits(), original.cp().value().to_bits());
        if t < 1000.0 {
            assert_eq!(p, original);
        }
        assert!((p.h().value() - original.h().value()).abs() <= 1e-3);
        assert!((p.s().value() - original.s().value()).abs() <= 1e-5);
        assert_eq!(
            p.u().value(),
            p.h().value() - UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K * t
        );
        assert_eq!(p.g().value(), p.h().value() - t * p.s().value());
        assert_eq!(
            p,
            corrected.evaluate_properties(Temperature::new(t)).unwrap()
        );
        assert_eq!(
            full.receipt().evaluator_version(),
            NASA9_CONTINUOUS_CALORIC_EVALUATOR_VERSION_V1
        );
        assert!(std::ptr::eq(
            policy,
            full.receipt().caloric_reference().unwrap()
        ));
        assert_eq!(
            full.receipt().card_identity(),
            raw.regions()[full.receipt().region_index()].card_identity()
        );
    }
    let left = raw.evaluate_source_region(0, 1000.0).unwrap();
    let right = corrected
        .evaluate_properties(Temperature::new(1000.0))
        .unwrap();
    assert_eq!(left.h().value().to_bits(), right.h().value().to_bits());
    assert_eq!(left.s().value().to_bits(), right.s().value().to_bits());
    assert_eq!(corrected, continuous(&corrected, 298.15));
    let high_anchor = continuous(&raw, 2000.0);
    assert_eq!(high_anchor.caloric_reference().unwrap().anchor_region(), 1);
    assert_eq!(
        high_anchor
            .evaluate_properties(Temperature::new(2000.0))
            .unwrap(),
        raw.evaluate_properties(Temperature::new(2000.0)).unwrap()
    );
    assert_ne!(
        high_anchor
            .evaluate_properties(Temperature::new(300.0))
            .unwrap()
            .h(),
        raw.evaluate_properties(Temperature::new(300.0))
            .unwrap()
            .h()
    );
}

#[test]
fn g1_reconciled_caloric_differences_and_gibbs_share_the_source_derivatives() {
    let model = continuous(&nitrogen_source(), 298.15);
    for t in [900.0, 1100.0, 2500.0] {
        let epsilon = 1e-3;
        let lo = Temperature::new(t - epsilon);
        let hi = Temperature::new(t + epsilon);
        let center = model.evaluate_properties(Temperature::new(t)).unwrap();
        let dh = model.molar_enthalpy_change(lo, hi).unwrap().value() / (2.0 * epsilon);
        let ds = model.molar_entropy_change(lo, hi).unwrap().value() / (2.0 * epsilon);
        assert!((dh - center.cp().value()).abs() < 1e-6);
        assert!((ds - center.cp().value() / t).abs() < 1e-9);
        let dg = (model.evaluate_properties(hi).unwrap().g().value()
            - model.evaluate_properties(lo).unwrap().g().value())
            / (2.0 * epsilon);
        assert!((dg + center.s().value()).abs() < 1e-5);
    }
    let lo = Temperature::new(300.0);
    let hi = Temperature::new(2500.0);
    let p0 = model.evaluate_properties(lo).unwrap();
    let p1 = model.evaluate_properties(hi).unwrap();
    let dh = model.molar_enthalpy_change(lo, hi).unwrap().value();
    let ds = model.molar_entropy_change(lo, hi).unwrap().value();
    assert!((dh - (p1.h().value() - p0.h().value())).abs() < 1e-7);
    assert!((ds - (p1.s().value() - p0.s().value())).abs() < 1e-10);
    assert_eq!(model.molar_enthalpy_change(hi, lo).unwrap().value(), -dh);
    assert_eq!(model.molar_entropy_change(hi, lo).unwrap().value(), -ds);
    let mixture = crate::mixture::FrozenIdealGasMixtureModelV1::new(
        vec![model.clone()],
        Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
    )
    .unwrap();
    let mixed = mixture.evaluate(hi, Pressure::new(100000.0)).unwrap();
    assert_eq!(mixed.properties().h(), p1.h());
    assert!(
        (mixture.specific_enthalpy_change(lo, hi).unwrap().value() - dh / 0.0280134).abs() < 1e-7
    );
    assert!(
        (mixture
            .specific_entropy_change(lo, Pressure::new(100000.0), hi, Pressure::new(100000.0))
            .unwrap()
            .value()
            - ds / 0.0280134)
            .abs()
            < 1e-10
    );
}

#[test]
fn g1_reconciled_relative_heat_survives_adjacent_temperatures_across_the_seam() {
    let model = continuous(&nitrogen_source(), 298.15);
    let seam = 1000.0_f64;
    let lo = seam.next_down();
    let hi = seam.next_up();
    let left_cp = model.evaluate_source_region(0, seam).unwrap().cp().value();
    let right_cp = model.evaluate_source_region(1, seam).unwrap().cp().value();
    let expected = left_cp * (seam - lo) + right_cp * (hi - seam);
    let dh = model
        .molar_enthalpy_change(Temperature::new(lo), Temperature::new(hi))
        .unwrap()
        .value();
    let ds = model
        .molar_entropy_change(Temperature::new(lo), Temperature::new(hi))
        .unwrap()
        .value();
    assert!(dh > 0.0 && ds > 0.0);
    assert!((dh / expected - 1.0).abs() < 1e-12);
    assert!((ds / (expected / seam) - 1.0).abs() < 1e-12);
}

#[test]
fn g0_caloric_reconciliation_refuses_gaps_and_excess_cumulative_adjustment() {
    let raw = nitrogen_source();
    for (h, s) in [
        (0.0, 1e-5),
        (1e-3, 0.0),
        (-1.0, 1.0),
        (f64::INFINITY, 1.0),
        (1.0, f64::NAN),
    ] {
        assert!(
            raw.with_caloric_continuity(
                Temperature::new(298.15),
                MolarEnergyQuantityV1::new(h),
                MolarThermalQuantityV1::new(s)
            )
            .is_err()
        );
    }
    for t in [0.0, 199.0, 6001.0, f64::NAN] {
        assert!(
            raw.with_caloric_continuity(
                Temperature::new(t),
                MolarEnergyQuantityV1::new(1.0),
                MolarThermalQuantityV1::new(1.0)
            )
            .is_err()
        );
    }
    let make = |gap: f64| {
        let mut regions = Vec::new();
        for (index, (lo, hi)) in [(200.0, 1000.0), (1000.0 + gap, 2000.0), (2000.0, 3000.0)]
            .into_iter()
            .enumerate()
        {
            let mut a = [0.0; 9];
            a[2] = 3.5;
            a[7] = index as f64;
            regions.push(Nasa9RegionV1::from_card(card(lo, hi, 100000.0, a)).unwrap());
        }
        Nasa9StandardStateModelV1::new(
            SpeciesId::new("N2").unwrap(),
            MolarMass::new(0.028),
            convention(100000.0),
            regions,
        )
        .unwrap()
    };
    assert!(matches!(
        make(0.0).with_caloric_continuity(
            Temperature::new(300.0),
            MolarEnergyQuantityV1::new(10.0),
            MolarThermalQuantityV1::new(1.0)
        ),
        Err(ThermochemErrorV1::CaloricReferenceRefused { region: 2, .. })
    ));
    assert!(
        make(1.0)
            .with_caloric_continuity(
                Temperature::new(300.0),
                MolarEnergyQuantityV1::new(100.0),
                MolarThermalQuantityV1::new(1.0)
            )
            .is_err()
    );
    assert!(raw.caloric_reference().is_none());
}
