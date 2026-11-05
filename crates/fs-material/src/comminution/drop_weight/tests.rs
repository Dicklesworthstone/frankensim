use super::*;
use fs_blake3::hash_bytes;
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, Provenance, UncertaintyModel};
use fs_qty::{Dims, QuantitySpec, Temperature};

fn point() -> QueryPoint {
    QueryPoint::new()
        .with_quantity("sample", QuantitySpec::dimensional(Dims::NONE), 1.0)
        .unwrap()
}
fn domain(low: f64, high: f64) -> ImpactEnergyDomain {
    ImpactEnergyDomain::new([
        SpecificGrindingEnergy::new(low),
        SpecificGrindingEnergy::new(high),
    ])
    .unwrap()
}
fn fixture(
    a: f64,
    b: f64,
    edit: impl FnOnce(&mut ClaimSet, &mut [PropertyClaim; 2]),
) -> NormalizedPack {
    let provenance = Provenance {
        source: "authored mathematical test fixture; not measured ore".into(),
        license: "project test fixture".into(),
        artifact: Some(hash_bytes(b"synthetic JKDWT row")),
    };
    let mut set = ClaimSet::new();
    let observation = set
        .register_observation(ObservationDataset {
            specimen: "authored sample".into(),
            method: "synthetic JKDWT fit fixture".into(),
            artifact: hash_bytes(b"synthetic JKDWT row"),
            caveats: "no plant qualification".into(),
            provenance: provenance.clone(),
        })
        .unwrap();
    let [a_key, b_key] = parameter_keys();
    let mut claims = [(a_key, a), (b_key, b)].map(|(key, value)| PropertyClaim {
        value: PropertyValue::Scalar {
            value,
            dims: key.dims(),
        },
        key,
        validity: ValidityDomain::unconstrained().with_quantity(
            "sample",
            QuantitySpec::dimensional(Dims::NONE),
            1.0,
            1.0,
        ),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![observation],
        provenance: provenance.clone(),
    });
    edit(&mut set, &mut claims);
    for claim in claims {
        set.insert_claim(claim).unwrap();
    }
    NormalizedPack::new(
        "synthetic-jkdwt",
        "unit-test",
        hash_bytes(b"synthetic JKDWT row"),
        "project fixture",
        set,
        vec![],
        vec![],
    )
    .unwrap()
}
fn resolve(
    pack: &NormalizedPack,
    domain: ImpactEnergyDomain,
) -> Result<DropWeightCard, DropWeightError> {
    resolve_jk_drop_weight(pack, &point(), SelectionPolicy::SingleClaimOnly, domain)
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= expected.abs().max(f64::MIN_POSITIVE) * 8e-14,
        "{actual} != {expected}"
    );
}

#[test]
fn fraction_and_si_units_preserve_reported_curve_and_original_evidence() {
    // Authored equivalent of A=60 percent and b=0.5/(kWh/metric tonne).
    let pack = fixture(0.6, 0.5 / 3600.0, |_, _| {});
    let card = resolve(&pack, domain(0.0, 72_000.0)).unwrap();
    let energy = SpecificGrindingEnergy::new(7200.0);
    let response = card.response(energy).unwrap();
    close(
        response.passing_one_tenth().value(),
        0.6 * (1.0 - (-1.0_f64).exp()),
    );
    close(
        response.retained_one_tenth().value(),
        0.4 + 0.6 * (-1.0_f64).exp(),
    );
    close(
        card.energy_for_t10(response.passing_one_tenth())
            .unwrap()
            .value(),
        energy.value(),
    );
    let source = card.calibration_source();
    assert_eq!(source.pack_hash(), pack.content_hash());
    assert_eq!(source.pack_id(), pack.pack_id());
    assert_eq!(source.query(), &point());
    assert!(
        source
            .observation()
            .caveats
            .contains("no plant qualification")
    );
    for (claim, answer) in source.claims().iter().zip(source.answers()) {
        assert_eq!(pack.claims().claim(answer.receipt.selected), Some(claim));
        pack.claims().verify_receipt(&answer.receipt).unwrap();
        assert_eq!(
            answer.evidence.value.uncertainty,
            UncertaintyModel::Unstated
        );
    }
    assert_eq!(card.clone(), card);
}

#[test]
fn response_handles_zero_small_energy_and_near_asymptote_without_loss_of_tails() {
    let card = resolve(&fixture(1.0, 1.0, |_, _| {}), domain(0.0, 1000.0)).unwrap();
    let zero = card.response(SpecificGrindingEnergy::new(0.0)).unwrap();
    assert_eq!(zero.passing_one_tenth().value(), 0.0);
    assert_eq!(zero.retained_one_tenth().value(), 1.0);
    assert_eq!(
        card.energy_for_t10(Dimensionless::new(0.0))
            .unwrap()
            .value(),
        0.0
    );
    for e in [1e-300, 1e-20, 1e-8, 0.01, 1.0] {
        let t = card.response(SpecificGrindingEnergy::new(e)).unwrap();
        close(t.passing_one_tenth().value(), -(-e).exp_m1());
        close(
            card.energy_for_t10(t.passing_one_tenth()).unwrap().value(),
            e,
        );
    }
    // Near saturation the inverse is ill-conditioned in t10: test its
    // forward residual instead of expecting lost input bits to be recovered.
    let near = card.response(SpecificGrindingEnergy::new(12.0)).unwrap();
    let inverse = card.energy_for_t10(near.passing_one_tenth()).unwrap();
    close(
        card.response(inverse).unwrap().passing_one_tenth().value(),
        near.passing_one_tenth().value(),
    );
    let saturated = card.response(SpecificGrindingEnergy::new(40.0)).unwrap();
    assert_eq!(saturated.passing_one_tenth().value(), 1.0);
    assert!(saturated.retained_one_tenth().value() > 0.0);
    close(saturated.retained_one_tenth().value(), (-40.0_f64).exp());
    assert!(card.energy_for_t10(saturated.passing_one_tenth()).is_err());
    let almost_one = f64::from_bits(1.0_f64.to_bits() - 1);
    close(
        card.energy_for_t10(Dimensionless::new(almost_one))
            .unwrap()
            .value(),
        -(1.0 - almost_one).ln(),
    );
}

#[test]
fn same_a_times_b_does_not_identify_the_high_energy_response() {
    let first = resolve(&fixture(0.5, 2.0, |_, _| {}), domain(0.0, 10.0)).unwrap();
    let second = resolve(&fixture(1.0, 1.0, |_, _| {}), domain(0.0, 10.0)).unwrap();
    assert_eq!(
        first.asymptotic_t10().value() * first.energy_sensitivity().value(),
        second.asymptotic_t10().value() * second.energy_sensitivity().value()
    );
    assert!(
        second
            .response(SpecificGrindingEnergy::new(2.0))
            .unwrap()
            .passing_one_tenth()
            .value()
            > first
                .response(SpecificGrindingEnergy::new(2.0))
                .unwrap()
                .passing_one_tenth()
                .value()
                + 0.3
    );
}

#[test]
fn domain_and_nonphysical_parameters_refuse_without_extrapolation() {
    for (a, b) in [
        (0.0, 1.0),
        (-0.1, 1.0),
        (60.0, 1.0),
        (0.6, 0.0),
        (0.6, -1.0),
    ] {
        assert!(resolve(&fixture(a, b, |_, _| {}), domain(0.0, 100.0)).is_err());
    }
    for bounds in [
        [-1.0, 1.0],
        [2.0, 1.0],
        [0.0, f64::INFINITY],
        [f64::NAN, 1.0],
    ] {
        assert!(ImpactEnergyDomain::new(bounds.map(SpecificGrindingEnergy::new)).is_err());
    }
    let card = resolve(&fixture(0.6, 0.001, |_, _| {}), domain(360.0, 9000.0)).unwrap();
    for e in [-1.0, 0.0, 359.0, 9001.0, f64::NAN, f64::INFINITY] {
        assert!(card.response(SpecificGrindingEnergy::new(e)).is_err());
    }
    for t in [-0.1, 0.0, 0.6, 0.61, f64::NAN, f64::INFINITY] {
        assert!(card.energy_for_t10(Dimensionless::new(t)).is_err());
    }
}

#[test]
fn mixed_observations_support_and_unsupported_conditions_refuse() {
    let mixed = fixture(0.6, 1.0, |set, claims| {
        let first = set.observation(claims[0].observations[0]).unwrap();
        let mut other = first.clone();
        other.specimen = "different specimen".into();
        claims[1].observations = vec![set.register_observation(other).unwrap()];
    });
    assert!(resolve(&mixed, domain(0.0, 10.0)).is_err());
    let support = fixture(0.6, 1.0, |_, claims| {
        claims[1].validity = ValidityDomain::unconstrained().with_quantity(
            "sample",
            QuantitySpec::dimensional(Dims::NONE),
            1.0,
            2.0,
        );
    });
    assert!(resolve(&support, domain(0.0, 10.0)).is_err());
    let pack = fixture(0.6, 1.0, |_, _| {});
    let extra = point()
        .with_quantity(
            "temperature",
            QuantitySpec::dimensional(Temperature::DIMS),
            298.15,
        )
        .unwrap();
    assert!(
        resolve_jk_drop_weight(
            &pack,
            &extra,
            SelectionPolicy::SingleClaimOnly,
            domain(0.0, 10.0)
        )
        .is_err()
    );
    assert!(
        resolve_jk_drop_weight(
            &pack,
            &QueryPoint::new(),
            SelectionPolicy::SingleClaimOnly,
            domain(0.0, 10.0)
        )
        .is_err()
    );
}
