//! Source-qualified, fixed-volume thermal storage and its bounded inverse.
#![deny(unsafe_code)]

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimId, ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    NormalizedMaterialCardPack, NormalizedPack, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::thermal::transport::{
    ThermalTransportError,
    volumetric::{
        PreparedVolumetricThermalTransport as ThermalPath, VOLUMETRIC_HEAT_CAPACITY_DIMS as C,
        VolumetricThermalPathQuery as PathQuery,
    },
};
use fs_qty::{
    Dims, QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use std::{path::Path, process::Command};

const K: Dims = Dims([1, 1, -3, -1, 0, 0]);
const POLICY: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
const CAPACITY: [(f64, f64); 4] = [(275.0, 1e6), (300.0, 2e6), (340.0, 1.2e6), (350.0, 1.5e6)];
const CONDUCTIVITY: [(f64, f64); 4] = [(275.0, 0.2), (290.0, 0.4), (320.0, 0.1), (350.0, 0.7)];

fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn keys() -> [PropertyKey; 2] {
    [
        PropertyKey::new("volumetric-heat-capacity", C),
        PropertyKey::with_quantity(
            "thermal-conductivity",
            quantity(QuantityKind::ThermalConductivity),
        ),
    ]
}
fn at(p: QueryPoint, kelvin: f64) -> QueryPoint {
    p.with_quantity(
        "temperature",
        quantity(QuantityKind::AbsoluteTemperature),
        kelvin,
    )
    .unwrap()
}
fn point(kelvin: f64) -> QueryPoint {
    at(
        QueryPoint::new()
            .with_quantity("specimen", QuantitySpec::dimensional(Dims::NONE), 1.0)
            .unwrap()
            .with_quantity(
                "plate-difference",
                quantity(QuantityKind::TemperatureDifference),
                20.0,
            )
            .unwrap(),
        kelvin,
    )
}
fn state(name: &str) -> MaterialStateId {
    MaterialStateId {
        chemistry: name.into(),
        phase: "solid".into(),
        process: "explicit source context".into(),
        revision: 0,
    }
}
fn curve(index: usize, knots: &[(f64, f64)]) -> PropertyClaim {
    PropertyClaim {
        key: keys()[index].clone(),
        value: PropertyValue::Curve {
            abscissa: "temperature".into(),
            abscissa_dims: Temperature::DIMS,
            knots: knots.to_vec(),
            dims: [C, K][index],
        },
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                "temperature",
                quantity(QuantityKind::AbsoluteTemperature),
                knots[0].0,
                knots.last().unwrap().0,
            )
            .with_quantity(
                "plate-difference",
                quantity(QuantityKind::TemperatureDifference),
                20.0,
                20.0,
            )
            .with_quantity("specimen", QuantitySpec::dimensional(Dims::NONE), 1.0, 1.0),
        interpolation: InterpolationPolicy::LinearInside,
        uncertainty: UncertaintyModel::RelativeHalfWidth {
            fraction: 0.02,
            confidence: 0.95,
        },
        observations: Vec::new(),
        provenance: Provenance {
            source: "authored analytic fixture; not measured material".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}
fn card(claims: impl IntoIterator<Item = PropertyClaim>) -> MaterialCard {
    let mut set = ClaimSet::new();
    for claim in claims {
        set.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(state("analytic volumetric fixture"), set, Vec::new()).unwrap()
}
fn prepare(
    c: &MaterialCard,
    lo: &QueryPoint,
    hi: &QueryPoint,
    selections: [ClaimSelection; 2],
) -> Result<ThermalPath, ThermalTransportError> {
    let k = keys();
    ThermalPath::try_new(
        c,
        PathQuery {
            volumetric_heat_capacity: &k[0],
            conductivity: &k[1],
            temperature_axis: "temperature",
            lower: lo,
            upper: hi,
            selections,
        },
    )
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-12 * expected.abs().max(1e-10),
        "{actual:e} != {expected:e}"
    );
}
fn primitive(t: f64) -> f64 {
    if t <= 300.0 {
        let x = t - 275.0;
        1e6 * x + 20000.0 * x * x
    } else if t <= 340.0 {
        let x = t - 300.0;
        37.5e6 + 2e6 * x - 10000.0 * x * x
    } else {
        let x = t - 340.0;
        101.5e6 + 1.2e6 * x + 15000.0 * x * x
    }
}
fn interpolate(knots: &[(f64, f64)], t: f64) -> f64 {
    let pair = knots
        .windows(2)
        .find(|w| w[0].0 <= t && t <= w[1].0)
        .unwrap();
    pair[0].1 + (t - pair[0].0) * (pair[1].1 - pair[0].1) / (pair[1].0 - pair[0].0)
}

#[test]
fn g1_g3_piecewise_storage_inverse_and_transport_keep_the_selected_source_curves() {
    let c = curve(0, &CAPACITY);
    let k = curve(1, &CONDUCTIVITY);
    let pins = [
        ClaimSelection::Pinned(ClaimId(c.content_hash())),
        ClaimSelection::Pinned(ClaimId(k.content_hash())),
    ];
    let mut alternative = c.clone();
    alternative.provenance.source = "independent alternate fixture, identical numbers".into();
    let alternate_pin = ClaimSelection::Pinned(ClaimId(alternative.content_hash()));
    let card = card([c, k, alternative]);
    let (lo, hi) = (point(275.0), point(350.0));
    assert!(prepare(&card, &lo, &hi, [POLICY; 2]).is_err());
    let path = prepare(&card, &lo, &hi, pins).unwrap();
    let original = path.clone();
    assert_eq!(
        path.temperatures_k().collect::<Vec<_>>(),
        [275.0, 290.0, 300.0, 320.0, 340.0, 350.0]
    );
    assert_eq!(path.temperature_span_k(), [275.0, 350.0]);
    assert_eq!(path.temperature_axis(), "temperature");
    assert_eq!(path.query_endpoints(), &[lo.clone(), hi.clone()]);
    assert_eq!(path.material_card_identity(), card.content_hash());
    assert_eq!(
        path.identity(),
        prepare(&card, &lo, &hi, pins).unwrap().identity()
    );
    assert_ne!(
        path.identity(),
        prepare(&card, &lo, &hi, [alternate_pin, pins[1]])
            .unwrap()
            .identity()
    );
    for (i, source) in path.source_lines().iter().enumerate() {
        assert_eq!(
            source,
            &card
                .claims()
                .query_line_typed(&keys()[i], "temperature", &lo, &hi, pins[i])
                .unwrap()
        );
        for (_, answer) in &source.knots {
            assert_eq!(
                answer.evidence.value.uncertainty,
                UncertaintyModel::RelativeHalfWidth {
                    fraction: 0.02,
                    confidence: 0.95
                }
            );
            card.claims().verify_receipt(&answer.receipt).unwrap();
        }
    }
    let temperatures = [
        275.0, 282.5, 290.0, 300.0, 307.25, 320.0, 333.0, 340.0, 345.0, 350.0,
    ];
    for t in temperatures {
        let s = path.at_temperature(t).unwrap();
        let c = interpolate(&CAPACITY, t);
        let k = interpolate(&CONDUCTIVITY, t);
        assert_eq!(s.path_identity(), path.identity());
        assert_eq!(s.temperature_k(), t);
        close(s.volumetric_heat_capacity_j_m3_k(), c);
        close(s.conductivity_w_m_k(), k);
        close(s.diffusivity_m2_s(), k / c);
        close(s.effusivity_j_m2_k_sqrt_s(), (k * c).sqrt());
        for end in temperatures {
            let expected = primitive(end) - primitive(t);
            close(path.storage_change_j_m3(t, end).unwrap(), expected);
            close(
                path.temperature_after_storage_change(t, expected).unwrap(),
                end,
            );
        }
    }
    let mid = path.at_temperature(282.5).unwrap();
    let wrong_alpha = 0.5 * (0.2 / 1e6 + 0.4 / 1.6e6);
    assert!((mid.diffusivity_m2_s() - wrong_alpha).abs() > 1e-9);
    assert_eq!(path.storage_change_j_m3(275.0, 350.0).unwrap(), 115e6);
    assert_eq!(path.storage_change_j_m3(350.0, 275.0).unwrap(), -115e6);
    assert_eq!(
        path.temperature_after_storage_change(275.0, 115e6).unwrap(),
        350.0
    );
    assert_eq!(
        path.temperature_after_storage_change(350.0, -115e6)
            .unwrap(),
        275.0
    );
    assert_eq!(path, original);
}

#[test]
fn g0_g3_source_support_quantity_roles_and_unrepresentable_storage_refuse() {
    let (c, k) = (curve(0, &CAPACITY), curve(1, &CONDUCTIVITY));
    let good = card([c.clone(), k.clone()]);
    let (lo, hi) = (point(275.0), point(350.0));
    let path = prepare(&good, &lo, &hi, [POLICY; 2]).unwrap();
    let before = path.clone();
    let wrong = PropertyKey::new("specific-heat-capacity", Dims([2, 0, -2, -1, 0, 0]));
    assert!(
        ThermalPath::try_new(
            &good,
            PathQuery {
                volumetric_heat_capacity: &wrong,
                conductivity: &keys()[1],
                temperature_axis: "temperature",
                lower: &lo,
                upper: &hi,
                selections: [POLICY; 2]
            }
        )
        .is_err()
    );
    for bad in [
        point(275.0),
        point(0.0),
        hi.clone().with("specimen", 2.0).unwrap(),
        hi.clone().with("extra", 1.0).unwrap(),
    ] {
        assert!(prepare(&good, &lo, &bad, [POLICY; 2]).is_err());
    }
    assert!(
        QueryPoint::new()
            .with_quantity(
                "temperature",
                quantity(QuantityKind::AbsoluteTemperature),
                -1.0
            )
            .is_err()
    );
    let mut discrete = c.clone();
    discrete.interpolation = InterpolationPolicy::TabulatedOnly;
    let mut single = c.clone();
    single.validity = single.validity.with_quantity(
        "temperature",
        quantity(QuantityKind::AbsoluteTemperature),
        300.0,
        300.0,
    );
    let mut unbounded = c.clone();
    unbounded.value = PropertyValue::Scalar {
        value: 1e6,
        dims: C,
    };
    unbounded.interpolation = InterpolationPolicy::ConstantWithinValidity;
    unbounded.validity = ValidityDomain::unconstrained();
    for invalid in [
        discrete,
        single,
        unbounded,
        curve(0, &[(275.0, 1e6), (300.0, 0.0), (350.0, 1e6)]),
    ] {
        assert!(prepare(&card([invalid, k.clone()]), &lo, &hi, [POLICY; 2]).is_err());
    }
    let gap = card([
        curve(0, &[(275.0, 1e6), (300.0, 1e6)]),
        curve(0, &[(310.0, 1e6), (350.0, 1e6)]),
        k.clone(),
    ]);
    assert!(prepare(&gap, &lo, &hi, [POLICY; 2]).is_err());
    let mut second = c.clone();
    second.validity = second.validity.with_quantity(
        "plate-difference",
        quantity(QuantityKind::AbsoluteTemperature),
        20.0,
        20.0,
    );
    let other_temperature = |p: QueryPoint| {
        p.with_quantity(
            "plate-difference",
            quantity(QuantityKind::AbsoluteTemperature),
            20.0,
        )
        .unwrap()
    };
    assert!(
        prepare(
            &card([second, k.clone()]),
            &other_temperature(lo.clone()),
            &other_temperature(hi.clone()),
            [POLICY; 2]
        )
        .is_err()
    );
    let mut hidden = c.clone();
    if let PropertyValue::Curve { abscissa, .. } = &mut hidden.value {
        *abscissa = "bath-temperature".into();
    }
    let bath = |p: QueryPoint| {
        p.with_quantity(
            "bath-temperature",
            quantity(QuantityKind::AbsoluteTemperature),
            300.0,
        )
        .unwrap()
    };
    assert!(
        prepare(
            &card([hidden, k.clone()]),
            &bath(lo.clone()),
            &bath(hi.clone()),
            [POLICY; 2]
        )
        .is_err()
    );
    for t in [274.99, 350.01, f64::NAN, f64::INFINITY] {
        assert!(path.at_temperature(t).is_err());
        assert!(path.storage_change_j_m3(300.0, t).is_err());
        assert!(path.temperature_after_storage_change(t, 0.0).is_err());
    }
    for q in [
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
        116e6,
        -116e6,
        1e-300,
        -1e-300,
    ] {
        assert!(path.temperature_after_storage_change(300.0, q).is_err());
    }
    assert_eq!(path, before);
    let extreme = |value| {
        card([
            curve(0, &[(275.0, value), (350.0, value)]),
            curve(1, &[(275.0, value), (350.0, value)]),
        ])
    };
    let large = prepare(&extreme(1e308), &lo, &hi, [POLICY; 2]).unwrap();
    assert!(large.storage_change_j_m3(275.0, 350.0).is_err());
    let tiny = prepare(&extreme(1e-320), &lo, &hi, [POLICY; 2]).unwrap();
    assert!(
        tiny.storage_change_j_m3(300.0, 300.0_f64.next_up())
            .is_err()
    );
    // The requested heat is representable even though the full available interval is not.
    close(
        large
            .temperature_after_storage_change(300.0, 1e307)
            .unwrap(),
        300.1,
    );
}

fn compiled(id: &str) -> NormalizedMaterialCardPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let slug = format!("fsri-{id}-unconditioned-2022");
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let output = std::env::temp_dir().join(format!(
        "volumetric-curves-{}-{nonce}-{id}.fsmatpk",
        std::process::id()
    ));
    let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(root.join(format!("data/matdb/seed-v1/{slug}/manifest.tsv")))
        .arg("--out")
        .arg(&output)
        .env("CARGO_WORKSPACE_DIR", root)
        .current_dir(root)
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "{slug}: {}",
        String::from_utf8_lossy(&run.stderr)
    );
    let bytes = std::fs::read(output).unwrap();
    let pack = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(pack.content_hash(), &bytes).unwrap();
    assert_eq!(pack.claims().claims_ordered().count(), 8);
    let card = NormalizedMaterialCardPack::new(state(&slug), pack).unwrap();
    NormalizedMaterialCardPack::from_bytes_verified(card.content_hash(), &card.to_bytes()).unwrap()
}
fn source_point(card: &MaterialCard, t: f64) -> QueryPoint {
    let mut p = QueryPoint::new();
    for key in keys() {
        let (_, claim) = card
            .claims()
            .claims_ordered()
            .find(|(_, c)| c.key == key)
            .unwrap();
        for (name, &(lo, hi)) in claim.validity.bounds() {
            let value = if name == "temperature" {
                t
            } else {
                assert_eq!(lo, hi);
                lo
            };
            p = p
                .with_quantity(name, claim.validity.axis_quantities()[name], value)
                .unwrap();
        }
    }
    p
}

#[test]
fn g1_g3_four_real_foam_cards_drive_fixed_volume_heating_and_cooling() {
    // Independent source means, not extracted from the model under test.
    let sources = [
        (
            "pe-foam-pipe",
            [76344.0, 80182.33333333333, 84950.0, 91188.66666666667],
            [0.049, 0.058],
        ),
        (
            "rubber-foam-pipe",
            [
                68817.33333333333,
                70242.66666666667,
                73210.0,
                75710.33333333333,
            ],
            [0.034, 0.038],
        ),
        (
            "memory-foam-pad",
            [237307.0, 244425.66666666666, 248445.33333333334, 251435.0],
            [0.039, 0.042],
        ),
        (
            "xps-board",
            [35334.666666666664, 36374.666666666664, 37865.0, 38192.0],
            [0.028, 0.032],
        ),
    ];
    let mut warm = Vec::new();
    for (id, capacity, conductivity) in sources {
        let pack = compiled(id);
        let card = pack.card();
        let (lo, hi) = (source_point(card, 288.15), source_point(card, 313.15));
        let path = prepare(card, &lo, &hi, [POLICY; 2]).unwrap();
        assert_eq!(
            path.temperatures_k().collect::<Vec<_>>(),
            [288.15, 293.15, 303.15, 313.15]
        );
        for (i, source) in path.source_lines().iter().enumerate() {
            for (t, answer) in &source.knots {
                assert!(answer.receipt.observation_backed);
                assert_eq!(
                    answer.evidence.value.uncertainty,
                    UncertaintyModel::Unstated
                );
                card.claims().verify_receipt(&answer.receipt).unwrap();
                assert_eq!(
                    answer,
                    &card
                        .claims()
                        .query_typed(
                            &keys()[i],
                            &source_point(card, *t),
                            SelectionPolicy::SingleClaimOnly
                        )
                        .unwrap()
                );
            }
        }
        let s = path.at_temperature(303.15).unwrap();
        close(s.volumetric_heat_capacity_j_m3_k(), capacity[2]);
        close(
            s.conductivity_w_m_k(),
            (conductivity[0] + conductivity[1]) / 2.0,
        );
        close(
            s.diffusivity_m2_s(),
            (conductivity[0] + conductivity[1]) / (2.0 * capacity[2]),
        );
        let c15 = (capacity[0] + capacity[1]) / 2.0;
        let whole = 2.5 * (c15 + capacity[1])
            + 5.0 * (capacity[1] + capacity[2])
            + 5.0 * (capacity[2] + capacity[3]);
        close(path.storage_change_j_m3(288.15, 313.15).unwrap(), whole);
        close(path.storage_change_j_m3(313.15, 288.15).unwrap(), -whole);
        // Uniform, fixed gross volume2L, initially20°C, receives600J; no density is needed.
        let energy_per_volume = 600.0 / 0.002;
        let slope = (capacity[2] - capacity[1]) / 10.0;
        let delta = ((capacity[1] * capacity[1] + 2.0 * slope * energy_per_volume).sqrt()
            - capacity[1])
            / slope;
        let actual = path
            .temperature_after_storage_change(293.15, energy_per_volume)
            .unwrap();
        close(actual, 293.15 + delta);
        close(
            path.storage_change_j_m3(293.15, actual).unwrap() * 0.002,
            600.0,
        );
        close(
            path.temperature_after_storage_change(actual, -energy_per_volume)
                .unwrap(),
            293.15,
        );
        warm.push(actual);
        let density = PropertyKey::new("density", Dims([-3, 1, 0, 0, 0, 0]));
        assert!(
            card.claims()
                .query_typed(&density, &lo, SelectionPolicy::SingleClaimOnly)
                .is_err()
        );
        for axis in [
            "mean-test-temperature-as-local-model",
            "source-unconditioned",
        ] {
            let bad_lo = lo.clone().with(axis, 0.0).unwrap();
            let bad_hi = hi.clone().with(axis, 0.0).unwrap();
            assert!(prepare(card, &bad_lo, &bad_hi, [POLICY; 2]).is_err());
        }
        assert!(prepare(card, &lo, &at(hi.clone(), 318.15), [POLICY; 2]).is_err());
    }
    assert!(warm[2] < warm[0] && warm[0] < warm[1] && warm[1] < warm[3]);
}
