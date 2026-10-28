//! G0/G1/G3: prepared transport preserves source models and refuses unsupported paths.
//! Authored numerical fixtures are not experimental material validation.
#![deny(unsafe_code)]

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    NormalizedMaterialCardPack, NormalizedPack, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::thermal::transport::{
    PreparedThermalTransport, ThermalTransportError, ThermalTransportSourceQuery,
};
use fs_qty::{
    Dims, QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use std::{path::Path, process::Command};

const RHO: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const CP: Dims = Dims([2, 0, -2, -1, 0, 0]);
const K: Dims = Dims([1, 1, -3, -1, 0, 0]);
const INVERSE_TIME: Dims = Dims([0, 0, -1, 0, 0, 0]);
const POLICY: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);

fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

fn keys() -> [PropertyKey; 3] {
    [
        PropertyKey::new("density", RHO),
        PropertyKey::new("specific-heat-capacity", CP),
        PropertyKey::with_quantity(
            "thermal-conductivity",
            quantity(QuantityKind::ThermalConductivity),
        ),
    ]
}

fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "temperature",
            quantity(QuantityKind::AbsoluteTemperature),
            t,
        )
        .unwrap()
}

fn domain(low: f64, high: f64) -> ValidityDomain {
    ValidityDomain::unconstrained().with_quantity(
        "temperature",
        quantity(QuantityKind::AbsoluteTemperature),
        low,
        high,
    )
}

fn curve(index: usize, knots: Vec<(f64, f64)>) -> PropertyClaim {
    PropertyClaim {
        key: keys()[index].clone(),
        validity: domain(knots[0].0, knots[knots.len() - 1].0),
        value: PropertyValue::Curve {
            abscissa: "temperature".into(),
            abscissa_dims: Temperature::DIMS,
            knots,
            dims: [RHO, CP, K][index],
        },
        uncertainty: UncertaintyModel::RelativeHalfWidth {
            fraction: 0.02,
            confidence: 0.95,
        },
        interpolation: InterpolationPolicy::LinearInside,
        observations: Vec::new(),
        provenance: Provenance {
            source: "authored transport regression fixture; not measured material data".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    }
}

fn source_claims() -> [PropertyClaim; 3] {
    [
        curve(0, vec![(280.0, 1000.0), (300.0, 1200.0), (360.0, 1600.0)]),
        curve(1, vec![(280.0, 400.0), (320.0, 900.0), (360.0, 600.0)]),
        curve(2, vec![(280.0, 1.0), (340.0, 4.0), (360.0, 10.0)]),
    ]
}

fn card(claims: impl IntoIterator<Item = PropertyClaim>) -> MaterialCard {
    let mut set = ClaimSet::new();
    for claim in claims {
        set.insert_claim(claim).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "authored transport regression".into(),
            phase: "solid".into(),
            process: "numerical source models".into(),
            revision: 0,
        },
        set,
        Vec::new(),
    )
    .unwrap()
}

fn prepare(
    card: &MaterialCard,
    keys: &[PropertyKey; 3],
    lower: &QueryPoint,
    upper: &QueryPoint,
) -> Result<PreparedThermalTransport, ThermalTransportError> {
    PreparedThermalTransport::try_new(
        card,
        ThermalTransportSourceQuery {
            density: &keys[0],
            specific_heat: &keys[1],
            conductivity: &keys[2],
            temperature_axis: "temperature",
            lower,
            upper,
            selections: [POLICY; 3],
        },
    )
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-13 * expected.abs().max(1e-300),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_union_grid_preserves_each_input_before_deriving_transport() {
    let card = card(source_claims());
    let keys = keys();
    let prepared = prepare(&card, &keys, &point(280.0), &point(360.0)).unwrap();
    assert_eq!(
        prepared.temperatures_k().collect::<Vec<_>>(),
        [280.0, 300.0, 320.0, 340.0, 360.0]
    );
    for t in [
        280.0, 285.0, 295.0, 300.0, 310.0, 330.0, 350.0, 359.0, 360.0,
    ] {
        let expected = keys.each_ref().map(|key| {
            card.claims()
                .query_typed(key, &point(t), SelectionPolicy::SingleClaimOnly)
                .unwrap()
                .evidence
                .value
                .value
        });
        let result = prepared.at_temperature(t).unwrap();
        close(result.density_kg_m3(), expected[0]);
        close(result.specific_heat_j_kg_k(), expected[1]);
        close(result.conductivity_w_m_k(), expected[2]);
        close(
            result.volumetric_heat_capacity_j_m3_k(),
            expected[0] * expected[1],
        );
        close(
            result.diffusivity_m2_s(),
            expected[2] / (expected[0] * expected[1]),
        );
        close(
            result.effusivity_j_m2_k_sqrt_s(),
            (expected[2] * expected[0] * expected[1]).sqrt(),
        );
        assert_eq!(result.temperature_k(), t);
        assert_eq!(result.path_identity(), prepared.identity());
    }
    // On the union segment 300..320 K, all three inputs are affine but their
    // products, ratios and square roots are not. Interpolating derived endpoint
    // values would produce observably different answers at 310 K.
    let low = prepared.at_temperature(300.0).unwrap();
    let high = prepared.at_temperature(320.0).unwrap();
    let mid = prepared.at_temperature(310.0).unwrap();
    close(
        mid.volumetric_heat_capacity_j_m3_k(),
        (3800.0 / 3.0) * 775.0,
    );
    assert!(
        (mid.volumetric_heat_capacity_j_m3_k()
            - 0.5
                * (low.volumetric_heat_capacity_j_m3_k() + high.volumetric_heat_capacity_j_m3_k()))
        .abs()
            > 1000.0
    );
    assert!(
        (mid.diffusivity_m2_s() - 0.5 * (low.diffusivity_m2_s() + high.diffusivity_m2_s())).abs()
            > 1e-9
    );
    assert!(
        (mid.effusivity_j_m2_k_sqrt_s()
            - 0.5 * (low.effusivity_j_m2_k_sqrt_s() + high.effusivity_j_m2_k_sqrt_s()))
        .abs()
            > 1.0
    );
}

#[test]
fn g3_original_lines_uncertainties_receipts_and_selection_identity_survive() {
    let card = card(source_claims());
    let keys = keys();
    let lower = point(280.0);
    let upper = point(360.0);
    let prepared = prepare(&card, &keys, &lower, &upper).unwrap();
    assert_eq!(prepared.material_card_identity(), card.content_hash());
    assert_eq!(prepared.query_endpoints(), [&lower, &upper]);
    assert_eq!(prepared.temperature_span_k(), [280.0, 360.0]);
    assert_eq!(prepared.clone(), prepared);
    for (key, retained) in keys.iter().zip(prepared.source_lines()) {
        let original = card
            .claims()
            .query_line_typed(key, "temperature", &lower, &upper, POLICY)
            .unwrap();
        assert_eq!(*retained, original);
        for (_, answer) in &retained.knots {
            card.claims().verify_receipt(&answer.receipt).unwrap();
            assert_eq!(
                answer.evidence.value.uncertainty,
                UncertaintyModel::RelativeHalfWidth {
                    fraction: 0.02,
                    confidence: 0.95,
                }
            );
        }
    }
    let again = prepare(&card, &keys, &lower, &upper).unwrap();
    assert_eq!(again.identity(), prepared.identity());
    let selected = prepared
        .source_lines()
        .each_ref()
        .map(|line| line.knots[0].1.receipt.selected);
    let pinned = PreparedThermalTransport::try_new(
        &card,
        ThermalTransportSourceQuery {
            density: &keys[0],
            specific_heat: &keys[1],
            conductivity: &keys[2],
            temperature_axis: "temperature",
            lower: &lower,
            upper: &upper,
            selections: selected.map(ClaimSelection::Pinned),
        },
    )
    .unwrap();
    assert_ne!(pinned.identity(), prepared.identity());
    for (line, expected) in pinned.source_lines().iter().zip(selected) {
        assert_eq!(line.knots[0].1.receipt.selected, expected);
    }
    close(
        pinned.at_temperature(310.0).unwrap().diffusivity_m2_s(),
        prepared.at_temperature(310.0).unwrap().diffusivity_m2_s(),
    );
    let changed_lower = lower
        .with_quantity(
            "fixed-context-tag",
            QuantitySpec::dimensional(Dims::NONE),
            1.0,
        )
        .unwrap();
    let changed_upper = upper
        .with_quantity(
            "fixed-context-tag",
            QuantitySpec::dimensional(Dims::NONE),
            1.0,
        )
        .unwrap();
    let changed = prepare(&card, &keys, &changed_lower, &changed_upper).unwrap();
    assert_ne!(changed.identity(), prepared.identity());
    assert_eq!(
        changed.material_card_identity(),
        prepared.material_card_identity()
    );
}

#[test]
fn g0_unsupported_discrete_gapped_or_nonpositive_sources_refuse() {
    let keys = keys();
    let lower = point(280.0);
    let upper = point(360.0);
    let supported = card(source_claims());
    assert!(matches!(
        prepare(&supported, &keys, &point(279.0), &upper),
        Err(ThermalTransportError::Source { .. })
    ));
    assert!(prepare(&supported, &keys, &upper, &lower).is_err());
    assert!(prepare(&supported, &keys, &lower, &lower).is_err());
    let changed_lower = lower.clone().with("condition", 0.0).unwrap();
    let changed_upper = upper.clone().with("condition", 1.0).unwrap();
    assert!(prepare(&supported, &keys, &changed_lower, &changed_upper).is_err());

    let mut claims = source_claims();
    claims[2].interpolation = InterpolationPolicy::TabulatedOnly;
    let discrete = card(claims);
    assert!(
        matches!(prepare(&discrete, &keys, &lower, &upper), Err(ThermalTransportError::Source { property, .. }) if property == "thermal-conductivity")
    );

    let mut claims = source_claims();
    claims[2].value = PropertyValue::Scalar {
        value: 2.0,
        dims: K,
    };
    claims[2].validity = ValidityDomain::unconstrained();
    claims[2].interpolation = InterpolationPolicy::ConstantWithinValidity;
    let unknown = card(claims);
    assert!(matches!(
        prepare(&unknown, &keys, &lower, &upper),
        Err(ThermalTransportError::Invalid(_))
    ));

    for value in [0.0, -1.0] {
        let mut claims = source_claims();
        claims[2] = curve(2, vec![(280.0, 1.0), (340.0, value), (360.0, 10.0)]);
        assert!(prepare(&card(claims), &keys, &lower, &upper).is_err());
    }

    // Both endpoints lie in valid triangles of the SAME source. The absent
    // 300..340 K cells must not be bridged by endpoint interpolation.
    let mut claims = source_claims();
    claims[2].value = PropertyValue::Surface2D {
        axes: ["temperature".into(), "fixed-moisture".into()],
        axis_dims: [Temperature::DIMS, Dims::NONE],
        samples: vec![
            ([280.0, 0.0], 1.0),
            ([300.0, 0.0], 2.0),
            ([280.0, 1.0], 1.5),
            ([300.0, 1.0], 2.5),
            ([340.0, 0.0], 4.0),
            ([360.0, 0.0], 5.0),
            ([340.0, 1.0], 4.5),
            ([360.0, 1.0], 5.5),
        ],
        triangles: vec![[0, 1, 2], [1, 3, 2], [4, 5, 6], [5, 7, 6]],
        dims: K,
    };
    claims[2].validity = domain(280.0, 360.0).with_quantity(
        "fixed-moisture",
        QuantitySpec::dimensional(Dims::NONE),
        0.0,
        1.0,
    );
    let gapped = card(claims);
    let lower = lower
        .with_quantity("fixed-moisture", QuantitySpec::dimensional(Dims::NONE), 0.5)
        .unwrap();
    let upper = upper
        .with_quantity("fixed-moisture", QuantitySpec::dimensional(Dims::NONE), 0.5)
        .unwrap();
    for point in [&lower, &upper] {
        assert!(
            gapped
                .claims()
                .query_typed(&keys[2], point, SelectionPolicy::SingleClaimOnly)
                .is_ok()
        );
    }
    assert!(
        matches!(prepare(&gapped, &keys, &lower, &upper), Err(ThermalTransportError::Source { property, .. }) if property == "thermal-conductivity")
    );
}

#[test]
fn g0_temperature_roles_and_hidden_frequency_dependencies_cannot_be_relabelled() {
    let keys = keys();
    let supported = card(source_claims());
    let lower = point(280.0);
    let upper = point(360.0);
    let mut wrong_keys = keys.clone();
    wrong_keys[0] = keys[1].clone();
    assert!(prepare(&supported, &wrong_keys, &lower, &upper).is_err());
    let difference = |t| {
        QueryPoint::new()
            .with_quantity(
                "temperature",
                quantity(QuantityKind::TemperatureDifference),
                t,
            )
            .unwrap()
    };
    assert!(prepare(&supported, &keys, &difference(280.0), &difference(360.0)).is_err());

    // Reproduces the formerly admitted pressure curve presented as kelvin:
    // legacy source/query axes are untyped, while the value axis declares Pa.
    let pressure_dims = Dims([-1, 1, -2, 0, 0, 0]);
    let mut claims = source_claims();
    for claim in &mut claims {
        if let PropertyValue::Curve { abscissa_dims, .. } = &mut claim.value {
            *abscissa_dims = pressure_dims;
        }
        claim.validity = ValidityDomain::unconstrained().with("temperature", 280.0, 360.0);
    }
    let pressure = card(claims);
    let untyped = |t| QueryPoint::new().with("temperature", t).unwrap();
    assert!(
        pressure
            .claims()
            .query_line_typed(
                &keys[0],
                "temperature",
                &untyped(280.0),
                &untyped(360.0),
                POLICY
            )
            .is_ok()
    );
    assert!(matches!(
        prepare(&pressure, &keys, &untyped(280.0), &untyped(360.0)),
        Err(ThermalTransportError::Invalid(_))
    ));

    // A differently named legacy frequency curve is supported by the general
    // source evaluator, but a nonzero-frequency cut is not static transport.
    let mut claims = source_claims();
    claims[2].value = PropertyValue::Curve {
        abscissa: "f-test".into(),
        abscissa_dims: INVERSE_TIME,
        knots: vec![(0.0, 1.0), (10000.0, 10.0)],
        dims: K,
    };
    claims[2].validity = domain(280.0, 360.0).with("f-test", 0.0, 10000.0);
    let dynamic = card(claims.clone());
    let low_f = lower.clone().with("f-test", 1000.0).unwrap();
    let high_f = upper.clone().with("f-test", 1000.0).unwrap();
    assert!(
        dynamic
            .claims()
            .query_line_typed(&keys[2], "temperature", &low_f, &high_f, POLICY)
            .is_ok()
    );
    assert!(matches!(
        prepare(&dynamic, &keys, &low_f, &high_f),
        Err(ThermalTransportError::Invalid(_))
    ));
    claims[2].validity = domain(280.0, 360.0).with_quantity(
        "f-test",
        QuantitySpec::dimensional(INVERSE_TIME),
        0.0,
        10000.0,
    );
    let typed_frequency = card(claims);
    let at_frequency = |t, f| {
        point(t)
            .with_quantity("f-test", QuantitySpec::dimensional(INVERSE_TIME), f)
            .unwrap()
    };
    assert!(
        prepare(
            &typed_frequency,
            &keys,
            &at_frequency(280.0, 1000.0),
            &at_frequency(360.0, 1000.0)
        )
        .is_err()
    );
    let static_path = prepare(
        &typed_frequency,
        &keys,
        &at_frequency(280.0, 0.0),
        &at_frequency(360.0, 0.0),
    )
    .unwrap();
    assert_eq!(
        static_path
            .at_temperature(310.0)
            .unwrap()
            .conductivity_w_m_k(),
        1.0
    );

    let mut claims = source_claims();
    claims[2].value = PropertyValue::Curve {
        abscissa: "second-temperature".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(400.0, 1.0), (600.0, 10.0)],
        dims: K,
    };
    claims[2].validity = domain(280.0, 360.0).with_quantity(
        "second-temperature",
        quantity(QuantityKind::AbsoluteTemperature),
        400.0,
        600.0,
    );
    let second_temperature = card(claims);
    let low = lower
        .with_quantity(
            "second-temperature",
            quantity(QuantityKind::AbsoluteTemperature),
            500.0,
        )
        .unwrap();
    let high = upper
        .with_quantity(
            "second-temperature",
            quantity(QuantityKind::AbsoluteTemperature),
            500.0,
        )
        .unwrap();
    assert!(prepare(&second_temperature, &keys, &low, &high).is_err());
}

#[test]
fn g0_extreme_derived_values_fail_safely_without_poisoning_the_prepared_path() {
    let keys = keys();
    let extreme = card([
        curve(0, vec![(280.0, 1e300), (360.0, 1e-300)]),
        curve(1, vec![(280.0, 1e-300), (360.0, 1e300)]),
        curve(2, vec![(280.0, 1.0), (360.0, 1.0)]),
    ]);
    let prepared = prepare(&extreme, &keys, &point(280.0), &point(360.0)).unwrap();
    let before = prepared.clone();
    // Endpoint products are representable; the interior product is not.
    close(
        prepared
            .at_temperature(280.0)
            .unwrap()
            .volumetric_heat_capacity_j_m3_k(),
        1.0,
    );
    assert!(matches!(
        prepared.at_temperature(320.0),
        Err(ThermalTransportError::Invalid(_))
    ));
    close(
        prepared
            .at_temperature(360.0)
            .unwrap()
            .volumetric_heat_capacity_j_m3_k(),
        1.0,
    );
    for invalid in [279.0, 361.0, 0.0, f64::NAN, f64::INFINITY] {
        assert!(prepared.at_temperature(invalid).is_err());
    }
    assert_eq!(prepared, before);

    // k*(rho*c) overflows, although its square root is representable. The
    // production factorized calculation must still deliver finite effusivity.
    let large = card([
        curve(0, vec![(280.0, 1e200), (360.0, 1e200)]),
        curve(1, vec![(280.0, 1.0), (360.0, 1.0)]),
        curve(2, vec![(280.0, 1e200), (360.0, 1e200)]),
    ]);
    let large = prepare(&large, &keys, &point(280.0), &point(360.0))
        .unwrap()
        .at_temperature(320.0)
        .unwrap();
    close(large.effusivity_j_m2_k_sqrt_s(), 1e200);
    close(large.diffusivity_m2_s(), 1.0);
}

#[test]
fn g0_g1_g3_compiled_nasa_316_pack_drives_source_bound_transport() {
    let workspace = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let directory = std::env::temp_dir().join(format!(
        "frankensim-thermal-transport-{}-{nonce}",
        std::process::id()
    ));
    std::fs::create_dir(&directory).unwrap();
    let output = directory.join("nasa-316.fsmatpk");
    let result = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(workspace.join("data/matdb/seed-v1/stainless-316-nasa-tp216435/manifest.tsv"))
        .arg("--out")
        .arg(&output)
        .env("CARGO_WORKSPACE_DIR", workspace)
        .current_dir(workspace)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let bytes = std::fs::read(output).unwrap();
    let decoded = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap();
    assert_eq!(pack.to_bytes(), bytes);
    let normalized = NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: "generic 316 stainless steel; NASA Table 23".into(),
            phase: "solid".into(),
            process: "original source conditions; transport reduction".into(),
            revision: 0,
        },
        pack,
    )
    .unwrap();
    let card = normalized.card();
    let context = |t| {
        [
            ("source-generic-316-identity", 1.0),
            ("source-uns-designation-known", 0.0),
            ("source-heat-treatment-known", 0.0),
            ("source-product-form-known", 0.0),
            ("source-pressure-known", 0.0),
        ]
        .iter()
        .fold(point(t), |p, &(axis, value)| {
            p.with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
                .unwrap()
        })
    };
    let keys = keys();
    let prepared = prepare(card, &keys, &context(273.0), &context(400.0)).unwrap();
    assert_eq!(
        prepared.temperatures_k().collect::<Vec<_>>(),
        [273.0, 300.0, 400.0]
    );
    assert_eq!(prepared.material_card_identity(), card.content_hash());
    for (t, rho, cp, k) in [
        (273.0, 7972.70, 456.39, 12.99),
        (300.0, 7960.74, 467.69, 13.46),
        (350.0, 7938.60, 485.28, 14.30),
        (400.0, 7916.46, 502.87, 15.14),
    ] {
        let sample = prepared.at_temperature(t).unwrap();
        close(sample.density_kg_m3(), rho);
        close(sample.specific_heat_j_kg_k(), cp);
        close(sample.conductivity_w_m_k(), k);
        close(sample.volumetric_heat_capacity_j_m3_k(), rho * cp);
        close(sample.diffusivity_m2_s(), k / (rho * cp));
        close(sample.effusivity_j_m2_k_sqrt_s(), (k * rho * cp).sqrt());
    }
    for line in prepared.source_lines() {
        for (_, answer) in &line.knots {
            assert!(answer.receipt.observation_backed);
            assert!(!answer.receipt.source_hashes.is_empty());
            assert_eq!(
                answer.evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
            card.claims().verify_receipt(&answer.receipt).unwrap();
        }
    }
    assert!(prepare(card, &keys, &context(272.0), &context(400.0)).is_err());
    assert!(
        prepare(
            card,
            &keys,
            &context(273.0).with("source-pressure-known", 1.0).unwrap(),
            &context(400.0).with("source-pressure-known", 1.0).unwrap()
        )
        .is_err()
    );
}
