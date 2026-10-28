//! G0/G1/G3 source-driven steady heat flow through a uniform slab.
//! Numerical checks of the admitted source interpolation, not experimental validation.
#![deny(unsafe_code)]

use fs_conduction::{
    ResistanceOrigin, ResistanceUncertainty, SeriesThermalResistance, SlabConductivityQuery,
    ThermalResistanceTerm,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, NormalizedMaterialCardPack,
    NormalizedPack, PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint,
    SelectionPolicy, UncertaintyModel,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec, Temperature};
use std::path::Path;
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

const K: Dims = Dims([1, 1, -3, -1, 0, 0]);
fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn key() -> PropertyKey {
    PropertyKey::with_quantity(
        "thermal-conductivity",
        quantity(QuantityKind::ThermalConductivity),
    )
}
fn point(t: f64, flags: &[(&str, f64)]) -> QueryPoint {
    let mut p = QueryPoint::new()
        .with_quantity(
            "temperature",
            quantity(QuantityKind::AbsoluteTemperature),
            t,
        )
        .unwrap();
    for &(axis, value) in flags {
        p = p
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    p
}
fn slab(
    card: &MaterialCard,
    lower: &QueryPoint,
    upper: &QueryPoint,
    length: f64,
    area: f64,
) -> Result<ThermalResistanceTerm, fs_conduction::ConductionError> {
    ThermalResistanceTerm::slab_over_temperature_span_from_card(
        "source slab",
        length,
        area,
        card,
        SlabConductivityQuery {
            property: &key(),
            temperature_axis: "temperature",
            lower,
            upper,
            selection: SelectionPolicy::SingleClaimOnly,
        },
    )
}
fn manufactured(
    value: PropertyValue,
    validity: ValidityDomain,
    interpolation: InterpolationPolicy,
) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: key(),
            value,
            validity,
            interpolation,
            uncertainty: UncertaintyModel::Unstated,
            observations: Vec::new(),
            provenance: Provenance {
                source: "manufactured conductivity for numerical integration; not measured data"
                    .into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured conductivity".into(),
            phase: "solid".into(),
            process: "analytic regression".into(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .unwrap()
}
fn temperature_domain(low: f64, high: f64) -> ValidityDomain {
    ValidityDomain::unconstrained().with_quantity(
        "temperature",
        quantity(QuantityKind::AbsoluteTemperature),
        low,
        high,
    )
}

fn curve(knots: Vec<(f64, f64)>, interpolation: InterpolationPolicy) -> MaterialCard {
    manufactured(
        PropertyValue::Curve {
            abscissa: "temperature".into(),
            abscissa_dims: Temperature::DIMS,
            knots,
            dims: K,
        },
        temperature_domain(280.0, 360.0),
        interpolation,
    )
}

#[test]
fn g1_g3_piecewise_conductivity_controls_heat_flow_and_retains_sources() {
    let card = curve(
        vec![(280.0, 2.0), (300.0, 10.0), (360.0, 4.0)],
        InterpolationPolicy::LinearInside,
    );
    let lower = point(280.0, &[]);
    let upper = point(360.0, &[]);
    let term = slab(&card, &lower, &upper, 0.02, 0.004).unwrap();
    // Independent exact areas:20*(2+10)/2 +60*(10+4)/2 =540 W/m.
    // A/L gives108 W; substituting midpoint k instead gives128 W.
    assert!((term.source_span_heat_rate_w().unwrap() - 108.0).abs() < 1e-12);
    assert!((term.value_k_per_w() - 80.0 / 108.0).abs() < 1e-14);
    assert!((term.source_span_heat_rate_w().unwrap() - 128.0).abs() > 10.0);
    assert_eq!(term.uncertainty(), &ResistanceUncertainty::Unstated);
    let ResistanceOrigin::BulkMaterialTemperatureSpan {
        card_identity,
        conductivity,
        length_m,
        area_m2,
    } = term.origin()
    else {
        panic!("source-integrated provenance required");
    };
    assert_eq!(*card_identity, card.content_hash());
    assert_eq!((*length_m, *area_m2), (0.02, 0.004));
    assert_eq!(conductivity.axis, "temperature");
    assert_eq!(
        conductivity.knots.iter().map(|k| k.0).collect::<Vec<_>>(),
        [280.0, 300.0, 360.0]
    );
    for (_, answer) in &conductivity.knots {
        card.claims().verify_receipt(&answer.receipt).unwrap();
        assert_eq!(
            answer.evidence.value.uncertainty,
            UncertaintyModel::Unstated
        );
    }
    let longer = slab(&card, &lower, &upper, 0.04, 0.004).unwrap();
    let wider = slab(&card, &lower, &upper, 0.02, 0.008).unwrap();
    assert_eq!(
        longer.source_span_heat_rate_w().unwrap(),
        term.source_span_heat_rate_w().unwrap() / 2.0
    );
    assert_eq!(
        wider.source_span_heat_rate_w().unwrap(),
        term.source_span_heat_rate_w().unwrap() * 2.0
    );
    // Partial segments:10*(6+10)/2 +30*(10+7)/2 =335 W/m.
    let partial = slab(&card, &point(290.0, &[]), &point(330.0, &[]), 0.02, 0.004).unwrap();
    assert!((partial.source_span_heat_rate_w().unwrap() - 67.0).abs() < 1e-12);
    let network = SeriesThermalResistance::new(vec![term.clone()]).unwrap();
    assert_eq!(network.budget().value_k_per_w, term.value_k_per_w());
    assert_eq!(network.budget().complete_half_width_k_per_w(), None);
}

#[test]
fn g0_span_refuses_missing_support_discrete_values_and_changed_context() {
    let lower = point(280.0, &[]);
    let upper = point(360.0, &[]);
    let card = curve(
        vec![(280.0, 2.0), (300.0, 10.0), (360.0, 4.0)],
        InterpolationPolicy::LinearInside,
    );
    for bad in [0.0, -1.0, f64::INFINITY, f64::NAN] {
        assert!(slab(&card, &lower, &upper, bad, 0.004).is_err());
        assert!(slab(&card, &lower, &upper, 0.02, bad).is_err());
    }
    assert!(slab(&card, &upper, &lower, 0.02, 0.004).is_err());
    assert!(slab(&card, &lower, &lower, 0.02, 0.004).is_err());
    assert!(slab(&card, &point(279.0, &[]), &upper, 0.02, 0.004).is_err());
    assert!(slab(&card, &lower, &point(361.0, &[]), 0.02, 0.004).is_err());
    assert!(
        slab(
            &card,
            &point(280.0, &[("condition", 0.0)]),
            &point(360.0, &[("condition", 1.0)]),
            0.02,
            0.004
        )
        .is_err()
    );
    let difference = QueryPoint::new()
        .with_quantity(
            "temperature",
            quantity(QuantityKind::TemperatureDifference),
            280.0,
        )
        .unwrap();
    assert!(slab(&card, &difference, &upper, 0.02, 0.004).is_err());
    let discrete = curve(
        vec![(280.0, 2.0), (360.0, 4.0)],
        InterpolationPolicy::TabulatedOnly,
    );
    assert!(slab(&discrete, &lower, &upper, 0.02, 0.004).is_err());
    let nonpositive = curve(
        vec![(280.0, 2.0), (300.0, 0.0), (360.0, 4.0)],
        InterpolationPolicy::LinearInside,
    );
    assert!(slab(&nonpositive, &lower, &upper, 0.02, 0.004).is_err());
    let unknown = manufactured(
        PropertyValue::Scalar {
            value: 2.0,
            dims: K,
        },
        ValidityDomain::unconstrained(),
        InterpolationPolicy::ConstantWithinValidity,
    );
    assert!(slab(&unknown, &lower, &upper, 0.02, 0.004).is_err());
    let point_only = manufactured(
        PropertyValue::Scalar {
            value: 2.0,
            dims: K,
        },
        temperature_domain(280.0, 280.0),
        InterpolationPolicy::ConstantWithinValidity,
    );
    assert!(slab(&point_only, &lower, &upper, 0.02, 0.004).is_err());
    let constant = manufactured(
        PropertyValue::Scalar {
            value: 2.0,
            dims: K,
        },
        temperature_domain(280.0, 360.0),
        InterpolationPolicy::ConstantWithinValidity,
    );
    assert!(
        (slab(&constant, &lower, &upper, 0.02, 0.004)
            .unwrap()
            .source_span_heat_rate_w()
            .unwrap()
            - 32.0)
            .abs()
            < 1e-12
    );
    let wrong = PropertyKey::new("specific-heat-capacity", Dims([2, 0, -2, -1, 0, 0]));
    assert!(
        ThermalResistanceTerm::slab_over_temperature_span_from_card(
            "wrong role",
            0.02,
            0.004,
            &card,
            SlabConductivityQuery {
                property: &wrong,
                temperature_axis: "temperature",
                lower: &lower,
                upper: &upper,
                selection: SelectionPolicy::SingleClaimOnly,
            }
        )
        .is_err()
    );
    let reference = ThermalResistanceTerm::slab_from_card(
        "point only",
        0.02,
        0.004,
        &card,
        &key(),
        &lower,
        SelectionPolicy::SingleClaimOnly,
    )
    .unwrap();
    assert!(reference.source_span_heat_rate_w().is_err());
}

fn compile(slug: &str) -> NormalizedMaterialCardPack {
    static NEXT: AtomicU64 = AtomicU64::new(0);
    let workspace = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let directory = std::env::temp_dir().join(format!(
        "frankensim-slab-conductivity-{}-{nonce}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir(&directory).unwrap();
    let output = directory.join("material.fsmatpk");
    let result = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(workspace.join(format!("data/matdb/seed-v1/{slug}/manifest.tsv")))
        .arg("--out")
        .arg(&output)
        .current_dir(workspace)
        .env("CARGO_WORKSPACE_DIR", workspace)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{slug}: {}",
        String::from_utf8_lossy(&result.stderr)
    );
    let bytes = std::fs::read(&output).unwrap();
    let pack = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(pack.content_hash(), &bytes).unwrap();
    assert_eq!(pack.to_bytes(), bytes);
    NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: slug.into(),
            phase: "solid".into(),
            process: "original source conditions; steady slab reduction".into(),
            revision: 0,
        },
        pack,
    )
    .unwrap()
}

#[test]
fn g1_g3_existing_metal_and_insulation_curves_drive_slab_heat_flow() {
    let cases: &[(&str, &[(&str, f64)], &[(f64, f64)], f64)] = &[
        (
            "stainless-316-nasa-tp216435",
            &[
                ("source-generic-316-identity", 1.0),
                ("source-uns-designation-known", 0.0),
                ("source-heat-treatment-known", 0.0),
                ("source-product-form-known", 0.0),
                ("source-pressure-known", 0.0),
            ],
            &[(273.0, 12.99), (300.0, 13.46), (400.0, 15.14)],
            1787.075,
        ),
        (
            "aluminum-2024-t3-first-heating-km2024",
            &[("first-heating", 1.0), ("source-pressure-known", 0.0)],
            &[(298.15, 121.1), (323.15, 125.2), (373.15, 133.0)],
            9533.75,
        ),
        (
            "rockwool-prorox-sl960-en-2020",
            &[
                ("source-grade-prorox-sl960-en-2020", 1.0),
                ("source-temperature-known", 1.0),
                ("source-moisture-condition-known", 0.0),
                ("source-pressure-known", 0.0),
            ],
            &[
                (323.15, 0.040),
                (373.15, 0.045),
                (423.15, 0.052),
                (473.15, 0.060),
                (523.15, 0.071),
                (573.15, 0.081),
                (623.15, 0.094),
                (673.15, 0.107),
                (773.15, 0.140),
                (873.15, 0.175),
                (933.15, 0.200),
            ],
            63.175,
        ),
    ];
    for &(slug, flags, expected, integral_w_m) in cases {
        let normalized = compile(slug);
        let card = normalized.card();
        let lower = point(expected[0].0, flags);
        let upper = point(expected[expected.len() - 1].0, flags);
        let term = slab(card, &lower, &upper, 0.02, 0.004).unwrap();
        let expected_heat = 0.2 * integral_w_m;
        assert!(
            (term.source_span_heat_rate_w().unwrap() / expected_heat - 1.0).abs() < 1e-13,
            "{slug}"
        );
        let ResistanceOrigin::BulkMaterialTemperatureSpan {
            card_identity,
            conductivity,
            ..
        } = term.origin()
        else {
            panic!("source origin required")
        };
        assert_eq!(*card_identity, card.content_hash());
        assert_eq!(conductivity.knots.len(), expected.len());
        for ((temperature, answer), &(t, k)) in conductivity.knots.iter().zip(expected) {
            assert!((temperature - t).abs() < 1e-10);
            assert!((answer.evidence.value.value - k).abs() < 1e-13);
            assert!(answer.receipt.observation_backed);
            assert!(!answer.receipt.source_hashes.is_empty());
            card.claims().verify_receipt(&answer.receipt).unwrap();
        }
        let frozen =
            0.2 * expected[0].1 * (upper.axes()["temperature"] - lower.axes()["temperature"]);
        assert!((term.source_span_heat_rate_w().unwrap() - frozen).abs() > 0.01);
        assert!(
            slab(
                card,
                &lower
                    .clone()
                    .with("temperature", expected[0].0 - 0.01)
                    .unwrap(),
                &upper,
                0.02,
                0.004
            )
            .is_err()
        );
        assert!(
            slab(
                card,
                &lower,
                &upper.clone().with("source-pressure-known", 1.0).unwrap(),
                0.02,
                0.004
            )
            .is_err()
        );
        assert_eq!(term, slab(card, &lower, &upper, 0.02, 0.004).unwrap());
        eprintln!(
            "{slug}: source-integrated heat={:.9} W, resistance={:.9} K/W, retained-knots={}",
            term.source_span_heat_rate_w().unwrap(),
            term.value_k_per_w(),
            conductivity.knots.len()
        );
    }
}
