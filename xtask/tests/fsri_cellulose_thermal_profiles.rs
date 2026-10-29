//! Primary cellulose-insulation measurements on the shared thermal consumers.
#![deny(unsafe_code)]

use fs_matdb::{
    ClaimSelection, MaterialCard, MaterialStateId, NormalizedMaterialCardPack, NormalizedPack,
    PropertyKey, PropertyValue, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::thermal::{
    conduction::{ConductivityCurve, LayerConductionConfig, TwoLayerConduction},
    transport::volumetric::{PreparedVolumetricThermalTransport, VolumetricThermalPathQuery},
};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use std::{path::Path, process::Command};

const SLUG: &str = "fsri-cellulose-insulation-unconditioned-2022";
const C: Dims = Dims([-1, 1, -2, -1, 0, 0]);
const CP: Dims = Dims([2, 0, -2, -1, 0, 0]);
const K: Dims = Dims([1, 1, -3, -1, 0, 0]);
const RHO: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const POLICY: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);
// Independent literal final Results Table values from the three pinned runs.
const CAPACITY: [[f64; 4]; 3] = [
    [124118.0, 129140.0, 135311.0, 140681.0],
    [112244.0, 117979.0, 123255.0, 128034.0],
    [76292.0, 78908.0, 82994.0, 86489.0],
];
const READINGS: [[f64; 9]; 2] = [
    [
        0.04651, 0.04546, 0.04543, 0.04576, 0.04513, 0.04513, 0.04352, 0.04305, 0.04300,
    ],
    [
        0.05052, 0.05005, 0.04991, 0.05006, 0.04985, 0.04980, 0.04883, 0.04866, 0.04861,
    ],
];

fn compiled() -> NormalizedMaterialCardPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let output =
        std::env::temp_dir().join(format!("cellulose-{}-{nonce}.fsmatpk", std::process::id()));
    let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(root.join(format!("data/matdb/seed-v1/{SLUG}/manifest.tsv")))
        .arg("--out")
        .arg(&output)
        .env("CARGO_WORKSPACE_DIR", root)
        .current_dir(root)
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "{}",
        String::from_utf8_lossy(&run.stderr)
    );
    let bytes = std::fs::read(output).unwrap();
    let raw = NormalizedPack::from_bytes(&bytes).unwrap();
    let raw = NormalizedPack::from_bytes_verified(raw.content_hash(), &bytes).unwrap();
    assert_eq!(raw.to_bytes(), bytes);
    assert_eq!(raw.claims().claims_ordered().count(), 8);
    let values: usize = raw
        .claims()
        .claims_ordered()
        .map(|(_, claim)| match &claim.value {
            PropertyValue::Scalar { .. } => 1,
            PropertyValue::Curve { knots, .. } => knots.len(),
            _ => panic!("unexpected source payload"),
        })
        .sum();
    assert_eq!(values, 23);
    let pack = NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: SLUG.into(),
            phase: "solid".into(),
            process: "unconditioned treated paper insulation; unknown packing and coupon match"
                .into(),
            revision: 0,
        },
        raw,
    )
    .unwrap();
    NormalizedMaterialCardPack::from_bytes_verified(pack.content_hash(), &pack.to_bytes()).unwrap()
}
fn semantic(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn context(statistic: f64, local: bool) -> QueryPoint {
    let mut point = QueryPoint::new();
    for (axis, value) in [
        ("source-fsri-cellulose-insulation", 1.0),
        ("source-unconditioned", 1.0),
        ("source-grade-known", 0.0),
        ("source-formulation-known", 0.0),
        ("source-moisture-content-known", 0.0),
        ("source-contact-pressure-known", 0.0),
        ("source-material-axis-known", 0.0),
        ("source-packing-state-known", 0.0),
        ("source-treatment-details-known", 0.0),
        ("source-specimens-matched", 0.0),
        ("source-specimen-count", 3.0),
        ("source-statistic-is-mean", statistic),
        ("source-hfm-through-thickness", 1.0),
        ("source-hfm-upward-heat-flow", 1.0),
        ("source-density-temperature-known", 0.0),
    ] {
        point = point
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    if local {
        point = point
            .with_quantity(
                "mean-test-temperature-as-local-model",
                QuantitySpec::dimensional(Dims::NONE),
                1.0,
            )
            .unwrap();
    }
    point
        .with_quantity(
            "source-hfm-temperature-difference",
            semantic(QuantityKind::TemperatureDifference),
            20.0,
        )
        .unwrap()
}
fn at(point: QueryPoint, axis: &str, kelvin: f64) -> QueryPoint {
    point
        .with_quantity(axis, semantic(QuantityKind::AbsoluteTemperature), kelvin)
        .unwrap()
}
fn point(kelvin: f64) -> QueryPoint {
    at(context(1.0, true), "temperature", kelvin)
}
fn key(name: &str, dims: Dims) -> PropertyKey {
    if name == "thermal-conductivity" || name == "apparent-thermal-conductivity" {
        PropertyKey::with_quantity(name, semantic(QuantityKind::ThermalConductivity))
    } else {
        PropertyKey::new(name, dims)
    }
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-11 * expected.abs().max(1e-10),
        "{actual} != {expected}"
    );
}
fn mean(values: &[f64]) -> f64 {
    values.iter().sum::<f64>() / values.len() as f64
}
fn sd(values: &[f64]) -> f64 {
    let m = mean(values);
    (values.iter().map(|v| (v - m) * (v - m)).sum::<f64>() / (values.len() - 1) as f64).sqrt()
}
fn check(card: &MaterialCard, name: &str, dims: Dims, point: &QueryPoint, expected: f64) {
    let answer = card
        .claims()
        .query_typed(&key(name, dims), point, SelectionPolicy::SingleClaimOnly)
        .unwrap();
    close(answer.evidence.value.value, expected);
    assert_eq!(answer.evidence.value.dims, dims);
    assert_eq!(
        answer.evidence.value.uncertainty,
        UncertaintyModel::Unstated
    );
    assert!(answer.receipt.observation_backed && !answer.receipt.source_hashes.is_empty());
    card.claims().verify_receipt(&answer.receipt).unwrap();
}
fn refuses(card: &MaterialCard, name: &str, dims: Dims, point: &QueryPoint) {
    assert!(
        card.claims()
            .query_typed(&key(name, dims), point, SelectionPolicy::SingleClaimOnly)
            .is_err()
    );
}

#[test]
fn g1_g3_native_population_and_published_reductions_survive_the_real_compiler() {
    let pack = compiled();
    let card = pack.card();
    let density = mean(&[80.0, 71.0, 53.0]);
    check(
        card,
        "reference-normalization-density",
        RHO,
        &context(1.0, false),
        density,
    );
    for j in 0..4 {
        let samples = CAPACITY.map(|row| row[j]);
        let (capacity, scatter) = (mean(&samples), sd(&samples));
        let t = 283.15 + 10.0 * j as f64;
        check(card, "volumetric-heat-capacity", C, &point(t), capacity);
        check(
            card,
            "specific-heat-capacity",
            CP,
            &point(t),
            [1533.0, 1598.0, 1674.0, 1741.0][j],
        );
        check(
            card,
            "computed-volumetric-heat-capacity-sample-standard-deviation",
            C,
            &at(context(0.0, false), "temperature", t),
            scatter,
        );
        check(
            card,
            "reported-specific-heat-capacity-standard-deviation",
            CP,
            &at(context(0.0, false), "temperature", t),
            [366.0, 388.0, 403.0, 417.0][j],
        );
        assert_eq!(
            (capacity / density).round_ties_even(),
            [1533.0, 1598.0, 1674.0, 1741.0][j]
        );
        assert_eq!(
            (scatter / density).round_ties_even(),
            [366.0, 388.0, 403.0, 417.0][j]
        );
        assert!(scatter / capacity > 0.23); // Keep the large reported population variability.
    }
    for j in 0..2 {
        let t = [288.15, 318.15][j];
        let k = [0.045, 0.050][j];
        check(card, "thermal-conductivity", K, &point(t), k);
        check(
            card,
            "apparent-thermal-conductivity",
            K,
            &at(context(1.0, false), "test-mean-temperature", t),
            k,
        );
        check(
            card,
            "reported-apparent-thermal-conductivity-standard-deviation",
            K,
            &at(context(0.0, false), "test-mean-temperature", t),
            0.001,
        );
        assert_eq!((1000.0 * mean(&READINGS[j])).round_ties_even() / 1000.0, k);
        assert_eq!(
            (1000.0 * sd(&READINGS[j])).round_ties_even() / 1000.0,
            0.001
        );
    }
}

#[test]
fn g0_packing_coupons_conditioning_and_temperature_semantics_remain_explicit() {
    let pack = compiled();
    let card = pack.card();
    for (axis, value) in [
        ("source-unconditioned", 0.0),
        ("source-packing-state-known", 1.0),
        ("source-treatment-details-known", 1.0),
        ("source-specimens-matched", 1.0),
        ("source-contact-pressure-known", 1.0),
        ("source-moisture-content-known", 1.0),
        ("source-specimen-count", 9.0),
    ] {
        refuses(
            card,
            "volumetric-heat-capacity",
            C,
            &point(298.15).with(axis, value).unwrap(),
        );
    }
    for (axis, value) in [
        ("mean-test-temperature-as-local-model", 0.0),
        ("source-hfm-through-thickness", 0.0),
        ("source-hfm-temperature-difference", 25.0),
    ] {
        refuses(
            card,
            "thermal-conductivity",
            K,
            &point(298.15).with(axis, value).unwrap(),
        );
    }
    refuses(
        card,
        "thermal-conductivity",
        K,
        &at(context(1.0, false), "temperature", 298.15),
    );
    refuses(card, "apparent-thermal-conductivity", K, &point(298.15));
    refuses(card, "volumetric-heat-capacity", CP, &point(298.15));
    refuses(
        card,
        "volumetric-heat-capacity",
        C,
        &point(298.15)
            .with_quantity(
                "temperature",
                semantic(QuantityKind::TemperatureDifference),
                298.15,
            )
            .unwrap(),
    );
    for t in [283.14, 313.16] {
        refuses(card, "volumetric-heat-capacity", C, &point(t));
    }
    for t in [288.14, 318.16] {
        refuses(card, "thermal-conductivity", K, &point(t));
    }
    refuses(
        card,
        "computed-volumetric-heat-capacity-sample-standard-deviation",
        C,
        &at(context(0.0, false), "temperature", 298.15),
    );
    refuses(
        card,
        "reference-normalization-density",
        RHO,
        &context(1.0, false)
            .with("source-density-temperature-known", 1.0)
            .unwrap(),
    );
    for name in ["density", "young-modulus", "thermal-expansion-coefficient"] {
        assert!(card.claims().claims_for(name).is_empty());
    }
}

#[test]
fn g1_g3_native_capacity_and_sourced_conduction_drive_bounded_thermal_response() {
    let pack = compiled();
    let card = pack.card();
    let ckey = key("volumetric-heat-capacity", C);
    let kkey = key("thermal-conductivity", K);
    let lower = point(288.15);
    let upper = point(313.15);
    let query = VolumetricThermalPathQuery {
        volumetric_heat_capacity: &ckey,
        conductivity: &kkey,
        temperature_axis: "temperature",
        lower: &lower,
        upper: &upper,
        selections: [POLICY; 2],
    };
    let path = PreparedVolumetricThermalTransport::try_new(card, query).unwrap();
    let original = path.clone();
    let c20 = (129140.0 + 117979.0 + 78908.0) / 3.0;
    let c30 = (135311.0 + 123255.0 + 82994.0) / 3.0;
    let slope = (c30 - c20) / 10.0;
    let sample = path.at_temperature(298.15).unwrap();
    close(sample.volumetric_heat_capacity_j_m3_k(), 0.5 * (c20 + c30));
    let k25 = 0.045 + (0.050 - 0.045) * 10.0 / 30.0;
    close(sample.conductivity_w_m_k(), k25);
    close(sample.diffusivity_m2_s(), k25 / (0.5 * (c20 + c30)));
    let storage = 5.0 * (c20 + c30);
    close(path.storage_change_j_m3(293.15, 303.15).unwrap(), storage);
    close(path.storage_change_j_m3(303.15, 293.15).unwrap(), -storage);
    // Fixed source reference volume and compatible sensible-caloric interpretation
    // are declared by this scenario. No installed packing/density law is supplied.
    let delivered = 600.0 / 0.002;
    let delta = 2.0 * delivered / (c20 + (c20 * c20 + 2.0 * slope * delivered).sqrt());
    close(
        path.temperature_after_storage_change(293.15, delivered)
            .unwrap(),
        293.15 + delta,
    );
    let entropy = slope * 10.0 + (c20 - slope * 293.15) * (303.15_f64 / 293.15).ln();
    close(
        path.sensible_entropy_change_j_m3_k(293.15, 303.15).unwrap(),
        entropy,
    );
    let k = ConductivityCurve::new(
        path.source_lines()[1]
            .knots
            .iter()
            .map(|(t, a)| [*t, a.evidence.value.value])
            .collect(),
    )
    .unwrap();
    let layers = TwoLayerConduction::new([&k, &k], [0.01; 2], 0.01, 0.0).unwrap();
    let result = layers
        .evaluate(
            [293.15, 303.15],
            LayerConductionConfig {
                max_iterations: 128,
                absolute_heat_tolerance_w: 1e-12,
                relative_heat_tolerance: 1e-12,
            },
            || Ok(()),
        )
        .unwrap();
    // Independent integral of the affine 15/45 C conductivity over 20..30 C.
    close(result.power_a_to_b_w, -0.01 / 0.02 * 10.0 * k25);
    assert!(
        PreparedVolumetricThermalTransport::try_new(
            card,
            VolumetricThermalPathQuery {
                lower: &point(283.15),
                ..query
            }
        )
        .is_err()
    );
    assert!(
        PreparedVolumetricThermalTransport::try_new(
            card,
            VolumetricThermalPathQuery {
                upper: &point(318.15),
                ..query
            }
        )
        .is_err()
    );
    assert!(path.temperature_after_storage_change(293.15, 1e9).is_err());
    for line in path.source_lines() {
        for (_, answer) in &line.knots {
            card.claims().verify_receipt(&answer.receipt).unwrap();
        }
    }
    assert_eq!(path, original);
}
