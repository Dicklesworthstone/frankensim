//! G0/G1/G3: primary acrylic/lumber HFM curves reach bounded thermal consumers.
#![deny(unsafe_code)]

use fs_matdb::{
    ClaimSelection, MaterialCard, MaterialStateId, NormalizedMaterialCardPack, NormalizedPack,
    PropertyKey, PropertyValue, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::thermal::transport::{
    PreparedThermalTransport, ThermalTransportSourceQuery,
    volumetric::{PreparedVolumetricThermalTransport, VolumetricThermalPathQuery},
};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use std::{path::Path, process::Command};

const C: Dims = Dims([-1, 1, -2, -1, 0, 0]);
const CP: Dims = Dims([2, 0, -2, -1, 0, 0]);
const K: Dims = Dims([1, 1, -3, -1, 0, 0]);
const RHO: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const POLICY: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);

struct Material {
    id: &'static str,
    year: u16,
    stock_m: Option<f64>,
    wood: bool,
    densities: [f64; 3],
    capacity: [[f64; 4]; 3],
    cp: [f64; 4],
    cp_sd: [f64; 4],
    k_readings: [[f64; 9]; 2],
    k: [f64; 2],
    k_sd: [f64; 2],
    k_high_c: f64,
}

// Literal final Results Table rows from the pinned primary .tst files and
// reported density CSV columns. Repeated conductivity cycles are not specimens.
const MATERIALS: [Material; 4] = [
    Material {
        id: "black-pmma",
        year: 2021,
        stock_m: Some(0.0048),
        wood: false,
        densities: [1188.0, 1177.0, 1181.0],
        capacity: [
            [1578187.0, 1677591.0, 1696192.0, 1743962.0],
            [1582106.0, 1632241.0, 1682553.0, 1735768.0],
            [1613127.0, 1643995.0, 1703303.0, 1761086.0],
        ],
        cp: [1346.0, 1397.0, 1433.0, 1478.0],
        cp_sd: [16.0, 20.0, 9.0, 11.0],
        k_readings: [
            [
                0.1471, 0.1482, 0.1484, 0.1485, 0.1493, 0.1487, 0.1538, 0.1543, 0.1541,
            ],
            [
                0.1528, 0.1528, 0.1530, 0.1536, 0.1531, 0.1533, 0.1583, 0.1579, 0.1582,
            ],
        ],
        k: [0.150, 0.155],
        k_sd: [0.003, 0.003],
        k_high_c: 45.0,
    },
    Material {
        id: "wood-stud",
        year: 2021,
        stock_m: None,
        wood: true,
        densities: [359.0, 370.0, 408.0],
        capacity: [
            [463087.0, 497836.0, 516761.0, 536194.0],
            [484980.0, 509238.0, 534077.0, 554226.0],
            [510984.0, 537051.0, 559843.0, 580115.0],
        ],
        cp: [1283.0, 1358.0, 1417.0, 1469.0],
        cp_sd: [63.0, 53.0, 57.0, 58.0],
        k_readings: [
            [
                0.08965, 0.09023, 0.09032, 0.09280, 0.09291, 0.09296, 0.1028, 0.1032, 0.1032,
            ],
            [
                0.09493, 0.09500, 0.09507, 0.09735, 0.09761, 0.09776, 0.1083, 0.1085, 0.1085,
            ],
        ],
        k: [0.095, 0.100],
        k_sd: [0.006, 0.006],
        k_high_c: 45.0,
    },
    Material {
        id: "pallet-wood",
        year: 2022,
        stock_m: None,
        wood: true,
        densities: [610.9, 610.0, 544.5],
        capacity: [
            [901123.0, 921574.0, 987232.0, 1004170.0],
            [790563.0, 830533.0, 866661.0, 902352.0],
            [620305.0, 646630.0, 675877.0, 703680.0],
        ],
        cp: [1310.0, 1359.0, 1433.0, 1479.0],
        cp_sd: [240.0, 238.0, 267.0, 260.0],
        k_readings: [
            [
                0.1193, 0.1205, 0.1204, 0.1174, 0.1183, 0.1185, 0.1023, 0.1028, 0.1028,
            ],
            [
                0.1331, 0.1333, 0.1331, 0.1296, 0.1298, 0.1296, 0.1143, 0.1140, 0.1138,
            ],
        ],
        k: [0.114, 0.126],
        k_sd: [0.008, 0.009],
        k_high_c: 65.0,
    },
    Material {
        id: "basswood-panel",
        year: 2021,
        stock_m: Some(0.019),
        wood: true,
        densities: [398.0, 409.0, 405.0],
        capacity: [
            [546690.0, 567279.0, 593691.0, 616682.0],
            [561349.0, 588942.0, 612195.0, 635139.0],
            [553250.0, 588180.0, 603472.0, 630820.0],
        ],
        cp: [1371.0, 1439.0, 1493.0, 1553.0],
        cp_sd: [18.0, 30.0, 23.0, 24.0],
        k_readings: [
            [
                0.09484, 0.09532, 0.09551, 0.09537, 0.09570, 0.09576, 0.09541, 0.09582, 0.09596,
            ],
            [
                0.1006, 0.1008, 0.1009, 0.1010, 0.1010, 0.1011, 0.1008, 0.1010, 0.1011,
            ],
        ],
        k: [0.096, 0.101],
        k_sd: [0.0, 0.0],
        k_high_c: 45.0,
    },
];

fn compile(m: &Material) -> NormalizedMaterialCardPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let slug = format!("fsri-{}-unconditioned-{}", m.id, m.year);
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let output = std::env::temp_dir().join(format!(
        "fsri-acrylic-lumber-{}-{nonce}-{}.fsmatpk",
        std::process::id(),
        m.id
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
    let raw = NormalizedPack::from_bytes(&bytes).unwrap();
    let raw = NormalizedPack::from_bytes_verified(raw.content_hash(), &bytes).unwrap();
    assert_eq!(raw.to_bytes(), bytes);
    assert_eq!(raw.pack_id(), slug);
    assert_eq!(raw.claims().claims_ordered().count(), 8);
    let values: usize = raw
        .claims()
        .claims_ordered()
        .map(|(_, c)| match &c.value {
            PropertyValue::Scalar { .. } => 1,
            PropertyValue::Curve { knots, .. } => knots.len(),
            _ => panic!("unexpected source payload"),
        })
        .sum();
    assert_eq!(values, 23);
    let pack = NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: slug,
            phase: "solid".into(),
            process:
                "FSRI unconditioned HFM source population with unknown grade and material axes"
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
fn context(m: &Material, statistic: f64, local: bool) -> QueryPoint {
    let mut p = QueryPoint::new();
    for (axis, value) in [
        (format!("source-fsri-{}", m.id), 1.0),
        ("source-unconditioned".into(), 1.0),
        ("source-grade-known".into(), 0.0),
        ("source-formulation-known".into(), 0.0),
        ("source-moisture-content-known".into(), 0.0),
        ("source-contact-pressure-known".into(), 0.0),
        ("source-material-axis-known".into(), 0.0),
        (
            "source-metadata-stock-thickness-known".into(),
            if m.stock_m.is_some() { 1.0 } else { 0.0 },
        ),
        ("source-statistic-is-mean".into(), statistic),
        ("source-specimen-count".into(), 3.0),
        ("source-hfm-through-thickness".into(), 1.0),
        ("source-hfm-upward-heat-flow".into(), 1.0),
        ("source-density-temperature-known".into(), 0.0),
        (
            "source-density-header-consistent".into(),
            if m.id == "pallet-wood" { 0.0 } else { 1.0 },
        ),
    ] {
        p = p
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    if local {
        p = p
            .with_quantity(
                "mean-test-temperature-as-local-model",
                QuantitySpec::dimensional(Dims::NONE),
                1.0,
            )
            .unwrap();
    }
    if m.wood {
        for axis in ["source-botanical-species-known", "source-treatment-known"] {
            p = p
                .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), 0.0)
                .unwrap();
        }
    }
    if let Some(thickness) = m.stock_m {
        p = p
            .with_quantity(
                "source-metadata-stock-thickness",
                QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0])),
                thickness,
            )
            .unwrap();
    }
    p.with_quantity(
        "source-hfm-temperature-difference",
        semantic(QuantityKind::TemperatureDifference),
        20.0,
    )
    .unwrap()
}
fn at(p: QueryPoint, axis: &str, kelvin: f64) -> QueryPoint {
    p.with_quantity(axis, semantic(QuantityKind::AbsoluteTemperature), kelvin)
        .unwrap()
}
fn point(m: &Material, kelvin: f64) -> QueryPoint {
    at(context(m, 1.0, true), "temperature", kelvin)
}
fn key(name: &str, dims: Dims) -> PropertyKey {
    if matches!(
        name,
        "thermal-conductivity" | "apparent-thermal-conductivity"
    ) {
        PropertyKey::with_quantity(name, semantic(QuantityKind::ThermalConductivity))
    } else {
        PropertyKey::new(name, dims)
    }
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-12 * expected.abs().max(1e-10),
        "{actual} != {expected}"
    );
}
fn mean(a: &[f64]) -> f64 {
    a.iter().sum::<f64>() / a.len() as f64
}
fn sd(a: &[f64]) -> f64 {
    let average = mean(a);
    (a.iter().map(|x| (x - average) * (x - average)).sum::<f64>() / (a.len() - 1) as f64).sqrt()
}
fn check(card: &MaterialCard, name: &str, dims: Dims, p: &QueryPoint, expected: f64) {
    let answer = card
        .claims()
        .query_typed(&key(name, dims), p, SelectionPolicy::SingleClaimOnly)
        .unwrap();
    close(answer.evidence.value.value, expected);
    assert_eq!(answer.evidence.value.dims, dims);
    assert_eq!(
        answer.evidence.value.uncertainty,
        UncertaintyModel::Unstated
    );
    assert!(answer.receipt.observation_backed);
    assert!(!answer.receipt.source_hashes.is_empty());
    card.claims().verify_receipt(&answer.receipt).unwrap();
}
fn refuses(card: &MaterialCard, name: &str, dims: Dims, p: &QueryPoint) {
    assert!(
        card.claims()
            .query_typed(&key(name, dims), p, SelectionPolicy::SingleClaimOnly)
            .is_err(),
        "must refuse {name} at {p:?}"
    );
}

#[test]
fn g1_g3_native_data_and_published_reductions_preserve_units_and_receipts() {
    for m in &MATERIALS {
        let pack = compile(m);
        let card = pack.card();
        let density = mean(&m.densities);
        check(
            card,
            "reference-normalization-density",
            RHO,
            &context(m, 1.0, false),
            density,
        );
        for j in 0..4 {
            let t = 283.15 + 10.0 * j as f64;
            let capacity = m.capacity.map(|r| r[j]);
            let (volume, scatter) = (mean(&capacity), sd(&capacity));
            check(card, "volumetric-heat-capacity", C, &point(m, t), volume);
            check(card, "specific-heat-capacity", CP, &point(m, t), m.cp[j]);
            let statistics = at(context(m, 0.0, false), "temperature", t);
            check(
                card,
                "computed-volumetric-heat-capacity-sample-standard-deviation",
                C,
                &statistics,
                scatter,
            );
            check(
                card,
                "reported-specific-heat-capacity-standard-deviation",
                CP,
                &statistics,
                m.cp_sd[j],
            );
            assert_eq!((volume / density).round_ties_even(), m.cp[j]);
            assert_eq!((scatter / density).round_ties_even(), m.cp_sd[j]);
        }
        for j in 0..2 {
            let t = 273.15 + [15.0, m.k_high_c][j];
            check(card, "thermal-conductivity", K, &point(m, t), m.k[j]);
            check(
                card,
                "apparent-thermal-conductivity",
                K,
                &at(context(m, 1.0, false), "test-mean-temperature", t),
                m.k[j],
            );
            check(
                card,
                "reported-apparent-thermal-conductivity-standard-deviation",
                K,
                &at(context(m, 0.0, false), "test-mean-temperature", t),
                m.k_sd[j],
            );
            assert_eq!(
                (1000.0 * mean(&m.k_readings[j])).round_ties_even() / 1000.0,
                m.k[j]
            );
            assert_eq!(
                (1000.0 * sd(&m.k_readings[j])).round_ties_even() / 1000.0,
                m.k_sd[j]
            );
            assert!(sd(&m.k_readings[j]) > 0.0); // A published zero is rounding.
        }
        if m.id == "wood-stud" {
            assert_ne!(density, 355.0);
        }
        if m.id == "pallet-wood" {
            assert_ne!(density, 480.0);
        }
    }
}

#[test]
fn g0_identity_unknown_conditions_and_apparent_temperature_semantics_remain_required() {
    for m in &MATERIALS {
        let pack = compile(m);
        let card = pack.card();
        for (axis, value) in [
            ("source-unconditioned", 0.0),
            ("source-grade-known", 1.0),
            ("source-formulation-known", 1.0),
            ("source-moisture-content-known", 1.0),
            ("source-contact-pressure-known", 1.0),
            ("source-material-axis-known", 1.0),
            ("source-specimen-count", 9.0),
        ] {
            refuses(
                card,
                "volumetric-heat-capacity",
                C,
                &point(m, 298.15).with(axis, value).unwrap(),
            );
        }
        if m.wood {
            for axis in ["source-botanical-species-known", "source-treatment-known"] {
                refuses(
                    card,
                    "volumetric-heat-capacity",
                    C,
                    &point(m, 298.15).with(axis, 1.0).unwrap(),
                );
            }
        }
        match m.stock_m {
            Some(t) => refuses(
                card,
                "volumetric-heat-capacity",
                C,
                &point(m, 298.15)
                    .with("source-metadata-stock-thickness", 1.1 * t)
                    .unwrap(),
            ),
            None => refuses(
                card,
                "volumetric-heat-capacity",
                C,
                &point(m, 298.15)
                    .with("source-metadata-stock-thickness-known", 1.0)
                    .unwrap(),
            ),
        }
        for (axis, value) in [
            ("mean-test-temperature-as-local-model", 0.0),
            ("source-hfm-through-thickness", 0.0),
            ("source-hfm-upward-heat-flow", 0.0),
            ("source-hfm-temperature-difference", 25.0),
        ] {
            refuses(
                card,
                "thermal-conductivity",
                K,
                &point(m, 298.15).with(axis, value).unwrap(),
            );
        }
        refuses(
            card,
            "thermal-conductivity",
            K,
            &at(context(m, 1.0, false), "temperature", 298.15),
        );
        refuses(card, "apparent-thermal-conductivity", K, &point(m, 298.15));
        for t in [283.14, 313.16] {
            refuses(card, "volumetric-heat-capacity", C, &point(m, t));
        }
        for t in [288.14, 273.16 + m.k_high_c] {
            refuses(card, "thermal-conductivity", K, &point(m, t));
        }
        refuses(card, "volumetric-heat-capacity", CP, &point(m, 298.15));
        refuses(
            card,
            "volumetric-heat-capacity",
            C,
            &point(m, 298.15)
                .with_quantity(
                    "temperature",
                    semantic(QuantityKind::TemperatureDifference),
                    298.15,
                )
                .unwrap(),
        );
        refuses(
            card,
            "computed-volumetric-heat-capacity-sample-standard-deviation",
            C,
            &at(context(m, 0.0, false), "temperature", 298.15),
        );
        refuses(
            card,
            "reference-normalization-density",
            RHO,
            &context(m, 1.0, false)
                .with("source-density-temperature-known", 1.0)
                .unwrap(),
        );
        if m.id == "pallet-wood" {
            // The original inconsistent CSV headers are not silently repaired.
            refuses(
                card,
                "specific-heat-capacity",
                CP,
                &point(m, 293.15)
                    .with("source-density-header-consistent", 1.0)
                    .unwrap(),
            );
            refuses(
                card,
                "reference-normalization-density",
                RHO,
                &context(m, 1.0, false)
                    .with("source-density-header-consistent", 1.0)
                    .unwrap(),
            );
        }
        for name in ["density", "young-modulus", "thermal-expansion-coefficient"] {
            assert!(card.claims().claims_for(name).is_empty());
        }
        assert!(
            PreparedThermalTransport::try_new(
                card,
                ThermalTransportSourceQuery {
                    density: &key("reference-normalization-density", RHO),
                    specific_heat: &key("specific-heat-capacity", CP),
                    conductivity: &key("thermal-conductivity", K),
                    temperature_axis: "temperature",
                    lower: &point(m, 288.15),
                    upper: &point(m, 313.15),
                    selections: [POLICY; 3],
                }
            )
            .is_err()
        );
    }
}

#[test]
fn g1_sourced_native_capacity_drives_transport_heat_storage_inverse_and_sensible_entropy() {
    for m in &MATERIALS {
        let pack = compile(m);
        let card = pack.card();
        let capacity = key("volumetric-heat-capacity", C);
        let conductivity = key("thermal-conductivity", K);
        let lower = point(m, 288.15);
        let upper = point(m, 313.15);
        let query = VolumetricThermalPathQuery {
            volumetric_heat_capacity: &capacity,
            conductivity: &conductivity,
            temperature_axis: "temperature",
            lower: &lower,
            upper: &upper,
            selections: [POLICY; 2],
        };
        let prepared = PreparedVolumetricThermalTransport::try_new(card, query).unwrap();
        let original = prepared.clone();
        assert_eq!(prepared.temperature_span_k(), [288.15, 313.15]);
        let c: [f64; 4] = std::array::from_fn(|j| mean(&m.capacity.map(|r| r[j])));
        let c25 = 0.5 * (c[1] + c[2]);
        let k25 = m.k[0] + (m.k[1] - m.k[0]) * 10.0 / (m.k_high_c - 15.0);
        let sample = prepared.at_temperature(298.15).unwrap();
        close(sample.volumetric_heat_capacity_j_m3_k(), c25);
        close(sample.conductivity_w_m_k(), k25);
        close(sample.diffusivity_m2_s(), k25 / c25);
        close(sample.effusivity_j_m2_k_sqrt_s(), (k25 * c25).sqrt());
        let c15 = 0.5 * (c[0] + c[1]);
        let storage = 2.5 * (c15 + c[1]) + 5.0 * (c[1] + c[2]) + 5.0 * (c[2] + c[3]);
        close(
            prepared.storage_change_j_m3(288.15, 313.15).unwrap(),
            storage,
        );
        close(
            prepared.storage_change_j_m3(313.15, 288.15).unwrap(),
            -storage,
        );
        close(
            prepared
                .temperature_after_storage_change(288.15, storage)
                .unwrap(),
            313.15,
        );
        close(
            prepared
                .temperature_after_storage_change(313.15, -storage)
                .unwrap(),
            288.15,
        );

        // Independent affine-capacity quadratic for 600 J into 2 litres from 20 C.
        // No density normalization or rounded Cp enters this physical response.
        let slope = (c[2] - c[1]) / 10.0;
        let per_volume = 600.0 / 0.002;
        let delta = 2.0 * per_volume / (c[1] + (c[1] * c[1] + 2.0 * slope * per_volume).sqrt());
        close(
            prepared
                .temperature_after_storage_change(293.15, per_volume)
                .unwrap(),
            293.15 + delta,
        );
        let entropy = slope * 10.0 + (c[1] - slope * 293.15) * (303.15_f64 / 293.15).ln();
        close(
            prepared
                .sensible_entropy_change_j_m3_k(293.15, 303.15)
                .unwrap(),
            entropy,
        );
        close(
            prepared
                .sensible_entropy_change_j_m3_k(303.15, 293.15)
                .unwrap(),
            -entropy,
        );
        assert!(
            prepared
                .temperature_after_storage_change(288.15, 1.01 * storage)
                .is_err()
        );
        assert!(
            PreparedVolumetricThermalTransport::try_new(
                card,
                VolumetricThermalPathQuery {
                    lower: &point(m, 283.15),
                    ..query
                }
            )
            .is_err()
        );
        assert!(
            PreparedVolumetricThermalTransport::try_new(
                card,
                VolumetricThermalPathQuery {
                    upper: &point(m, 318.15),
                    ..query
                }
            )
            .is_err()
        );
        assert!(
            PreparedVolumetricThermalTransport::try_new(
                card,
                VolumetricThermalPathQuery {
                    conductivity: &key("apparent-thermal-conductivity", K),
                    ..query
                }
            )
            .is_err()
        );
        assert_eq!(prepared, original);
        for source in prepared.source_lines() {
            for (_, answer) in &source.knots {
                card.claims().verify_receipt(&answer.receipt).unwrap();
            }
        }
    }
}
