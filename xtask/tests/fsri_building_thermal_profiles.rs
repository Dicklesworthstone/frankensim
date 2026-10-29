//! Native FSRI construction thermal curves, including density-free fiberglass.
//! G0/G1/G3: original reductions, typed source support and caloric consumer use.
#![deny(unsafe_code)]

use fs_matdb::{
    ClaimSelection, MaterialCard, MaterialStateId, NormalizedMaterialCardPack, NormalizedPack,
    PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::thermal::transport::volumetric::{
    PreparedVolumetricThermalTransport, VolumetricThermalPathQuery,
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

struct Building {
    id: &'static str,
    year: u16,
    stock: f64,
    capacity: [[f64; 4]; 3],
    density: Option<[f64; 3]>,
    cp: Option<([f64; 4], [f64; 4])>,
    k_readings: [[f64; 9]; 2],
    k: [f64; 2],
    k_sd: [f64; 2],
    hot_k: f64,
}

// Literal primary result rows. Mineral wool keeps the final R1/R2 files selected
// by FSRI's sorted-column replacement, not a five-specimen pooled estimate.
// Gypsum normalization uses the published density column despite its documented
// metadata/header/mass inconsistencies. Native C_vol does not use that density.
const BUILDINGS: [Building; 4] = [
    Building {
        id: "fiber-cement-board",
        year: 2021,
        stock: 0.0085,
        capacity: [
            [1171294.0, 1204310.0, 1239544.0, 1269595.0],
            [1216015.0, 1268076.0, 1304073.0, 1345828.0],
            [1237398.0, 1271615.0, 1312065.0, 1340740.0],
        ],
        density: Some([1359.0, 1338.0, 1301.0]),
        cp: Some(([907.0, 936.0, 964.0, 990.0], [25.0, 28.0, 30.0, 32.0])),
        k_readings: [
            [
                0.1138, 0.1148, 0.1151, 0.1453, 0.1467, 0.1469, 0.1484, 0.1504, 0.1506,
            ],
            [
                0.1291, 0.1294, 0.1295, 0.1626, 0.1627, 0.1622, 0.1662, 0.1674, 0.1675,
            ],
        ],
        k: [0.137, 0.153],
        k_sd: [0.017, 0.018],
        hot_k: 318.15,
    },
    Building {
        id: "gypsum-wallboard",
        year: 2021,
        stock: 0.0128,
        capacity: [
            [616163.0, 633908.0, 651825.0, 656377.0],
            [623206.0, 638616.0, 655301.0, 662155.0],
            [592216.0, 611031.0, 624168.0, 631953.0],
        ],
        density: Some([595.0, 591.0, 588.0]),
        cp: Some(([1032.0, 1062.0, 1089.0, 1099.0], [27.0, 25.0, 29.0, 27.0])),
        k_readings: [
            [
                0.1454, 0.1469, 0.1470, 0.1584, 0.1594, 0.1589, 0.1497, 0.1513, 0.1510,
            ],
            [
                0.1485, 0.1484, 0.1484, 0.1595, 0.1592, 0.1588, 0.1519, 0.1517, 0.1514,
            ],
        ],
        k: [0.152, 0.153],
        k_sd: [0.006, 0.005],
        hot_k: 318.15,
    },
    Building {
        id: "mineral-wool",
        year: 2021,
        stock: 0.089,
        capacity: [
            [38899.0, 41780.0, 42989.0, 43334.0],
            [40529.0, 43315.0, 44883.0, 45161.0],
            [42288.0, 43594.0, 44713.0, 45601.0],
        ],
        density: Some([50.0, 54.0, 60.0]),
        cp: Some(([742.0, 785.0, 808.0, 818.0], [31.0, 18.0, 19.0, 22.0])),
        k_readings: [
            [
                0.03371, 0.03362, 0.03364, 0.03330, 0.03331, 0.03330, 0.03296, 0.03291, 0.03292,
            ],
            [
                0.03862, 0.03864, 0.03861, 0.03827, 0.03824, 0.03820, 0.03778, 0.03775, 0.03775,
            ],
        ],
        k: [0.033, 0.038],
        k_sd: [0.0, 0.0],
        hot_k: 318.15,
    },
    Building {
        id: "frp-panel",
        year: 2022,
        stock: 0.0023,
        capacity: [
            [1472292.0, 1546098.0, 1653830.0, 1745092.0],
            [1450575.0, 1537905.0, 1628139.0, 1741910.0],
            [1482404.0, 1552146.0, 1666615.0, 1766131.0],
        ],
        density: None,
        cp: None,
        k_readings: [
            [
                0.09268, 0.09417, 0.09364, 0.08412, 0.08819, 0.08866, 0.09301, 0.09597, 0.09592,
            ],
            [
                0.1149, 0.1148, 0.1145, 0.1059, 0.1058, 0.1058, 0.1155, 0.1155, 0.1154,
            ],
        ],
        k: [0.092, 0.112],
        k_sd: [0.004, 0.005],
        hot_k: 338.15,
    },
];

fn compile(b: &Building) -> NormalizedMaterialCardPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let slug = format!("fsri-{}-unconditioned-{}", b.id, b.year);
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let output = std::env::temp_dir().join(format!(
        "fsri-construction-{}-{nonce}-{}.fsmatpk",
        std::process::id(),
        b.id
    ));
    let result = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(root.join(format!("data/matdb/seed-v1/{slug}/manifest.tsv")))
        .arg("--out")
        .arg(&output)
        .env("CARGO_WORKSPACE_DIR", root)
        .current_dir(root)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{slug}: {}",
        String::from_utf8_lossy(&result.stderr)
    );
    let bytes = std::fs::read(output).unwrap();
    let pack = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(pack.content_hash(), &bytes).unwrap();
    assert_eq!(pack.to_bytes(), bytes);
    assert_eq!(pack.pack_id(), slug);
    assert_eq!(
        pack.claims().claims_ordered().count(),
        if b.density.is_some() { 8 } else { 5 }
    );
    let count: usize = pack
        .claims()
        .claims_ordered()
        .map(|(_, c)| match &c.value {
            fs_matdb::PropertyValue::Scalar { .. } => 1,
            fs_matdb::PropertyValue::Curve { knots, .. } => knots.len(),
            _ => panic!("unexpected source representation"),
        })
        .sum();
    assert_eq!(count, if b.density.is_some() { 23 } else { 14 });
    let card = NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: slug,
            phase: "solid".into(),
            process: "FSRI unconditioned construction product; exact HFM source context".into(),
            revision: 0,
        },
        pack,
    )
    .unwrap();
    let bytes = card.to_bytes();
    let restored =
        NormalizedMaterialCardPack::from_bytes_verified(card.content_hash(), &bytes).unwrap();
    assert_eq!(restored.to_bytes(), bytes);
    restored
}
fn semantic(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn context(b: &Building, statistic: f64, local: bool) -> QueryPoint {
    let mut q = QueryPoint::new();
    for (axis, value) in [
        (format!("source-fsri-{}", b.id), 1.0),
        ("source-unconditioned".into(), 1.0),
        ("source-grade-known".into(), 0.0),
        ("source-formulation-known".into(), 0.0),
        ("source-moisture-content-known".into(), 0.0),
        ("source-contact-pressure-known".into(), 0.0),
        ("source-statistic-is-mean".into(), statistic),
        ("source-specimen-count".into(), 3.0),
        ("source-hfm-through-thickness".into(), 1.0),
        ("source-hfm-upward-heat-flow".into(), 1.0),
        ("source-density-temperature-known".into(), 0.0),
    ] {
        q = q
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    if local {
        q = q
            .with_quantity(
                "mean-test-temperature-as-local-model",
                QuantitySpec::dimensional(Dims::NONE),
                1.0,
            )
            .unwrap();
    }
    q.with_quantity(
        "source-stock-thickness",
        QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0])),
        b.stock,
    )
    .unwrap()
    .with_quantity(
        "source-hfm-temperature-difference",
        semantic(QuantityKind::TemperatureDifference),
        20.0,
    )
    .unwrap()
}
fn at(q: QueryPoint, axis: &str, t: f64) -> QueryPoint {
    q.with_quantity(axis, semantic(QuantityKind::AbsoluteTemperature), t)
        .unwrap()
}
fn point(b: &Building, t: f64) -> QueryPoint {
    at(context(b, 1.0, true), "temperature", t)
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
fn mean(x: &[f64]) -> f64 {
    x.iter().sum::<f64>() / x.len() as f64
}
fn sd(x: &[f64]) -> f64 {
    let m = mean(x);
    (x.iter().map(|v| (v - m).powi(2)).sum::<f64>() / (x.len() - 1) as f64).sqrt()
}
fn sample(card: &MaterialCard, name: &str, dims: Dims, q: &QueryPoint, expected: f64) {
    let a = card
        .claims()
        .query_typed(&key(name, dims), q, SelectionPolicy::SingleClaimOnly)
        .unwrap();
    close(a.evidence.value.value, expected);
    assert_eq!(a.evidence.value.dims, dims);
    assert_eq!(a.evidence.value.uncertainty, UncertaintyModel::Unstated);
    assert!(a.receipt.observation_backed);
    assert!(!a.receipt.source_hashes.is_empty());
    card.claims().verify_receipt(&a.receipt).unwrap();
}
fn refuses(card: &MaterialCard, name: &str, dims: Dims, q: &QueryPoint) {
    assert!(
        card.claims()
            .query_typed(&key(name, dims), q, SelectionPolicy::SingleClaimOnly)
            .is_err(),
        "accepted {name} at {q:?}"
    );
}

#[test]
fn g1_g3_all_construction_values_match_primary_reductions() {
    for b in &BUILDINGS {
        let pack = compile(b);
        let card = pack.card();
        for i in 0..4 {
            let t = 283.15 + 10.0 * i as f64;
            let original = b.capacity.map(|r| r[i]);
            let volume = mean(&original);
            let scatter = sd(&original);
            sample(card, "volumetric-heat-capacity", C, &point(b, t), volume);
            let stats = at(context(b, 0.0, true), "temperature", t);
            sample(
                card,
                "computed-volumetric-heat-capacity-sample-standard-deviation",
                C,
                &stats,
                scatter,
            );
            if let Some((cp, cp_sd)) = b.cp {
                let density = mean(&b.density.unwrap());
                sample(card, "specific-heat-capacity", CP, &point(b, t), cp[i]);
                sample(
                    card,
                    "reported-specific-heat-capacity-standard-deviation",
                    CP,
                    &stats,
                    cp_sd[i],
                );
                assert_eq!((volume / density).round_ties_even(), cp[i]);
                assert_eq!((scatter / density).round_ties_even(), cp_sd[i]);
            }
        }
        for (i, t) in [288.15, b.hot_k].into_iter().enumerate() {
            sample(card, "thermal-conductivity", K, &point(b, t), b.k[i]);
            sample(
                card,
                "apparent-thermal-conductivity",
                K,
                &at(context(b, 1.0, false), "test-mean-temperature", t),
                b.k[i],
            );
            sample(
                card,
                "reported-apparent-thermal-conductivity-standard-deviation",
                K,
                &at(context(b, 0.0, false), "test-mean-temperature", t),
                b.k_sd[i],
            );
            assert_eq!(
                (mean(&b.k_readings[i]) * 1000.0).round_ties_even() / 1000.0,
                b.k[i]
            );
            assert_eq!(
                (sd(&b.k_readings[i]) * 1000.0).round_ties_even() / 1000.0,
                b.k_sd[i]
            );
            assert!(
                sd(&b.k_readings[i]) > 0.0,
                "rounded zero SD is not zero original scatter"
            );
        }
        if let Some(density) = b.density {
            sample(
                card,
                "reference-normalization-density",
                RHO,
                &context(b, 1.0, true),
                mean(&density),
            );
            let claim = card.claims().claims_for("reference-normalization-density");
            assert!(
                claim
                    .iter()
                    .all(|(_, c)| !c.validity.bounds().contains_key("temperature"))
            );
        } else {
            for name in [
                "specific-heat-capacity",
                "reference-normalization-density",
                "density",
            ] {
                assert!(card.claims().claims_for(name).is_empty());
            }
        }
    }
}

#[test]
fn g0_condition_direction_temperature_and_statistic_refusals() {
    for b in &BUILDINGS {
        let pack = compile(b);
        let card = pack.card();
        for (axis, value) in [
            ("source-unconditioned", 0.0),
            ("source-grade-known", 1.0),
            ("source-moisture-content-known", 1.0),
            ("source-specimen-count", 9.0),
            ("source-stock-thickness", b.stock * 2.0),
        ] {
            refuses(
                card,
                "volumetric-heat-capacity",
                C,
                &point(b, 298.15).with(axis, value).unwrap(),
            );
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
                &point(b, 298.15).with(axis, value).unwrap(),
            );
        }
        refuses(
            card,
            "thermal-conductivity",
            K,
            &at(context(b, 1.0, false), "temperature", 298.15),
        );
        refuses(card, "apparent-thermal-conductivity", K, &point(b, 298.15));
        for t in [283.14, 313.16] {
            refuses(card, "volumetric-heat-capacity", C, &point(b, t));
        }
        for t in [288.14, b.hot_k + 0.01] {
            refuses(card, "thermal-conductivity", K, &point(b, t));
        }
        refuses(card, "volumetric-heat-capacity", CP, &point(b, 298.15));
        refuses(
            card,
            "volumetric-heat-capacity",
            C,
            &point(b, 298.15)
                .with_quantity(
                    "temperature",
                    semantic(QuantityKind::TemperatureDifference),
                    298.15,
                )
                .unwrap(),
        );
        refuses(
            card,
            "thermal-conductivity",
            K,
            &point(b, 298.15)
                .with_quantity(
                    "source-hfm-temperature-difference",
                    semantic(QuantityKind::AbsoluteTemperature),
                    20.0,
                )
                .unwrap(),
        );
        refuses(
            card,
            "computed-volumetric-heat-capacity-sample-standard-deviation",
            C,
            &at(context(b, 0.0, true), "temperature", 298.15),
        );
        refuses(
            card,
            "computed-volumetric-heat-capacity-sample-standard-deviation",
            C,
            &point(b, 283.15),
        );
        for name in [
            "density",
            "thermal-expansion-coefficient",
            "young-modulus",
            "latent-heat",
        ] {
            assert!(card.claims().claims_for(name).is_empty());
        }
    }
}

#[test]
fn g1_g3_construction_curves_drive_density_free_storage_and_transport() {
    for b in &BUILDINGS {
        let pack = compile(b);
        let card = pack.card();
        let (capacity, conductivity) = (
            key("volumetric-heat-capacity", C),
            key("thermal-conductivity", K),
        );
        let (low, high) = (point(b, 288.15), point(b, 313.15));
        let path = PreparedVolumetricThermalTransport::try_new(
            card,
            VolumetricThermalPathQuery {
                volumetric_heat_capacity: &capacity,
                conductivity: &conductivity,
                temperature_axis: "temperature",
                lower: &low,
                upper: &high,
                selections: [POLICY; 2],
            },
        )
        .unwrap();
        assert_eq!(path.material_card_identity(), card.content_hash());
        let c20 = mean(&b.capacity.map(|v| v[1]));
        let c30 = mean(&b.capacity.map(|v| v[2]));
        let c40 = mean(&b.capacity.map(|v| v[3]));
        let c25 = 0.5 * (c20 + c30);
        let k25 = b.k[0] + (b.k[1] - b.k[0]) * 10.0 / (b.hot_k - 288.15);
        let sample = path.at_temperature(298.15).unwrap();
        close(sample.volumetric_heat_capacity_j_m3_k(), c25);
        close(sample.conductivity_w_m_k(), k25);
        close(sample.diffusivity_m2_s(), k25 / c25);
        close(sample.effusivity_j_m2_k_sqrt_s(), (k25 * c25).sqrt());
        // Independent exact affine-capacity integral over 20->30->35 C.
        let expected = 5.0 * (c20 + c30) + 2.5 * (c30 + 0.5 * (c30 + c40));
        close(path.storage_change_j_m3(293.15, 308.15).unwrap(), expected);
        close(
            path.temperature_after_storage_change(293.15, expected)
                .unwrap(),
            308.15,
        );
        close(
            path.temperature_after_storage_change(308.15, -expected)
                .unwrap(),
            293.15,
        );
        // Fixed gross volume 0.002 m3, supplied 600 J, no loss or volume change.
        let final_t = path
            .temperature_after_storage_change(293.15, 300000.0)
            .unwrap();
        let slope = (c30 - c20) / 10.0;
        let dt = 600000.0 / (c20 + (c20 * c20 + 600000.0 * slope).sqrt());
        close(final_t, 293.15 + dt);
        close(
            path.storage_change_j_m3(293.15, final_t).unwrap() * 0.002,
            600.0,
        );
        assert!(path.at_temperature(313.16).is_err());
        assert!(path.temperature_after_storage_change(313.15, 1.0).is_err());
        assert!(
            PreparedVolumetricThermalTransport::try_new(
                card,
                VolumetricThermalPathQuery {
                    volumetric_heat_capacity: &capacity,
                    conductivity: &conductivity,
                    temperature_axis: "temperature",
                    lower: &low,
                    upper: &high.with("source-unconditioned", 0.0).unwrap(),
                    selections: [POLICY; 2],
                }
            )
            .is_err()
        );
        println!(
            "{}: 600 J in 2 L from 20 C reaches {:.10} C",
            b.id,
            final_t - 273.15
        );
    }
}
