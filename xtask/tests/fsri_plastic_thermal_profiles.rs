//! FSRI HFM plastics preserve native volumetric storage, published rounding and source conditions.
#![deny(unsafe_code)]

use fs_matdb::{
    ClaimSelection, MaterialAnswer, MaterialCard, MaterialStateId, NormalizedMaterialCardPack,
    NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel,
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

struct Plastic {
    id: &'static str,
    stock: f64,
    densities: [f64; 3],
    capacity: [[f64; 4]; 3],
    cp: [f64; 4],
    cp_sd: [f64; 4],
    k_readings: [[f64; 9]; 2],
    k: [f64; 2],
    k_sd: [f64; 2],
}

// Native volume-capacity rows and final Average Cond readings from the pinned
// primary .tst files; density inputs are the reported CSV column, not mass/volume.
const PLASTICS: [Plastic; 4] = [
    Plastic {
        id: "hdpe",
        stock: 0.0032,
        densities: [974.0, 970.0, 970.0],
        capacity: [
            [1656208.0, 1729797.0, 1816900.0, 1875391.0],
            [1633251.0, 1716784.0, 1804399.0, 1868124.0],
            [1626676.0, 1702556.0, 1789628.0, 1859554.0],
        ],
        cp: [1687.0, 1767.0, 1857.0, 1923.0],
        cp_sd: [16.0, 14.0, 14.0, 8.0],
        k_readings: [
            [
                0.249, 0.2643, 0.2641, 0.246, 0.2606, 0.2641, 0.2494, 0.2724, 0.2723,
            ],
            [
                0.298, 0.2994, 0.2972, 0.2855, 0.287, 0.2892, 0.303, 0.3034, 0.3016,
            ],
        ],
        k: [0.26, 0.296],
        k_sd: [0.01, 0.007],
    },
    Plastic {
        id: "pvc",
        stock: 0.0033,
        densities: [1385.0, 1388.0, 1391.0],
        capacity: [
            [1350841.0, 1414902.0, 1434708.0, 1495362.0],
            [1339459.0, 1435836.0, 1439617.0, 1507433.0],
            [1355911.0, 1394661.0, 1462926.0, 1478290.0],
        ],
        cp: [972.0, 1020.0, 1042.0, 1076.0],
        cp_sd: [6.0, 15.0, 11.0, 11.0],
        k_readings: [
            [
                0.1287, 0.1301, 0.1293, 0.1313, 0.1322, 0.1333, 0.1273, 0.1302, 0.1301,
            ],
            [
                0.1423, 0.1432, 0.142, 0.1462, 0.144, 0.1443, 0.141, 0.1412, 0.1422,
            ],
        ],
        k: [0.13, 0.143],
        k_sd: [0.002, 0.002],
    },
    Plastic {
        id: "pc",
        stock: 0.0056,
        densities: [1213.0, 1223.0, 1224.0],
        capacity: [
            [1393874.0, 1448894.0, 1496894.0, 1546756.0],
            [1390507.0, 1434027.0, 1489976.0, 1532565.0],
            [1391615.0, 1432167.0, 1498494.0, 1544392.0],
        ],
        cp: [1141.0, 1179.0, 1226.0, 1263.0],
        cp_sd: [1.0, 8.0, 4.0, 6.0],
        k_readings: [
            [
                0.1671, 0.17, 0.1693, 0.1562, 0.157, 0.1571, 0.1627, 0.1655, 0.1646,
            ],
            [
                0.1851, 0.1848, 0.185, 0.1718, 0.1718, 0.1716, 0.1801, 0.1791, 0.1787,
            ],
        ],
        k: [0.163, 0.179],
        k_sd: [0.005, 0.006],
    },
    Plastic {
        id: "petg",
        stock: 0.0027,
        densities: [1319.0, 1312.0, 1311.0],
        capacity: [
            [1481289.0, 1538766.0, 1588049.0, 1643255.0],
            [1466324.0, 1520565.0, 1570005.0, 1619456.0],
            [1485290.0, 1538896.0, 1593949.0, 1645234.0],
        ],
        cp: [1125.0, 1166.0, 1205.0, 1245.0],
        cp_sd: [8.0, 8.0, 9.0, 11.0],
        k_readings: [
            [
                0.1426, 0.1455, 0.1433, 0.1363, 0.1392, 0.141, 0.1451, 0.1502, 0.1509,
            ],
            [
                0.1603, 0.1585, 0.1583, 0.1566, 0.1572, 0.157, 0.1628, 0.1627, 0.1633,
            ],
        ],
        k: [0.144, 0.16],
        k_sd: [0.005, 0.003],
    },
];

fn compile(p: &Plastic) -> NormalizedMaterialCardPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let slug = format!("fsri-{}-unconditioned-2021", p.id);
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let path = std::env::temp_dir().join(format!(
        "fsri-plastic-test-{}-{nonce}-{}.fsmatpk",
        std::process::id(),
        p.id
    ));
    let result = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(root.join(format!("data/matdb/seed-v1/{slug}/manifest.tsv")))
        .arg("--out")
        .arg(&path)
        .env("CARGO_WORKSPACE_DIR", root)
        .current_dir(root)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{slug}: {}",
        String::from_utf8_lossy(&result.stderr)
    );
    let bytes = std::fs::read(path).unwrap();
    let pack = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(pack.content_hash(), &bytes).unwrap();
    assert_eq!(pack.to_bytes(), bytes);
    assert_eq!(pack.pack_id(), slug);
    assert_eq!(pack.claims().claims_ordered().count(), 8);
    let values: usize = pack
        .claims()
        .claims_ordered()
        .map(|(_, claim)| match &claim.value {
            fs_matdb::PropertyValue::Scalar { .. } => 1,
            fs_matdb::PropertyValue::Curve { knots, .. } => knots.len(),
            _ => panic!("unexpected material representation"),
        })
        .sum();
    assert_eq!(values, 23);
    let card = NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: slug,
            phase: "solid".into(),
            process: "FSRI unconditioned commercial sheet; primary HFM experiments in 2021".into(),
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
fn context(p: &Plastic, statistic: f64, local_model: bool) -> QueryPoint {
    let mut point = QueryPoint::new();
    for (name, value) in [
        (format!("source-fsri-{}", p.id), 1.0),
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
        point = point
            .with_quantity(name, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    if local_model {
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
            "source-stock-thickness",
            QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0])),
            p.stock,
        )
        .unwrap()
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
fn point(p: &Plastic, kelvin: f64) -> QueryPoint {
    at(context(p, 1.0, true), "temperature", kelvin)
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
fn mean(values: &[f64]) -> f64 {
    values.iter().sum::<f64>() / values.len() as f64
}
fn sd(values: &[f64]) -> f64 {
    let average = mean(values);
    (values.iter().map(|v| (v - average).powi(2)).sum::<f64>() / (values.len() - 1) as f64).sqrt()
}
fn sample(
    card: &MaterialCard,
    name: &str,
    dims: Dims,
    point: &QueryPoint,
    expected: f64,
) -> MaterialAnswer {
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
    assert!(answer.receipt.observation_backed);
    assert!(!answer.receipt.source_hashes.is_empty());
    card.claims().verify_receipt(&answer.receipt).unwrap();
    answer
}
fn refuses(card: &MaterialCard, name: &str, dims: Dims, point: &QueryPoint) {
    assert!(
        card.claims()
            .query_typed(&key(name, dims), point, SelectionPolicy::SingleClaimOnly)
            .is_err(),
        "must refuse {name} at {point:?}"
    );
}

#[test]
fn g1_g3_all_native_and_published_values_preserve_reduction_receipts_and_units() {
    for p in &PLASTICS {
        let pack = compile(p);
        let card = pack.card();
        let density = mean(&p.densities);
        let answer = sample(
            card,
            "reference-normalization-density",
            RHO,
            &context(p, 1.0, true),
            density,
        );
        assert!(
            !card
                .claims()
                .claim(answer.receipt.selected)
                .unwrap()
                .validity
                .bounds()
                .contains_key("temperature")
        );
        if p.id == "petg" {
            assert_eq!(
                density, 1314.0,
                "the primary CSV reduction differs from the 1315 headline"
            );
        }
        for j in 0..4 {
            let t = 283.15 + 10.0 * j as f64;
            let values = p.capacity.map(|row| row[j]);
            let (volume, scatter) = (mean(&values), sd(&values));
            sample(card, "volumetric-heat-capacity", C, &point(p, t), volume);
            sample(card, "specific-heat-capacity", CP, &point(p, t), p.cp[j]);
            let statistic = at(context(p, 0.0, true), "temperature", t);
            sample(
                card,
                "computed-volumetric-heat-capacity-sample-standard-deviation",
                C,
                &statistic,
                scatter,
            );
            sample(
                card,
                "reported-specific-heat-capacity-standard-deviation",
                CP,
                &statistic,
                p.cp_sd[j],
            );
            assert_eq!((volume / density).round_ties_even(), p.cp[j]);
            assert_eq!((scatter / density).round_ties_even(), p.cp_sd[j]);
            assert_ne!(
                volume,
                density * p.cp[j],
                "published rounded Cp must remain rounded"
            );
        }
        for j in 0..2 {
            let t = 288.15 + 30.0 * j as f64;
            sample(card, "thermal-conductivity", K, &point(p, t), p.k[j]);
            sample(
                card,
                "apparent-thermal-conductivity",
                K,
                &at(context(p, 1.0, false), "test-mean-temperature", t),
                p.k[j],
            );
            sample(
                card,
                "reported-apparent-thermal-conductivity-standard-deviation",
                K,
                &at(context(p, 0.0, false), "test-mean-temperature", t),
                p.k_sd[j],
            );
            assert_eq!(
                (mean(&p.k_readings[j]) * 1000.0).round_ties_even() / 1000.0,
                p.k[j]
            );
            assert_eq!(
                (sd(&p.k_readings[j]) * 1000.0).round_ties_even() / 1000.0,
                p.k_sd[j]
            );
        }
    }
}

#[test]
fn g0_source_conditions_domains_and_scatter_cannot_be_silently_reinterpreted() {
    for p in &PLASTICS {
        let pack = compile(p);
        let card = pack.card();
        for (axis, value) in [
            ("source-unconditioned", 0.0),
            ("source-grade-known", 1.0),
            ("source-formulation-known", 1.0),
            ("source-moisture-content-known", 1.0),
            ("source-contact-pressure-known", 1.0),
            ("source-specimen-count", 9.0),
            ("source-stock-thickness", p.stock * 1.1),
        ] {
            refuses(
                card,
                "volumetric-heat-capacity",
                C,
                &point(p, 298.15).with(axis, value).unwrap(),
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
                &point(p, 298.15).with(axis, value).unwrap(),
            );
        }
        refuses(
            card,
            "thermal-conductivity",
            K,
            &at(context(p, 1.0, false), "temperature", 298.15),
        );
        refuses(card, "apparent-thermal-conductivity", K, &point(p, 298.15));
        for t in [283.14, 313.16] {
            refuses(card, "volumetric-heat-capacity", C, &point(p, t));
        }
        for t in [288.14, 318.16] {
            refuses(card, "thermal-conductivity", K, &point(p, t));
        }
        refuses(card, "volumetric-heat-capacity", CP, &point(p, 298.15));
        let wrong_temperature = point(p, 298.15)
            .with_quantity(
                "temperature",
                semantic(QuantityKind::TemperatureDifference),
                298.15,
            )
            .unwrap();
        refuses(card, "volumetric-heat-capacity", C, &wrong_temperature);
        let wrong_delta = point(p, 298.15)
            .with_quantity(
                "source-hfm-temperature-difference",
                semantic(QuantityKind::AbsoluteTemperature),
                20.0,
            )
            .unwrap();
        refuses(card, "thermal-conductivity", K, &wrong_delta);
        for name in ["density", "thermal-expansion-coefficient", "young-modulus"] {
            assert!(card.claims().claims_for(name).is_empty());
        }
        let sd_key = key(
            "computed-volumetric-heat-capacity-sample-standard-deviation",
            C,
        );
        let lo = at(context(p, 0.0, true), "temperature", 283.15);
        let hi = at(context(p, 0.0, true), "temperature", 313.15);
        let selected = card
            .claims()
            .query_typed(&sd_key, &lo, SelectionPolicy::SingleClaimOnly)
            .unwrap();
        refuses(
            card,
            sd_key.name(),
            C,
            &at(context(p, 0.0, true), "temperature", 298.15),
        );
        refuses(card, sd_key.name(), C, &point(p, 283.15));
        assert!(
            card.claims()
                .query_line_typed(
                    &sd_key,
                    "temperature",
                    &lo,
                    &hi,
                    ClaimSelection::Pinned(selected.receipt.selected)
                )
                .is_err()
        );
        refuses(
            card,
            "reference-normalization-density",
            RHO,
            &context(p, 1.0, true)
                .with("source-density-temperature-known", 1.0)
                .unwrap(),
        );
        // A numerical normalization input without any measured density temperature
        // must not become a continuous rho(T) input to the mass-specific adapter.
        assert!(
            PreparedThermalTransport::try_new(
                card,
                ThermalTransportSourceQuery {
                    density: &key("reference-normalization-density", RHO),
                    specific_heat: &key("specific-heat-capacity", CP),
                    conductivity: &key("thermal-conductivity", K),
                    temperature_axis: "temperature",
                    lower: &point(p, 288.15),
                    upper: &point(p, 313.15),
                    selections: [POLICY; 3],
                }
            )
            .is_err()
        );
    }
}

#[test]
fn g1_native_volume_curves_drive_prepared_transport_and_heat_storage_without_density() {
    for p in &PLASTICS {
        let pack = compile(p);
        let card = pack.card();
        let capacity = key("volumetric-heat-capacity", C);
        let conductivity = key("thermal-conductivity", K);
        let lower = point(p, 288.15);
        let upper = point(p, 313.15);
        let q = VolumetricThermalPathQuery {
            volumetric_heat_capacity: &capacity,
            conductivity: &conductivity,
            temperature_axis: "temperature",
            lower: &lower,
            upper: &upper,
            selections: [POLICY; 2],
        };
        let prepared = PreparedVolumetricThermalTransport::try_new(card, q).unwrap();
        assert_eq!(prepared.temperature_span_k(), [288.15, 313.15]);
        let c: [f64; 4] = std::array::from_fn(|j| mean(&p.capacity.map(|row| row[j])));
        let at25 = prepared.at_temperature(298.15).unwrap();
        let c25 = (c[1] + c[2]) / 2.0;
        let k25 = p.k[0] + (p.k[1] - p.k[0]) / 3.0;
        close(at25.volumetric_heat_capacity_j_m3_k(), c25);
        close(at25.conductivity_w_m_k(), k25);
        close(at25.diffusivity_m2_s(), k25 / c25);
        close(at25.effusivity_j_m2_k_sqrt_s(), (k25 * c25).sqrt());
        let c15 = (c[0] + c[1]) / 2.0;
        let energy =
            5.0 * (c15 + c[1]) / 2.0 + 10.0 * (c[1] + c[2]) / 2.0 + 10.0 * (c[2] + c[3]) / 2.0;
        let stored = prepared.storage_change_j_m3(288.15, 313.15).unwrap();
        close(stored, energy);
        close(
            prepared.storage_change_j_m3(313.15, 288.15).unwrap(),
            -energy,
        );
        close(
            prepared
                .temperature_after_storage_change(288.15, stored)
                .unwrap(),
            313.15,
        );
        close(
            prepared
                .temperature_after_storage_change(313.15, -stored)
                .unwrap(),
            288.15,
        );
        assert!(
            PreparedVolumetricThermalTransport::try_new(
                card,
                VolumetricThermalPathQuery {
                    lower: &point(p, 283.15),
                    ..q
                }
            )
            .is_err()
        );
        assert!(
            PreparedVolumetricThermalTransport::try_new(
                card,
                VolumetricThermalPathQuery {
                    upper: &point(p, 318.15),
                    ..q
                }
            )
            .is_err()
        );
        assert!(
            PreparedVolumetricThermalTransport::try_new(
                card,
                VolumetricThermalPathQuery {
                    conductivity: &key("apparent-thermal-conductivity", K),
                    ..q
                }
            )
            .is_err()
        );
        assert!(
            prepared
                .temperature_after_storage_change(288.15, 1.01 * energy)
                .is_err()
        );
    }
}
