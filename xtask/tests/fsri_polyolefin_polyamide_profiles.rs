//! PP, LDPE, PET and nylon use the primary FSRI HFM reductions and bounded thermal consumers.
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
        id: "pp",
        stock: 0.0033,
        densities: [884.0, 892.0, 884.0],
        capacity: [
            [1490701.0, 1555453.0, 1631050.0, 1696007.0],
            [1507089.0, 1570552.0, 1652433.0, 1702766.0],
            [1497365.0, 1565518.0, 1639811.0, 1708719.0],
        ],
        cp: [1690.0, 1764.0, 1851.0, 1920.0],
        cp_sd: [9.0, 9.0, 12.0, 7.0],
        k_readings: [
            [
                0.1368, 0.1420, 0.1424, 0.1386, 0.1416, 0.1421, 0.1293, 0.1320, 0.1322,
            ],
            [
                0.1544, 0.1549, 0.1551, 0.1538, 0.1541, 0.1550, 0.1450, 0.1454, 0.1450,
            ],
        ],
        k: [0.137, 0.151],
        k_sd: [0.005, 0.005],
    },
    Plastic {
        id: "ldpe",
        stock: 0.0033,
        densities: [925.0, 930.0, 927.0],
        capacity: [
            [1911039.0, 1952408.0, 2170804.0, 2346103.0],
            [1912463.0, 1946004.0, 2181253.0, 2421462.0],
            [1921474.0, 1995072.0, 2205595.0, 2464944.0],
        ],
        cp: [2065.0, 2118.0, 2357.0, 2600.0],
        cp_sd: [6.0, 29.0, 19.0, 65.0],
        k_readings: [
            [
                0.1895, 0.1950, 0.1963, 0.1834, 0.1858, 0.1881, 0.2120, 0.2201, 0.2249,
            ],
            [
                0.2106, 0.2134, 0.2136, 0.2032, 0.2057, 0.2060, 0.2332, 0.2336, 0.2367,
            ],
        ],
        k: [0.199, 0.217],
        k_sd: [0.016, 0.013],
    },
    Plastic {
        id: "pet",
        stock: 0.0069,
        densities: [1402.0, 1401.0, 1400.0],
        capacity: [
            [1486121.0, 1524330.0, 1589852.0, 1634221.0],
            [1501472.0, 1576477.0, 1593880.0, 1648320.0],
            [1495894.0, 1541567.0, 1593088.0, 1637435.0],
        ],
        cp: [1067.0, 1105.0, 1137.0, 1171.0],
        cp_sd: [6.0, 19.0, 2.0, 5.0],
        k_readings: [
            [
                0.1813, 0.1832, 0.1816, 0.1818, 0.1861, 0.1860, 0.1738, 0.1763, 0.1765,
            ],
            [
                0.1959, 0.1957, 0.1935, 0.1981, 0.1978, 0.1978, 0.1882, 0.1882, 0.1883,
            ],
        ],
        k: [0.181, 0.194],
        k_sd: [0.004, 0.004],
    },
    Plastic {
        id: "nylon",
        stock: 0.0035,
        densities: [1065.0, 1074.0, 1072.0],
        capacity: [
            [1596854.0, 1702472.0, 1831220.0, 1917091.0],
            [1644837.0, 1767301.0, 1916064.0, 1991341.0],
            [1636874.0, 1741989.0, 1887477.0, 1986433.0],
        ],
        cp: [1519.0, 1623.0, 1755.0, 1836.0],
        cp_sd: [24.0, 31.0, 40.0, 39.0],
        k_readings: [
            [
                0.1292, 0.1329, 0.1323, 0.1413, 0.1461, 0.1472, 0.1478, 0.1507, 0.1503,
            ],
            [
                0.1442, 0.1438, 0.1449, 0.1585, 0.1576, 0.1580, 0.1613, 0.1610, 0.1605,
            ],
        ],
        k: [0.142, 0.154],
        k_sd: [0.008, 0.008],
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
        "fsri-polyolefin-polyamide-test-{}-{nonce}-{}.fsmatpk",
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
            "source-metadata-stock-thickness",
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
        if let Some(headline) = match p.id {
            "pp" => Some(886.0),
            "ldpe" => Some(927.0),
            "nylon" => Some(1070.0),
            _ => None,
        } {
            assert_ne!(
                density, headline,
                "published Cp uses the primary CSV mean, not the rounded headline"
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
            ("source-metadata-stock-thickness", p.stock * 1.1),
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
