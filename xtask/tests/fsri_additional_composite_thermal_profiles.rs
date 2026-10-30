//! Additional FSRI product populations retain conditioning and actual run support.
#![deny(unsafe_code)]

use fs_matdb::{
    ClaimSelection, MaterialCard, MaterialStateId, NormalizedMaterialCardPack, NormalizedPack,
    PropertyKey, PropertyValue, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::thermal::{
    conduction::{ConductivityCurve, LayerConductionConfig, TwoLayerConduction},
    transport::{
        storage::{HeatStorageAmount, HeatStorageComponentQuery, PreparedHeatStorage},
        volumetric::{PreparedVolumetricThermalTransport, VolumetricThermalPathQuery},
    },
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
const SINGLE: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);

struct Source {
    slug: &'static str,
    family: &'static str,
    dried: f64,
    stock_m: f64,
    capacity: &'static [[f64; 4]],
    density: [f64; 3],
    conductivity: [[f64; 9]; 2],
    cp: [f64; 4],
    cp_sd: Option<[f64; 4]>,
    k: [f64; 2],
    k_sd: [f64; 2],
}

// Independent literal numerical facts from the pinned instrument final SI
// Results Tables and density CSVs. Nine conductivity cycles are not nine
// specimens. Memory foam has one capacity run; its source SD is undefined.
const SOURCES: &[Source] = &[
    Source {
        slug: "fsri-hips-unconditioned-2021",
        family: "source-fsri-hips",
        dried: 0.0,
        stock_m: 0.0029,
        capacity: &[
            [1322563.0, 1361253.0, 1412943.0, 1459136.0],
            [1314347.0, 1358218.0, 1424548.0, 1443467.0],
            [1299537.0, 1346727.0, 1390234.0, 1421387.0],
        ],
        density: [1068.0, 1069.0, 1063.0],
        conductivity: [
            [
                0.1152, 0.1159, 0.1164, 0.1173, 0.1194, 0.1193, 0.1119, 0.1160, 0.1163,
            ],
            [
                0.1275, 0.1269, 0.1264, 0.1304, 0.1303, 0.1300, 0.1258, 0.1265, 0.1267,
            ],
        ],
        cp: [1230.0, 1271.0, 1321.0, 1351.0],
        cp_sd: Some([11.0, 7.0, 16.0, 18.0]),
        k: [0.116, 0.128],
        k_sd: [0.002, 0.002],
    },
    Source {
        slug: "fsri-memory-foam-pad-dried-2022",
        family: "source-fsri-memory-foam-pad",
        dried: 1.0,
        stock_m: 0.0128,
        capacity: &[[219352.0, 225888.0, 230918.0, 232225.0]],
        density: [136.0, 136.0, 130.0],
        conductivity: [
            [
                0.03827, 0.03829, 0.03827, 0.03827, 0.03831, 0.03833, 0.03764, 0.03768, 0.03769,
            ],
            [
                0.04171, 0.04174, 0.04175, 0.04175, 0.04179, 0.04181, 0.04121, 0.04126, 0.04126,
            ],
        ],
        cp: [1637.0, 1686.0, 1723.0, 1733.0],
        cp_sd: None,
        k: [0.038, 0.042],
        k_sd: [0.0, 0.0],
    },
    Source {
        slug: "fsri-homasote-unconditioned-2021",
        family: "source-fsri-homasote",
        dried: 0.0,
        stock_m: 0.0116,
        capacity: &[
            [594245.0, 621600.0, 653914.0, 687676.0],
            [597230.0, 626270.0, 660505.0, 693939.0],
            [594741.0, 629484.0, 657661.0, 697812.0],
        ],
        density: [450.0, 450.0, 445.0],
        conductivity: [
            [
                0.07034, 0.06941, 0.06919, 0.07042, 0.06982, 0.06968, 0.07001, 0.06942, 0.06927,
            ],
            [
                0.07803, 0.07590, 0.07506, 0.07827, 0.07630, 0.07553, 0.07819, 0.07593, 0.07506,
            ],
        ],
        cp: [1328.0, 1396.0, 1466.0, 1546.0],
        cp_sd: Some([4.0, 9.0, 7.0, 11.0]),
        k: [0.070, 0.076],
        k_sd: [0.0, 0.001],
    },
    Source {
        slug: "fsri-homasote-dried-2021",
        family: "source-fsri-homasote",
        dried: 1.0,
        stock_m: 0.0116,
        capacity: &[
            [495544.0, 515142.0, 538663.0, 560211.0],
            [476095.0, 493134.0, 513719.0, 535155.0],
            [489061.0, 500812.0, 523034.0, 544583.0],
        ],
        density: [413.0, 409.0, 406.0],
        conductivity: [
            [
                0.06085, 0.06095, 0.06096, 0.05960, 0.05968, 0.05977, 0.06108, 0.06109, 0.06109,
            ],
            [
                0.06438, 0.06436, 0.06435, 0.06318, 0.06317, 0.06317, 0.06477, 0.06479, 0.06462,
            ],
        ],
        cp: [1189.0, 1229.0, 1283.0, 1335.0],
        cp_sd: Some([24.0, 27.0, 31.0, 31.0]),
        k: [0.061, 0.064],
        k_sd: [0.001, 0.001],
    },
    Source {
        slug: "fsri-engineered-flooring-unconditioned-2021",
        family: "source-fsri-engineered-flooring",
        dried: 0.0,
        stock_m: 0.0032,
        capacity: &[
            [1111808.0, 1170658.0, 1234319.0, 1306565.0],
            [1126276.0, 1195393.0, 1245476.0, 1316210.0],
            [1109921.0, 1166206.0, 1230879.0, 1299542.0],
        ],
        density: [812.0, 816.0, 817.0],
        conductivity: [
            [
                0.1266, 0.1275, 0.1265, 0.1222, 0.1224, 0.1219, 0.1250, 0.1248, 0.1239,
            ],
            [
                0.1366, 0.1356, 0.1347, 0.1321, 0.1318, 0.1308, 0.1348, 0.1333, 0.1324,
            ],
        ],
        cp: [1369.0, 1445.0, 1518.0, 1604.0],
        cp_sd: Some([11.0, 19.0, 9.0, 10.0]),
        k: [0.125, 0.134],
        k_sd: [0.002, 0.002],
    },
    Source {
        slug: "fsri-engineered-flooring-dried-2021",
        family: "source-fsri-engineered-flooring",
        dried: 1.0,
        stock_m: 0.0032,
        capacity: &[
            [828832.0, 865033.0, 897579.0, 945371.0],
            [808019.0, 840191.0, 881287.0, 922878.0],
            [829503.0, 855184.0, 895872.0, 943048.0],
        ],
        density: [758.0, 761.0, 762.0],
        conductivity: [
            [
                0.08780, 0.08813, 0.08804, 0.08276, 0.08319, 0.08331, 0.08721, 0.08757, 0.08764,
            ],
            [
                0.09452, 0.09447, 0.09447, 0.08978, 0.08977, 0.08988, 0.09393, 0.09396, 0.09404,
            ],
        ],
        cp: [1081.0, 1122.0, 1173.0, 1232.0],
        cp_sd: Some([16.0, 16.0, 12.0, 16.0]),
        k: [0.086, 0.093],
        k_sd: [0.002, 0.002],
    },
];

fn compiled(source: &Source) -> NormalizedMaterialCardPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let output = std::env::temp_dir().join(format!(
        "fsri-additional-{}-{nonce}-{}.fsmatpk",
        std::process::id(),
        source.slug
    ));
    let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(root.join(format!("data/matdb/seed-v1/{}/manifest.tsv", source.slug)))
        .arg("--out")
        .arg(&output)
        .current_dir(root)
        .env("CARGO_WORKSPACE_DIR", root)
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "{}: {}",
        source.slug,
        String::from_utf8_lossy(&run.stderr)
    );
    let bytes = std::fs::read(output).unwrap();
    let raw = NormalizedPack::from_bytes(&bytes).unwrap();
    let raw = NormalizedPack::from_bytes_verified(raw.content_hash(), &bytes).unwrap();
    assert_eq!(raw.to_bytes(), bytes);
    assert_eq!(raw.pack_id(), source.slug);
    assert_eq!(
        raw.claims().claims_ordered().count(),
        if source.cp_sd.is_some() { 8 } else { 6 }
    );
    let values: usize = raw
        .claims()
        .claims_ordered()
        .map(|(_, c)| match &c.value {
            PropertyValue::Scalar { .. } => 1,
            PropertyValue::Curve { knots, .. } => knots.len(),
            _ => panic!("unexpected payload"),
        })
        .sum();
    assert_eq!(values, if source.cp_sd.is_some() { 23 } else { 15 });
    let pack = NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: source.slug.into(),
            phase: "solid".into(),
            process: "FSRI source population; conditioning and separate run counts retained".into(),
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
fn context(source: &Source, statistic: f64, local: bool) -> QueryPoint {
    let mut p = QueryPoint::new();
    for (axis, value) in [
        (source.family, 1.0),
        ("source-unconditioned", 1.0 - source.dried),
        ("source-dried", source.dried),
        ("source-grade-known", 0.0),
        ("source-formulation-known", 0.0),
        ("source-moisture-content-known", 0.0),
        ("source-contact-pressure-known", 0.0),
        ("source-material-axis-known", 0.0),
        ("source-filler-fraction-known", 0.0),
        ("source-layer-construction-known", 0.0),
        ("source-specimens-matched", 0.0),
        ("source-stock-coupon-equivalence-known", 0.0),
        ("source-unique-specimen-count-confirmed", 0.0),
        ("source-header-labels-reconciled", 0.0),
        ("source-statistic-is-mean", statistic),
        ("source-capacity-run-count", source.capacity.len() as f64),
        ("source-conductivity-run-count", 3.0),
        ("source-density-run-count", 3.0),
        ("source-hfm-through-thickness", 1.0),
        ("source-hfm-upward-heat-flow", 1.0),
        ("source-density-temperature-known", 0.0),
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
    p.with_quantity(
        "source-stock-thickness",
        QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0])),
        source.stock_m,
    )
    .unwrap()
    .with_quantity(
        "source-hfm-temperature-difference",
        semantic(QuantityKind::TemperatureDifference),
        20.0,
    )
    .unwrap()
}
fn at(p: QueryPoint, axis: &str, t: f64) -> QueryPoint {
    p.with_quantity(axis, semantic(QuantityKind::AbsoluteTemperature), t)
        .unwrap()
}
fn point(source: &Source, t: f64) -> QueryPoint {
    at(context(source, 1.0, true), "temperature", t)
}
fn key(name: &str, dims: Dims) -> PropertyKey {
    if name == "thermal-conductivity" || name == "apparent-thermal-conductivity" {
        PropertyKey::with_quantity(name, semantic(QuantityKind::ThermalConductivity))
    } else {
        PropertyKey::new(name, dims)
    }
}
fn mean(v: &[f64]) -> f64 {
    v.iter().sum::<f64>() / v.len() as f64
}
fn sd(v: &[f64]) -> f64 {
    let m = mean(v);
    (v.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (v.len() - 1) as f64).sqrt()
}
fn native_capacity(source: &Source, j: usize) -> f64 {
    mean(&source.capacity.iter().map(|r| r[j]).collect::<Vec<_>>())
}
fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 3e-11 * b.abs().max(1e-8),
        "{a:.16e} != {b:.16e}"
    );
}
fn check(card: &MaterialCard, name: &str, dims: Dims, p: &QueryPoint, expected: f64) {
    let a = card
        .claims()
        .query_typed(&key(name, dims), p, SelectionPolicy::SingleClaimOnly)
        .unwrap();
    close(a.evidence.value.value, expected);
    assert_eq!(a.evidence.value.dims, dims);
    assert_eq!(a.evidence.value.uncertainty, UncertaintyModel::Unstated);
    assert!(a.receipt.observation_backed && !a.receipt.source_hashes.is_empty());
    card.claims().verify_receipt(&a.receipt).unwrap();
}
fn refuses(card: &MaterialCard, name: &str, dims: Dims, p: &QueryPoint) {
    assert!(
        card.claims()
            .query_typed(&key(name, dims), p, SelectionPolicy::SingleClaimOnly)
            .is_err(),
        "accepted {name}"
    );
}

#[test]
fn g1_g3_primary_run_reductions_native_units_and_evidence_survive_compilation() {
    for source in SOURCES {
        let pack = compiled(source);
        let card = pack.card();
        let rho = mean(&source.density);
        check(
            card,
            "reference-normalization-density",
            RHO,
            &context(source, 1.0, false),
            rho,
        );
        for j in 0..4 {
            let t = 283.15 + 10.0 * j as f64;
            let observations: Vec<_> = source.capacity.iter().map(|row| row[j]).collect();
            let c = mean(&observations);
            check(card, "volumetric-heat-capacity", C, &point(source, t), c);
            assert_eq!((c / rho).round(), source.cp[j]);
            check(
                card,
                "specific-heat-capacity",
                CP,
                &point(source, t),
                source.cp[j],
            );
            if let Some(cp_sd) = source.cp_sd {
                let scatter = sd(&observations);
                check(
                    card,
                    "computed-volumetric-heat-capacity-sample-standard-deviation",
                    C,
                    &at(context(source, 0.0, false), "temperature", t),
                    scatter,
                );
                assert_eq!((scatter / rho).round(), cp_sd[j]);
                check(
                    card,
                    "reported-specific-heat-capacity-standard-deviation",
                    CP,
                    &at(context(source, 0.0, false), "temperature", t),
                    cp_sd[j],
                );
            }
        }
        for j in 0..2 {
            let t = 288.15 + 30.0 * j as f64;
            assert_eq!(
                (mean(&source.conductivity[j]) * 1000.0).round() / 1000.0,
                source.k[j]
            );
            assert_eq!(
                (sd(&source.conductivity[j]) * 1000.0).round() / 1000.0,
                source.k_sd[j]
            );
            check(
                card,
                "apparent-thermal-conductivity",
                K,
                &at(context(source, 1.0, false), "test-mean-temperature", t),
                source.k[j],
            );
            check(
                card,
                "thermal-conductivity",
                K,
                &point(source, t),
                source.k[j],
            );
            check(
                card,
                "reported-apparent-thermal-conductivity-standard-deviation",
                K,
                &at(context(source, 0.0, false), "test-mean-temperature", t),
                source.k_sd[j],
            );
        }
        // Independent interpolation of the retained source reduction, with no
        // rewriting rounded Cp as an unrounded C/rho value.
        check(
            card,
            "volumetric-heat-capacity",
            C,
            &point(source, 298.15),
            0.5 * (native_capacity(source, 1) + native_capacity(source, 2)),
        );
        check(
            card,
            "specific-heat-capacity",
            CP,
            &point(source, 298.15),
            0.5 * (source.cp[1] + source.cp[2]),
        );
        check(
            card,
            "thermal-conductivity",
            K,
            &point(source, 298.15),
            source.k[0] + (source.k[1] - source.k[0]) / 3.0,
        );
    }
}

#[test]
fn g0_conditioning_single_run_scatter_labels_and_mean_temperature_cannot_be_erased() {
    for source in SOURCES {
        let pack = compiled(source);
        let card = pack.card();
        for (axis, value) in [
            (source.family, 0.0),
            ("source-dried", 1.0 - source.dried),
            ("source-unconditioned", source.dried),
            ("source-grade-known", 1.0),
            ("source-formulation-known", 1.0),
            ("source-moisture-content-known", 1.0),
            ("source-contact-pressure-known", 1.0),
            ("source-material-axis-known", 1.0),
            ("source-filler-fraction-known", 1.0),
            ("source-layer-construction-known", 1.0),
            ("source-specimens-matched", 1.0),
            ("source-stock-coupon-equivalence-known", 1.0),
            ("source-unique-specimen-count-confirmed", 1.0),
            ("source-header-labels-reconciled", 1.0),
            ("source-stock-thickness", source.stock_m * 2.0),
            ("source-statistic-is-mean", 0.0),
            ("source-capacity-run-count", 9.0),
        ] {
            refuses(
                card,
                "volumetric-heat-capacity",
                C,
                &point(source, 298.15).with(axis, value).unwrap(),
            );
        }
        refuses(
            card,
            "volumetric-heat-capacity",
            C,
            &at(QueryPoint::new(), "temperature", 298.15),
        );
        refuses(card, "volumetric-heat-capacity", CP, &point(source, 298.15));
        refuses(
            card,
            "volumetric-heat-capacity",
            C,
            &point(source, 298.15)
                .with_quantity(
                    "temperature",
                    semantic(QuantityKind::TemperatureDifference),
                    298.15,
                )
                .unwrap(),
        );
        for t in [283.15_f64.next_down(), 313.15_f64.next_up()] {
            refuses(card, "volumetric-heat-capacity", C, &point(source, t));
            refuses(card, "specific-heat-capacity", CP, &point(source, t));
        }
        for t in [288.15_f64.next_down(), 318.15_f64.next_up()] {
            refuses(card, "thermal-conductivity", K, &point(source, t));
        }
        for (axis, value) in [
            ("mean-test-temperature-as-local-model", 0.0),
            ("source-hfm-through-thickness", 0.0),
            ("source-hfm-upward-heat-flow", 0.0),
            ("source-hfm-temperature-difference", 25.0),
            ("source-conductivity-run-count", 9.0),
        ] {
            refuses(
                card,
                "thermal-conductivity",
                K,
                &point(source, 298.15).with(axis, value).unwrap(),
            );
        }
        refuses(
            card,
            "thermal-conductivity",
            K,
            &at(context(source, 1.0, false), "temperature", 298.15),
        );
        refuses(
            card,
            "apparent-thermal-conductivity",
            K,
            &point(source, 298.15),
        );
        refuses(
            card,
            "reference-normalization-density",
            RHO,
            &context(source, 1.0, false)
                .with("source-density-temperature-known", 1.0)
                .unwrap(),
        );
        refuses(
            card,
            "specific-heat-capacity",
            CP,
            &point(source, 298.15)
                .with("source-density-run-count", 1.0)
                .unwrap(),
        );
        for name in [
            "computed-volumetric-heat-capacity-sample-standard-deviation",
            "reported-specific-heat-capacity-standard-deviation",
        ] {
            refuses(
                card,
                name,
                if name.starts_with("computed") { C } else { CP },
                &at(context(source, 0.0, false), "temperature", 298.15),
            );
            if source.capacity.len() == 1 {
                assert!(card.claims().claims_for(name).is_empty());
            }
        }
        for name in [
            "density",
            "young-modulus",
            "linear-thermal-expansion-coefficient",
            "moisture-dependent-conductivity",
        ] {
            assert!(card.claims().claims_for(name).is_empty());
        }
    }
    // Same source family does not authorize crossing a conditioning boundary.
    let dry = compiled(&SOURCES[3]);
    refuses(
        dry.card(),
        "volumetric-heat-capacity",
        C,
        &point(&SOURCES[2], 298.15),
    );
    let single = compiled(&SOURCES[1]);
    refuses(
        single.card(),
        "volumetric-heat-capacity",
        C,
        &point(&SOURCES[1], 298.15)
            .with("source-capacity-run-count", 3.0)
            .unwrap(),
    );
}

#[test]
fn g1_g3_native_storage_and_integrated_flux_distinguish_conditioned_products() {
    let config = LayerConductionConfig {
        max_iterations: 128,
        absolute_heat_tolerance_w: 1e-12,
        relative_heat_tolerance: 1e-12,
    };
    let mut heated = Vec::new();
    let mut powers = Vec::new();
    let mut identities = Vec::new();
    for source in SOURCES {
        let pack = compiled(source);
        let card = pack.card();
        let capacity = key("volumetric-heat-capacity", C);
        let conductivity = key("thermal-conductivity", K);
        let lower = point(source, 288.15);
        let upper = point(source, 313.15);
        let query = VolumetricThermalPathQuery {
            volumetric_heat_capacity: &capacity,
            conductivity: &conductivity,
            temperature_axis: "temperature",
            lower: &lower,
            upper: &upper,
            selections: [SINGLE; 2],
        };
        // This explicitly authored local scalar/sensible model adopts two
        // source populations; it does not assert matched coupons or isotropy.
        let path = PreparedVolumetricThermalTransport::try_new(card, query).unwrap();
        let c20 = native_capacity(source, 1);
        let c30 = native_capacity(source, 2);
        let slope = (c30 - c20) / 10.0;
        let k25 = source.k[0] + (source.k[1] - source.k[0]) / 3.0;
        let sample = path.at_temperature(298.15).unwrap();
        close(sample.volumetric_heat_capacity_j_m3_k(), 0.5 * (c20 + c30));
        close(sample.conductivity_w_m_k(), k25);
        close(sample.diffusivity_m2_s(), k25 / (0.5 * (c20 + c30)));
        close(
            path.storage_change_j_m3(293.15, 303.15).unwrap(),
            5.0 * (c20 + c30),
        );
        let supplied = 300000.0;
        let delta = 2.0 * supplied / (c20 + (c20 * c20 + 2.0 * slope * supplied).sqrt());
        let final_temperature = path
            .temperature_after_storage_change(293.15, supplied)
            .unwrap();
        close(final_temperature, 293.15 + delta);
        heated.push(final_temperature);
        let entropy = slope * 10.0 + (c20 - slope * 293.15) * (303.15_f64 / 293.15).ln();
        close(
            path.sensible_entropy_change_j_m3_k(293.15, 303.15).unwrap(),
            entropy,
        );

        let curve = ConductivityCurve::new(
            path.source_lines()[1]
                .knots
                .iter()
                .map(|(t, a)| [*t, a.evidence.value.value])
                .collect(),
        )
        .unwrap();
        let layers = TwoLayerConduction::new([&curve; 2], [0.01; 2], 0.01, 0.0).unwrap();
        let flux = layers
            .evaluate([303.15, 293.15], config, || Ok(()))
            .unwrap();
        close(flux.power_a_to_b_w, 0.01 / 0.02 * 10.0 * k25);
        close(
            layers
                .evaluate([293.15, 303.15], config, || Ok(()))
                .unwrap()
                .power_a_to_b_w,
            -flux.power_a_to_b_w,
        );
        let contact = TwoLayerConduction::new([&curve; 2], [0.01; 2], 0.01, 0.1)
            .unwrap()
            .evaluate([303.15, 293.15], config, || Ok(()))
            .unwrap();
        assert!(contact.power_a_to_b_w < flux.power_a_to_b_w && contact.power_a_to_b_w > 0.0);
        powers.push(flux.power_a_to_b_w);

        // Storage can also use the native four-knot capacity span without k,
        // supplied density, or an assertion about the nominal stock geometry.
        let full_lower = point(source, 283.15);
        let full_upper = point(source, 313.15);
        let storage = PreparedHeatStorage::try_new(&[HeatStorageComponentQuery {
            id: 1,
            card,
            capacity: &capacity,
            amount: HeatStorageAmount::ReferenceVolumeM3(0.002),
            temperature_axis: "temperature",
            lower: &full_lower,
            upper: &full_upper,
            selection: SINGLE,
        }])
        .unwrap();
        let expected = 0.002
            * (0..3)
                .map(|j| 5.0 * (native_capacity(source, j) + native_capacity(source, j + 1)))
                .sum::<f64>();
        let energy = storage.energy_change_j(283.15, 313.15).unwrap();
        close(energy, expected);
        close(
            storage
                .temperature_after_energy_change(283.15, energy)
                .unwrap(),
            313.15,
        );
        close(storage.energy_change_j(313.15, 283.15).unwrap(), -energy);
        assert!(
            storage
                .temperature_after_energy_change(293.15, 1e9)
                .is_err()
        );
        assert!(
            PreparedVolumetricThermalTransport::try_new(
                card,
                VolumetricThermalPathQuery {
                    lower: &full_lower,
                    ..query
                }
            )
            .is_err()
        );
        assert!(
            PreparedVolumetricThermalTransport::try_new(
                card,
                VolumetricThermalPathQuery {
                    upper: &point(source, 318.15),
                    ..query
                }
            )
            .is_err()
        );
        for line in path.source_lines() {
            for (_, answer) in &line.knots {
                card.claims().verify_receipt(&answer.receipt).unwrap();
            }
        }
        identities.push(storage.identity());
    }
    assert!(heated[3] > heated[2] && heated[5] > heated[4]);
    assert!(powers[3] < powers[2] && powers[5] < powers[4]);
    assert!(heated[1] > heated[0] && powers[1] < powers[0]);
    assert_ne!(identities[2], identities[3]);
    assert_ne!(identities[4], identities[5]);
}
