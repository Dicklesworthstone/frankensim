//! Conditioned cement/gypsum boards and lumber retain their measured thermal response.
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

const C: Dims = Dims([-1, 1, -2, -1, 0, 0]);
const CP: Dims = Dims([2, 0, -2, -1, 0, 0]);
const K: Dims = Dims([1, 1, -3, -1, 0, 0]);
const RHO: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const POLICY: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);

struct Source {
    slug: &'static str,
    wet_slug: &'static str,
    family: &'static str,
    stock_m: Option<f64>,
    header_consistent: f64,
    // Independent literal instrument Results Table readings, R1/R2/R3.
    capacity: [[f64; 4]; 3],
    density: [f64; 3],
    // Three cycles from each of three specimens; these are nine readings.
    conductivity: [[f64; 9]; 2],
    cp: [f64; 4],
    cp_sd: [f64; 4],
    k: [f64; 2],
    k_sd: [f64; 2],
}

const SOURCES: [Source; 4] = [
    Source {
        slug: "fsri-fiber-cement-board-dried-2021",
        wet_slug: "fsri-fiber-cement-board-unconditioned-2021",
        family: "source-fsri-fiber-cement-board",
        stock_m: Some(0.0085),
        header_consistent: 1.0,
        capacity: [
            [993844.0, 1022843.0, 1053249.0, 1082535.0],
            [1008281.0, 1043708.0, 1072258.0, 1100260.0],
            [984440.0, 1026570.0, 1056503.0, 1086591.0],
        ],
        density: [1293.0, 1287.0, 1251.0],
        conductivity: [
            [
                0.1266, 0.1273, 0.1269, 0.1396, 0.1403, 0.1402, 0.1301, 0.1313, 0.1311,
            ],
            [
                0.1404, 0.1404, 0.1399, 0.1532, 0.1530, 0.1529, 0.1438, 0.1445, 0.1434,
            ],
        ],
        cp: [780.0, 807.0, 831.0, 853.0],
        cp_sd: [9.0, 9.0, 8.0, 7.0],
        k: [0.133, 0.146],
        k_sd: [0.006, 0.006],
    },
    Source {
        slug: "fsri-gypsum-wallboard-dried-2021",
        wet_slug: "fsri-gypsum-wallboard-unconditioned-2021",
        family: "source-fsri-gypsum-wallboard",
        stock_m: Some(0.0128),
        header_consistent: 1.0,
        capacity: [
            [389321.0, 397662.0, 409629.0, 415285.0],
            [394180.0, 403751.0, 414362.0, 420468.0],
            [398215.0, 415488.0, 417316.0, 432323.0],
        ],
        density: [499.0, 494.0, 495.0],
        conductivity: [
            [
                0.08780, 0.08794, 0.08815, 0.09103, 0.09127, 0.09118, 0.08231, 0.08245, 0.08241,
            ],
            [
                0.09111, 0.09106, 0.09092, 0.09437, 0.09417, 0.09398, 0.08601, 0.08601, 0.08583,
            ],
        ],
        cp: [794.0, 818.0, 834.0, 852.0],
        cp_sd: [9.0, 18.0, 8.0, 18.0],
        k: [0.087, 0.090],
        k_sd: [0.004, 0.004],
    },
    Source {
        slug: "fsri-wood-stud-dried-2021",
        wet_slug: "fsri-wood-stud-unconditioned-2021",
        family: "source-fsri-wood-stud",
        stock_m: None,
        header_consistent: 0.0,
        capacity: [
            [398690.0, 416970.0, 432961.0, 449454.0],
            [410704.0, 426893.0, 446008.0, 460719.0],
            [438628.0, 468616.0, 486449.0, 502146.0],
        ],
        density: [336.0, 346.0, 383.0],
        conductivity: [
            [
                0.08323, 0.08327, 0.08317, 0.08551, 0.08531, 0.08519, 0.09647, 0.09613, 0.09593,
            ],
            [
                0.08776, 0.08760, 0.08766, 0.08992, 0.08983, 0.08972, 0.1011, 0.1009, 0.1006,
            ],
        ],
        cp: [1172.0, 1232.0, 1282.0, 1326.0],
        cp_sd: [58.0, 77.0, 79.0, 78.0],
        k: [0.088, 0.093],
        k_sd: [0.006, 0.006],
    },
    Source {
        slug: "fsri-basswood-panel-dried-2022",
        wet_slug: "fsri-basswood-panel-unconditioned-2021",
        family: "source-fsri-basswood-panel",
        stock_m: Some(0.019),
        header_consistent: 1.0,
        capacity: [
            [448474.0, 464623.0, 482775.0, 498593.0],
            [460440.0, 477574.0, 495618.0, 510579.0],
            [462490.0, 474985.0, 491051.0, 507697.0],
        ],
        density: [372.0, 382.0, 379.0],
        conductivity: [
            [
                0.08535, 0.08547, 0.08547, 0.08500, 0.08509, 0.08512, 0.08633, 0.08638, 0.08637,
            ],
            [
                0.08959, 0.08953, 0.08962, 0.08932, 0.08935, 0.08938, 0.09046, 0.09046, 0.09051,
            ],
        ],
        cp: [1210.0, 1251.0, 1297.0, 1339.0],
        cp_sd: [20.0, 18.0, 17.0, 17.0],
        k: [0.086, 0.090],
        k_sd: [0.001, 0.001],
    },
];

fn compiled(slug: &str) -> NormalizedMaterialCardPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let output = std::env::temp_dir().join(format!(
        "fsri-conditioned-{}-{nonce}.fsmatpk",
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
    let raw = NormalizedPack::from_bytes(&bytes).unwrap();
    let raw = NormalizedPack::from_bytes_verified(raw.content_hash(), &bytes).unwrap();
    assert_eq!(raw.to_bytes(), bytes);
    assert_eq!(raw.claims().claims_ordered().count(), 8);
    let count: usize = raw
        .claims()
        .claims_ordered()
        .map(|(_, c)| match &c.value {
            PropertyValue::Scalar { .. } => 1,
            PropertyValue::Curve { knots, .. } => knots.len(),
            _ => panic!("unexpected source payload"),
        })
        .sum();
    assert_eq!(count, 23);
    NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: slug.into(),
            phase: "solid".into(),
            process: "explicit FSRI product conditioning population".into(),
            revision: 0,
        },
        raw,
    )
    .unwrap()
}
fn semantic(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn context(source: &Source, dried: bool, statistic: f64, local: bool) -> QueryPoint {
    let mut p = QueryPoint::new();
    for (axis, v) in [
        (source.family, 1.0),
        ("source-dried", f64::from(dried)),
        ("source-unconditioned", f64::from(!dried)),
        ("source-grade-known", 0.0),
        ("source-formulation-known", 0.0),
        ("source-moisture-content-known", 0.0),
        ("source-contact-pressure-known", 0.0),
        ("source-material-axis-known", 0.0),
        ("source-specimens-matched", 0.0),
        ("source-specimen-count", 3.0),
        (
            "source-metadata-stock-thickness-known",
            f64::from(source.stock_m.is_some()),
        ),
        ("source-botanical-species-known", 0.0),
        ("source-treatment-known", 0.0),
        ("source-statistic-is-mean", statistic),
        ("source-hfm-through-thickness", 1.0),
        ("source-hfm-upward-heat-flow", 1.0),
        ("source-density-temperature-known", 0.0),
        ("source-density-header-consistent", source.header_consistent),
    ] {
        p = p
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), v)
            .unwrap();
    }
    if let Some(stock) = source.stock_m {
        // Both extant manifest vocabularies identify stock, not test/installed thickness.
        for axis in ["source-metadata-stock-thickness", "source-stock-thickness"] {
            p = p
                .with_quantity(
                    axis,
                    QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0])),
                    stock,
                )
                .unwrap();
        }
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
fn point(s: &Source, dried: bool, t: f64) -> QueryPoint {
    at(context(s, dried, 1.0, true), "temperature", t)
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
fn mean(v: &[f64]) -> f64 {
    v.iter().sum::<f64>() / v.len() as f64
}
fn sd(v: &[f64]) -> f64 {
    let m = mean(v);
    (v.iter().map(|x| (x - m) * (x - m)).sum::<f64>() / (v.len() - 1) as f64).sqrt()
}
fn close(a: f64, b: f64) {
    assert!((a - b).abs() <= 3e-11 * b.abs().max(1e-10), "{a} != {b}");
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
            .is_err()
    );
}
fn prepared(
    card: &MaterialCard,
    source: &Source,
    dried: bool,
) -> PreparedVolumetricThermalTransport {
    PreparedVolumetricThermalTransport::try_new(
        card,
        VolumetricThermalPathQuery {
            volumetric_heat_capacity: &key("volumetric-heat-capacity", C),
            conductivity: &key("thermal-conductivity", K),
            temperature_axis: "temperature",
            lower: &point(source, dried, 288.15),
            upper: &point(source, dried, 313.15),
            selections: [POLICY; 2],
        },
    )
    .unwrap()
}

#[test]
fn g1_g3_conditioned_measurements_reductions_units_and_receipts_survive_compilation() {
    for s in &SOURCES {
        let pack = compiled(s.slug);
        let card = pack.card();
        let rho = mean(&s.density);
        check(
            card,
            "reference-normalization-density",
            RHO,
            &context(s, true, 1.0, false),
            rho,
        );
        for j in 0..4 {
            let v = s.capacity.map(|row| row[j]);
            let t = 283.15 + 10.0 * j as f64;
            close(s.cp[j], (mean(&v) / rho).round());
            close(s.cp_sd[j], (sd(&v) / rho).round());
            check(
                card,
                "volumetric-heat-capacity",
                C,
                &point(s, true, t),
                mean(&v),
            );
            check(
                card,
                "specific-heat-capacity",
                CP,
                &point(s, true, t),
                s.cp[j],
            );
            let p = at(context(s, true, 0.0, false), "temperature", t);
            check(
                card,
                "computed-volumetric-heat-capacity-sample-standard-deviation",
                C,
                &p,
                sd(&v),
            );
            check(
                card,
                "reported-specific-heat-capacity-standard-deviation",
                CP,
                &p,
                s.cp_sd[j],
            );
        }
        for j in 0..2 {
            let t = 288.15 + 30.0 * j as f64;
            close(s.k[j], (mean(&s.conductivity[j]) * 1000.0).round() / 1000.0);
            close(
                s.k_sd[j],
                (sd(&s.conductivity[j]) * 1000.0).round() / 1000.0,
            );
            check(card, "thermal-conductivity", K, &point(s, true, t), s.k[j]);
            check(
                card,
                "apparent-thermal-conductivity",
                K,
                &at(context(s, true, 1.0, false), "test-mean-temperature", t),
                s.k[j],
            );
            check(
                card,
                "reported-apparent-thermal-conductivity-standard-deviation",
                K,
                &at(context(s, true, 0.0, false), "test-mean-temperature", t),
                s.k_sd[j],
            );
        }
        check(
            card,
            "volumetric-heat-capacity",
            C,
            &point(s, true, 298.15),
            mean(&s.capacity.map(|v| 0.5 * (v[1] + v[2]))),
        );
    }
}

#[test]
fn g0_conditioning_density_label_conflicts_and_material_axes_cannot_be_erased() {
    for s in &SOURCES {
        let pack = compiled(s.slug);
        let card = pack.card();
        for (axis, value) in [
            (s.family, 0.0),
            ("source-dried", 0.0),
            ("source-unconditioned", 1.0),
            ("source-grade-known", 1.0),
            ("source-material-axis-known", 1.0),
            ("source-moisture-content-known", 1.0),
            ("source-specimens-matched", 1.0),
            ("source-specimen-count", 9.0),
        ] {
            refuses(
                card,
                "volumetric-heat-capacity",
                C,
                &point(s, true, 298.15).with(axis, value).unwrap(),
            );
        }
        let wrong_header = point(s, true, 298.15)
            .with(
                "source-density-header-consistent",
                1.0 - s.header_consistent,
            )
            .unwrap();
        refuses(card, "specific-heat-capacity", CP, &wrong_header);
        refuses(card, "reference-normalization-density", RHO, &wrong_header);
        // A normalization-label ambiguity is irrelevant to native volumetric calorimetry.
        check(
            card,
            "volumetric-heat-capacity",
            C,
            &wrong_header,
            mean(&s.capacity.map(|v| 0.5 * (v[1] + v[2]))),
        );
        refuses(
            card,
            "thermal-conductivity",
            K,
            &at(context(s, true, 1.0, false), "temperature", 298.15),
        );
        refuses(
            card,
            "apparent-thermal-conductivity",
            K,
            &point(s, true, 298.15),
        );
        refuses(
            card,
            "volumetric-heat-capacity",
            CP,
            &point(s, true, 298.15),
        );
        refuses(
            card,
            "volumetric-heat-capacity",
            C,
            &point(s, true, 298.15)
                .with_quantity(
                    "temperature",
                    semantic(QuantityKind::TemperatureDifference),
                    298.15,
                )
                .unwrap(),
        );
        for t in [283.14, 313.16] {
            refuses(card, "volumetric-heat-capacity", C, &point(s, true, t));
        }
        for t in [288.14, 318.16] {
            refuses(card, "thermal-conductivity", K, &point(s, true, t));
        }
        refuses(
            card,
            "computed-volumetric-heat-capacity-sample-standard-deviation",
            C,
            &at(context(s, true, 0.0, false), "temperature", 298.15),
        );
        refuses(
            card,
            "reference-normalization-density",
            RHO,
            &context(s, true, 1.0, false)
                .with("source-density-temperature-known", 1.0)
                .unwrap(),
        );
        for name in [
            "density",
            "young-modulus",
            "thermal-expansion-coefficient",
            "latent-heat",
        ] {
            assert!(card.claims().claims_for(name).is_empty());
        }
    }
}

#[test]
fn g1_g3_dried_populations_change_heat_storage_and_shared_layer_conduction() {
    let config = LayerConductionConfig {
        max_iterations: 128,
        absolute_heat_tolerance_w: 1e-12,
        relative_heat_tolerance: 1e-12,
    };
    for s in &SOURCES {
        let dry_pack = compiled(s.slug);
        let wet_pack = compiled(s.wet_slug);
        let dry = prepared(dry_pack.card(), s, true);
        let wet = prepared(wet_pack.card(), s, false);
        let original = dry.clone();
        let c20 = mean(&s.capacity.map(|v| v[1]));
        let c30 = mean(&s.capacity.map(|v| v[2]));
        let k25 = s.k[0] + (s.k[1] - s.k[0]) / 3.0;
        let sample = dry.at_temperature(298.15).unwrap();
        close(sample.volumetric_heat_capacity_j_m3_k(), 0.5 * (c20 + c30));
        close(sample.diffusivity_m2_s(), k25 / (0.5 * (c20 + c30)));
        close(
            dry.storage_change_j_m3(293.15, 303.15).unwrap(),
            5.0 * (c20 + c30),
        );
        close(
            dry.storage_change_j_m3(303.15, 293.15).unwrap(),
            -5.0 * (c20 + c30),
        );
        // One authored fixed volume receives 600 J. This compares source populations;
        // it does not transform a wet specimen into a dry one or assume equal masses.
        let delivered = 600.0 / 0.002;
        let slope = (c30 - c20) / 10.0;
        let delta = 2.0 * delivered / (c20 + (c20 * c20 + 2.0 * slope * delivered).sqrt());
        let heated = dry
            .temperature_after_storage_change(293.15, delivered)
            .unwrap();
        close(heated, 293.15 + delta);
        assert!(
            heated
                > wet
                    .temperature_after_storage_change(293.15, delivered)
                    .unwrap()
        );
        assert!(
            dry.storage_change_j_m3(293.15, 303.15).unwrap()
                < wet.storage_change_j_m3(293.15, 303.15).unwrap()
        );
        let entropy = slope * 10.0 + (c20 - slope * 293.15) * (303.15_f64 / 293.15).ln();
        close(
            dry.sensible_entropy_change_j_m3_k(293.15, 303.15).unwrap(),
            entropy,
        );
        let curve = |path: &PreparedVolumetricThermalTransport| {
            ConductivityCurve::new(
                path.source_lines()[1]
                    .knots
                    .iter()
                    .map(|(t, a)| [*t, a.evidence.value.value])
                    .collect(),
            )
            .unwrap()
        };
        let dk = curve(&dry);
        let wk = curve(&wet);
        let dry_layers = TwoLayerConduction::new([&dk; 2], [0.01; 2], 0.01, 0.0).unwrap();
        let wet_layers = TwoLayerConduction::new([&wk; 2], [0.01; 2], 0.01, 0.0).unwrap();
        let q = dry_layers
            .evaluate([293.15, 303.15], config, || Ok(()))
            .unwrap()
            .power_a_to_b_w;
        close(q, -0.01 / 0.02 * 10.0 * k25);
        assert!(
            q.abs()
                < wet_layers
                    .evaluate([293.15, 303.15], config, || Ok(()))
                    .unwrap()
                    .power_a_to_b_w
                    .abs()
        );
        close(
            dry_layers
                .evaluate([303.15, 293.15], config, || Ok(()))
                .unwrap()
                .power_a_to_b_w,
            -q,
        );
        let contact = TwoLayerConduction::new([&dk, &wk], [0.01; 2], 0.01, 0.1)
            .unwrap()
            .evaluate([293.15, 303.15], config, || Ok(()))
            .unwrap();
        let no_contact = TwoLayerConduction::new([&dk, &wk], [0.01; 2], 0.01, 0.0)
            .unwrap()
            .evaluate([293.15, 303.15], config, || Ok(()))
            .unwrap();
        assert!(contact.power_a_to_b_w.abs() < no_contact.power_a_to_b_w.abs());
        assert!(dry.temperature_after_storage_change(293.15, 1e9).is_err());
        assert!(dry.at_temperature(283.15).is_err());
        assert!(dry.at_temperature(318.15).is_err());
        for line in dry.source_lines() {
            for (_, answer) in &line.knots {
                dry_pack
                    .card()
                    .claims()
                    .verify_receipt(&answer.receipt)
                    .unwrap();
            }
        }
        assert_eq!(dry, original);
        assert_ne!(dry.identity(), wet.identity());
    }
}
