//! Primary FSRI panels: native volume capacity, actual specimen counts and source-conditioned heat use.
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

type ConductivityGroup = (f64, &'static [f64], f64, f64);
struct Panel {
    id: &'static str,
    year: u16,
    metadata_thickness: f64,
    density_inputs: [f64; 3],
    capacity_readings: &'static [[f64; 4]],
    cp: [f64; 4],
    cp_sd: Option<[f64; 4]>,
    conductivity: &'static [ConductivityGroup],
}

// Independent transcription of the pinned native HFM result rows. Conductivity
// groups carry mean Celsius, native Average Cond readings, published mean and SD.
// MDF has one capacity specimen and four k groups; OSB's upper k group is 65 C.
const PANELS: [Panel; 4] = [
    Panel {
        id: "mdf",
        year: 2021,
        metadata_thickness: 0.0159,
        density_inputs: [740.3654617, 734.5767314, 735.3357691],
        capacity_readings: &[[842070.0, 886846.0, 948692.0, 945432.0]],
        cp: [1143.0, 1204.0, 1288.0, 1283.0],
        cp_sd: None,
        conductivity: &[
            (30.0, &[0.1197, 0.1171, 0.1178], 0.118, 0.001),
            (35.0, &[0.1216, 0.1193, 0.1196], 0.120, 0.001),
            (40.0, &[0.1233, 0.1208, 0.1211], 0.122, 0.001),
            (45.0, &[0.1247, 0.1219, 0.1225], 0.123, 0.001),
        ],
    },
    Panel {
        id: "osb",
        year: 2022,
        metadata_thickness: 0.0159,
        density_inputs: [578.3864597, 602.4083981, 614.3630568],
        capacity_readings: &[
            [821546.0, 850111.0, 894755.0, 929026.0],
            [857658.0, 878342.0, 926534.0, 965576.0],
            [818517.0, 858945.0, 895473.0, 934432.0],
        ],
        cp: [1391.0, 1441.0, 1513.0, 1576.0],
        cp_sd: Some([36.0, 24.0, 30.0, 33.0]),
        conductivity: &[
            (
                15.0,
                &[
                    0.09831, 0.09831, 0.09811, 0.1114, 0.1123, 0.1116, 0.1101, 0.1102, 0.1099,
                ],
                0.107,
                0.006,
            ),
            (
                65.0,
                &[
                    0.1110, 0.1104, 0.1099, 0.1246, 0.1238, 0.1233, 0.1227, 0.1222, 0.1219,
                ],
                0.119,
                0.006,
            ),
        ],
    },
    Panel {
        id: "plywood",
        year: 2021,
        metadata_thickness: 0.0159,
        density_inputs: [512.7, 485.6, 496.4],
        capacity_readings: &[
            [698231.0, 720096.0, 761710.0, 802075.0],
            [624611.0, 654160.0, 687223.0, 719622.0],
            [674474.0, 693532.0, 723167.0, 756086.0],
        ],
        cp: [1336.0, 1383.0, 1453.0, 1524.0],
        cp_sd: Some([75.0, 67.0, 75.0, 83.0]),
        conductivity: &[
            (
                15.0,
                &[
                    0.1200, 0.1204, 0.1203, 0.1079, 0.1082, 0.1082, 0.1099, 0.1102, 0.1103,
                ],
                0.113,
                0.006,
            ),
            (
                45.0,
                &[
                    0.1283, 0.1282, 0.1280, 0.1153, 0.1153, 0.1153, 0.1175, 0.1177, 0.1177,
                ],
                0.120,
                0.006,
            ),
        ],
    },
    Panel {
        id: "particleboard",
        year: 2021,
        metadata_thickness: 0.0189,
        density_inputs: [684.0, 676.0, 671.0],
        capacity_readings: &[
            [941453.0, 995109.0, 1034608.0, 1086834.0],
            [929128.0, 971335.0, 1025477.0, 1074253.0],
            [908319.0, 952548.0, 1000509.0, 1042877.0],
        ],
        cp: [1368.0, 1437.0, 1507.0, 1578.0],
        cp_sd: Some([25.0, 32.0, 26.0, 33.0]),
        conductivity: &[
            (
                15.0,
                &[
                    0.1266, 0.1280, 0.1273, 0.1240, 0.1242, 0.1236, 0.1212, 0.1216, 0.1208,
                ],
                0.124,
                0.003,
            ),
            (
                45.0,
                &[
                    0.1386, 0.1379, 0.1364, 0.1355, 0.1342, 0.1327, 0.1326, 0.1314, 0.1300,
                ],
                0.134,
                0.003,
            ),
        ],
    },
];

fn compile(panel: &Panel) -> NormalizedMaterialCardPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let slug = format!("fsri-{}-unconditioned-{}", panel.id, panel.year);
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let output = std::env::temp_dir().join(format!(
        "fsri-panel-{}-{nonce}-{}.fsmatpk",
        std::process::id(),
        panel.id
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
        if panel.cp_sd.is_some() { 8 } else { 6 }
    );
    let values: usize = pack
        .claims()
        .claims_ordered()
        .map(|(_, claim)| match &claim.value {
            fs_matdb::PropertyValue::Scalar { .. } => 1,
            fs_matdb::PropertyValue::Curve { knots, .. } => knots.len(),
            _ => panic!("unexpected panel representation"),
        })
        .sum();
    assert_eq!(values, if panel.id == "mdf" { 21 } else { 23 });
    let pack = NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: slug,
            phase: "solid".into(),
            process: "FSRI unconditioned panel; exact grade and moisture content unknown".into(),
            revision: 0,
        },
        pack,
    )
    .unwrap();
    let bytes = pack.to_bytes();
    let restored =
        NormalizedMaterialCardPack::from_bytes_verified(pack.content_hash(), &bytes).unwrap();
    assert_eq!(restored.to_bytes(), bytes);
    restored
}

fn semantic(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
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
fn context(p: &Panel, statistic: f64, local: bool) -> QueryPoint {
    let mut q = QueryPoint::new();
    for (axis, value) in [
        (format!("source-fsri-{}", p.id), 1.0),
        ("source-unconditioned".into(), 1.0),
        ("source-grade-known".into(), 0.0),
        ("source-formulation-known".into(), 0.0),
        ("source-moisture-content-known".into(), 0.0),
        ("source-contact-pressure-known".into(), 0.0),
        ("source-statistic-is-mean".into(), statistic),
        (
            "source-capacity-specimen-count".into(),
            p.capacity_readings.len() as f64,
        ),
        ("source-conductivity-specimen-count".into(), 3.0),
        ("source-density-specimen-count".into(), 3.0),
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
        "source-metadata-stock-thickness",
        QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0])),
        p.metadata_thickness,
    )
    .unwrap()
    .with_quantity(
        "source-hfm-temperature-difference",
        semantic(QuantityKind::TemperatureDifference),
        20.0,
    )
    .unwrap()
}
fn at(q: QueryPoint, axis: &str, temperature: f64) -> QueryPoint {
    q.with_quantity(
        axis,
        semantic(QuantityKind::AbsoluteTemperature),
        temperature,
    )
    .unwrap()
}
fn point(p: &Panel, temperature: f64) -> QueryPoint {
    at(context(p, 1.0, true), "temperature", temperature)
}
fn mean(v: &[f64]) -> f64 {
    v.iter().sum::<f64>() / v.len() as f64
}
fn sd(v: &[f64]) -> f64 {
    let m = mean(v);
    (v.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (v.len() - 1) as f64).sqrt()
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-12 * expected.abs().max(1e-10),
        "{actual} != {expected}"
    );
}
fn sample(card: &MaterialCard, name: &str, dims: Dims, q: &QueryPoint, expected: f64) {
    let answer = card
        .claims()
        .query_typed(&key(name, dims), q, SelectionPolicy::SingleClaimOnly)
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
fn refuses(card: &MaterialCard, name: &str, dims: Dims, q: &QueryPoint) {
    assert!(
        card.claims()
            .query_typed(&key(name, dims), q, SelectionPolicy::SingleClaimOnly)
            .is_err(),
        "must refuse {name} at {q:?}"
    );
}

#[test]
fn g1_g3_native_values_counts_and_published_reductions_survive_compiler_and_card_roundtrip() {
    for p in &PANELS {
        let pack = compile(p);
        let card = pack.card();
        let rho = mean(&p.density_inputs);
        sample(
            card,
            "reference-normalization-density",
            RHO,
            &context(p, 1.0, true),
            rho,
        );
        for i in 0..4 {
            let t = 283.15 + 10.0 * i as f64;
            let values: Vec<_> = p.capacity_readings.iter().map(|row| row[i]).collect();
            let volume = mean(&values);
            sample(card, "volumetric-heat-capacity", C, &point(p, t), volume);
            sample(card, "specific-heat-capacity", CP, &point(p, t), p.cp[i]);
            assert_eq!((volume / rho).round_ties_even(), p.cp[i]);
            if let Some(cp_sd) = p.cp_sd {
                let q = at(context(p, 0.0, true), "temperature", t);
                sample(
                    card,
                    "computed-volumetric-heat-capacity-sample-standard-deviation",
                    C,
                    &q,
                    sd(&values),
                );
                sample(
                    card,
                    "reported-specific-heat-capacity-standard-deviation",
                    CP,
                    &q,
                    cp_sd[i],
                );
                assert_eq!((sd(&values) / rho).round_ties_even(), cp_sd[i]);
            } else {
                assert_eq!(values.len(), 1);
                assert!(
                    card.claims()
                        .claims_for("computed-volumetric-heat-capacity-sample-standard-deviation")
                        .is_empty()
                );
                assert!(
                    card.claims()
                        .claims_for("reported-specific-heat-capacity-standard-deviation")
                        .is_empty()
                );
            }
        }
        for &(t, readings, k, scatter) in p.conductivity {
            sample(card, "thermal-conductivity", K, &point(p, t + 273.15), k);
            sample(
                card,
                "apparent-thermal-conductivity",
                K,
                &at(context(p, 1.0, false), "test-mean-temperature", t + 273.15),
                k,
            );
            sample(
                card,
                "reported-apparent-thermal-conductivity-standard-deviation",
                K,
                &at(context(p, 0.0, false), "test-mean-temperature", t + 273.15),
                scatter,
            );
            assert_eq!((mean(readings) * 1000.0).round_ties_even() / 1000.0, k);
            assert_eq!((sd(readings) * 1000.0).round_ties_even() / 1000.0, scatter);
        }
    }
}

#[test]
fn g0_missing_replication_moisture_grade_direction_and_temperature_do_not_become_defaults() {
    for p in &PANELS {
        let pack = compile(p);
        let card = pack.card();
        for (axis, value) in [
            ("source-unconditioned", 0.0),
            ("source-grade-known", 1.0),
            ("source-formulation-known", 1.0),
            ("source-moisture-content-known", 1.0),
            ("source-contact-pressure-known", 1.0),
            ("source-capacity-specimen-count", 2.0),
            ("source-conductivity-specimen-count", 9.0),
            (
                "source-metadata-stock-thickness",
                p.metadata_thickness * 1.1,
            ),
        ] {
            refuses(
                card,
                "volumetric-heat-capacity",
                C,
                &point(p, 308.15).with(axis, value).unwrap(),
            );
        }
        for (axis, value) in [
            ("source-hfm-through-thickness", 0.0),
            ("source-hfm-upward-heat-flow", 0.0),
            ("source-hfm-temperature-difference", 25.0),
            ("mean-test-temperature-as-local-model", 0.0),
        ] {
            refuses(
                card,
                "thermal-conductivity",
                K,
                &point(p, 308.15).with(axis, value).unwrap(),
            );
        }
        refuses(
            card,
            "thermal-conductivity",
            K,
            &at(context(p, 1.0, false), "temperature", 308.15),
        );
        refuses(card, "apparent-thermal-conductivity", K, &point(p, 308.15));
        refuses(card, "volumetric-heat-capacity", CP, &point(p, 308.15));
        for t in [283.14, 313.16] {
            refuses(card, "volumetric-heat-capacity", C, &point(p, t));
        }
        for t in [
            p.conductivity[0].0 + 273.14,
            p.conductivity.last().unwrap().0 + 273.16,
        ] {
            refuses(card, "thermal-conductivity", K, &point(p, t));
        }
        let wrong_delta = point(p, 308.15)
            .with_quantity(
                "source-hfm-temperature-difference",
                semantic(QuantityKind::AbsoluteTemperature),
                20.0,
            )
            .unwrap();
        refuses(card, "thermal-conductivity", K, &wrong_delta);
        let wrong_temperature = point(p, 308.15)
            .with_quantity(
                "temperature",
                semantic(QuantityKind::TemperatureDifference),
                308.15,
            )
            .unwrap();
        refuses(card, "volumetric-heat-capacity", C, &wrong_temperature);
        for name in ["density", "young-modulus", "thermal-expansion-coefficient"] {
            assert!(card.claims().claims_for(name).is_empty());
        }
        let density = card
            .claims()
            .query_typed(
                &key("reference-normalization-density", RHO),
                &context(p, 1.0, true),
                SelectionPolicy::SingleClaimOnly,
            )
            .unwrap();
        assert!(
            card.claims()
                .claim(density.receipt.selected)
                .unwrap()
                .validity
                .bound("temperature")
                .is_none()
        );
        refuses(
            card,
            "reference-normalization-density",
            RHO,
            &context(p, 1.0, true)
                .with("source-density-temperature-known", 1.0)
                .unwrap(),
        );
        if p.cp_sd.is_some() {
            let sd_key = key(
                "computed-volumetric-heat-capacity-sample-standard-deviation",
                C,
            );
            let lo = at(context(p, 0.0, true), "temperature", 283.15);
            let hi = at(context(p, 0.0, true), "temperature", 313.15);
            assert!(
                card.claims()
                    .query_line_typed(&sd_key, "temperature", &lo, &hi, POLICY)
                    .is_err()
            );
            refuses(
                card,
                sd_key.name(),
                C,
                &at(context(p, 0.0, true), "temperature", 308.15),
            );
            refuses(card, sd_key.name(), C, &point(p, 303.15));
        }
    }
}

#[test]
fn g1_native_panel_curves_drive_heat_storage_transport_and_inverse_on_actual_common_support() {
    for p in &PANELS {
        let pack = compile(p);
        let card = pack.card();
        let capacity = key("volumetric-heat-capacity", C);
        let conductivity = key("thermal-conductivity", K);
        let lower = point(p, 303.15);
        let upper = point(p, 313.15);
        let query = VolumetricThermalPathQuery {
            volumetric_heat_capacity: &capacity,
            conductivity: &conductivity,
            temperature_axis: "temperature",
            lower: &lower,
            upper: &upper,
            selections: [POLICY; 2],
        };
        let path = PreparedVolumetricThermalTransport::try_new(card, query).unwrap();
        assert_eq!(path.material_card_identity(), card.content_hash());
        let c30 = mean(
            &p.capacity_readings
                .iter()
                .map(|row| row[2])
                .collect::<Vec<_>>(),
        );
        let c40 = mean(
            &p.capacity_readings
                .iter()
                .map(|row| row[3])
                .collect::<Vec<_>>(),
        );
        let c35 = (c30 + c40) / 2.0;
        let k35 = if p.id == "mdf" {
            0.120
        } else {
            let (t0, _, k0, _) = p.conductivity[0];
            let (t1, _, k1, _) = p.conductivity[1];
            k0 + (k1 - k0) * (35.0 - t0) / (t1 - t0)
        };
        let sample = path.at_temperature(308.15).unwrap();
        close(sample.volumetric_heat_capacity_j_m3_k(), c35);
        close(sample.conductivity_w_m_k(), k35);
        close(sample.diffusivity_m2_s(), k35 / c35);
        close(sample.effusivity_j_m2_k_sqrt_s(), (k35 * c35).sqrt());
        let expected = 10.0 * (c30 + c40) / 2.0;
        close(path.storage_change_j_m3(303.15, 313.15).unwrap(), expected);
        close(path.storage_change_j_m3(313.15, 303.15).unwrap(), -expected);
        let half = 5.0 * (c30 + c35) / 2.0;
        close(
            path.temperature_after_storage_change(303.15, half).unwrap(),
            308.15,
        );
        close(
            path.temperature_after_storage_change(308.15, -half)
                .unwrap(),
            303.15,
        );
        close(
            path.temperature_after_storage_change(303.15, expected)
                .unwrap(),
            313.15,
        );
        assert!(
            path.temperature_after_storage_change(303.15, expected * 1.001)
                .is_err()
        );
        assert!(
            PreparedVolumetricThermalTransport::try_new(
                card,
                VolumetricThermalPathQuery {
                    upper: &point(p, 318.15),
                    ..query
                }
            )
            .is_err()
        );
        if p.id == "mdf" {
            // Cp/Cvol begins at 10 C, but its conductivity begins at 30 C.
            // No hidden 15 C conductivity default may fill that missing interval.
            assert!(
                PreparedVolumetricThermalTransport::try_new(
                    card,
                    VolumetricThermalPathQuery {
                        lower: &point(p, 288.15),
                        ..query
                    }
                )
                .is_err()
            );
            assert!(
                c40 < c30,
                "retain the measured decreasing final capacity segment"
            );
        }
    }
}
