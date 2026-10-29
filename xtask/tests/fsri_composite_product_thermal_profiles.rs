//! Native thermal measurements for acrylic/mineral, wood/polymer and vinyl products.
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
const LAYER: LayerConductionConfig = LayerConductionConfig {
    max_iterations: 128,
    absolute_heat_tolerance_w: 1e-12,
    relative_heat_tolerance: 1e-12,
};

struct Source {
    slug: &'static str,
    family: &'static str,
    dried: f64,
    // Independent literal Results Table ordinates; rows are physical R1/R2/R3.
    capacity: [[f64; 4]; 3],
    density: [f64; 3],
    // Three cycles per specimen, pooled in the original source reduction.
    conductivity: [[f64; 9]; 2],
    cp: [f64; 4],
    cp_sd: [f64; 4],
    k: [f64; 2],
    k_sd: [f64; 2],
}
const SOURCES: [Source; 4] = [
    Source {
        slug: "fsri-acrylic-countertop-dried-2021",
        family: "source-fsri-acrylic-countertop",
        dried: 1.0,
        capacity: [
            [2051024.0, 2102883.0, 2162612.0, 2257228.0],
            [2049294.0, 2123869.0, 2200405.0, 2261385.0],
            [2053829.0, 2128380.0, 2195664.0, 2252849.0],
        ],
        density: [1671.0, 1688.0, 1672.0],
        conductivity: [
            [
                0.3421, 0.3492, 0.3505, 0.3755, 0.3809, 0.3798, 0.3487, 0.3588, 0.3577,
            ],
            [
                0.3717, 0.3722, 0.3709, 0.4047, 0.4048, 0.4047, 0.3823, 0.3814, 0.3798,
            ],
        ],
        cp: [1223.0, 1263.0, 1304.0, 1346.0],
        cp_sd: [1.0, 8.0, 12.0, 3.0],
        k: [0.360, 0.386],
        k_sd: [0.015, 0.015],
    },
    Source {
        slug: "fsri-composite-deck-board-unconditioned-2021",
        family: "source-fsri-composite-deck-board",
        dried: 0.0,
        capacity: [
            [1527547.0, 1583519.0, 1667851.0, 1736954.0],
            [1498702.0, 1554755.0, 1658106.0, 1707120.0],
            [1496233.0, 1576153.0, 1674223.0, 1733543.0],
        ],
        density: [1076.0, 1121.0, 1112.0],
        conductivity: [
            [
                0.2148, 0.2205, 0.2205, 0.2158, 0.2176, 0.2179, 0.2144, 0.2208, 0.2214,
            ],
            [
                0.2322, 0.2332, 0.2334, 0.2317, 0.2313, 0.2315, 0.2330, 0.2338, 0.2344,
            ],
        ],
        cp: [1367.0, 1425.0, 1511.0, 1565.0],
        cp_sd: [16.0, 14.0, 7.0, 15.0],
        k: [0.218, 0.233],
        k_sd: [0.003, 0.001],
    },
    Source {
        slug: "fsri-vinyl-plank-unconditioned-2022",
        family: "source-fsri-vinyl-plank",
        dried: 0.0,
        capacity: [
            [820321.0, 844983.0, 872592.0, 889622.0],
            [832980.0, 856730.0, 886500.0, 903362.0],
            [827957.0, 846168.0, 880067.0, 896791.0],
        ],
        density: [645.0, 676.0, 670.0],
        conductivity: [
            [
                0.06814, 0.06848, 0.06895, 0.07021, 0.07147, 0.07157, 0.06971, 0.07050, 0.07029,
            ],
            [
                0.07316, 0.07334, 0.07384, 0.07513, 0.07545, 0.07534, 0.07464, 0.07480, 0.07456,
            ],
        ],
        cp: [1246.0, 1280.0, 1326.0, 1351.0],
        cp_sd: [10.0, 10.0, 10.0, 10.0],
        k: [0.070, 0.074],
        k_sd: [0.001, 0.001],
    },
    Source {
        slug: "fsri-vinyl-plank-dried-2022",
        family: "source-fsri-vinyl-plank",
        dried: 1.0,
        capacity: [
            [760586.0, 787450.0, 809263.0, 829403.0],
            [775535.0, 799497.0, 825817.0, 849161.0],
            [758191.0, 788195.0, 821738.0, 834387.0],
        ],
        density: [611.0, 668.0, 651.0],
        conductivity: [
            [
                0.06778, 0.06884, 0.06899, 0.06986, 0.07083, 0.07122, 0.06857, 0.06953, 0.06997,
            ],
            [
                0.07260, 0.07301, 0.07311, 0.07429, 0.07460, 0.07461, 0.07306, 0.07303, 0.07322,
            ],
        ],
        cp: [1189.0, 1231.0, 1273.0, 1302.0],
        cp_sd: [15.0, 10.0, 13.0, 16.0],
        k: [0.070, 0.074],
        k_sd: [0.001, 0.001],
    },
];

fn compiled(source: &Source) -> NormalizedMaterialCardPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let output = std::env::temp_dir().join(format!(
        "fsri-product-{}-{nonce}.fsmatpk",
        std::process::id()
    ));
    let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(root.join(format!("data/matdb/seed-v1/{}/manifest.tsv", source.slug)))
        .arg("--out")
        .arg(&output)
        .env("CARGO_WORKSPACE_DIR", root)
        .current_dir(root)
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
            chemistry: source.slug.into(),
            phase: "solid".into(),
            process: format!(
                "source product population; dried={}; unknown formulation and layer structure",
                source.dried
            ),
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
    let mut point = QueryPoint::new();
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
fn point(source: &Source, kelvin: f64) -> QueryPoint {
    at(context(source, 1.0, true), "temperature", kelvin)
}
fn key(name: &str, dims: Dims) -> PropertyKey {
    if name == "thermal-conductivity" || name == "apparent-thermal-conductivity" {
        PropertyKey::with_quantity(name, semantic(QuantityKind::ThermalConductivity))
    } else {
        PropertyKey::new(name, dims)
    }
}
fn mean(values: &[f64]) -> f64 {
    values.iter().sum::<f64>() / values.len() as f64
}
fn sd(values: &[f64]) -> f64 {
    let m = mean(values);
    (values.iter().map(|v| (v - m) * (v - m)).sum::<f64>() / (values.len() - 1) as f64).sqrt()
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-11 * expected.abs().max(1e-10),
        "{actual} != {expected}"
    );
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
fn g1_g3_all_native_values_reductions_units_and_receipts_survive_compilation() {
    for source in &SOURCES {
        let pack = compiled(source);
        let card = pack.card();
        let density = mean(&source.density);
        check(
            card,
            "reference-normalization-density",
            RHO,
            &context(source, 1.0, false),
            density,
        );
        for j in 0..4 {
            let samples = source.capacity.map(|row| row[j]);
            let (capacity, scatter) = (mean(&samples), sd(&samples));
            let t = 283.15 + 10.0 * j as f64;
            check(
                card,
                "volumetric-heat-capacity",
                C,
                &point(source, t),
                capacity,
            );
            check(
                card,
                "specific-heat-capacity",
                CP,
                &point(source, t),
                source.cp[j],
            );
            check(
                card,
                "computed-volumetric-heat-capacity-sample-standard-deviation",
                C,
                &at(context(source, 0.0, false), "temperature", t),
                scatter,
            );
            check(
                card,
                "reported-specific-heat-capacity-standard-deviation",
                CP,
                &at(context(source, 0.0, false), "temperature", t),
                source.cp_sd[j],
            );
            assert_eq!((capacity / density).round_ties_even(), source.cp[j]);
            assert_eq!((scatter / density).round_ties_even(), source.cp_sd[j]);
        }
        for j in 0..2 {
            let t = [288.15, 318.15][j];
            check(
                card,
                "thermal-conductivity",
                K,
                &point(source, t),
                source.k[j],
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
                "reported-apparent-thermal-conductivity-standard-deviation",
                K,
                &at(context(source, 0.0, false), "test-mean-temperature", t),
                source.k_sd[j],
            );
            assert_eq!(
                (1000.0 * mean(&source.conductivity[j])).round_ties_even() / 1000.0,
                source.k[j]
            );
            assert_eq!(
                (1000.0 * sd(&source.conductivity[j])).round_ties_even() / 1000.0,
                source.k_sd[j]
            );
        }
    }
}

#[test]
fn g0_product_construction_conditioning_and_apparent_temperature_remain_required() {
    for source in &SOURCES {
        let pack = compiled(source);
        let card = pack.card();
        for (axis, value) in [
            (source.family, 0.0),
            ("source-dried", 1.0 - source.dried),
            ("source-unconditioned", source.dried),
            ("source-formulation-known", 1.0),
            ("source-filler-fraction-known", 1.0),
            ("source-layer-construction-known", 1.0),
            ("source-specimens-matched", 1.0),
            ("source-contact-pressure-known", 1.0),
            ("source-moisture-content-known", 1.0),
            ("source-specimen-count", 9.0),
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
        for (axis, value) in [
            ("mean-test-temperature-as-local-model", 0.0),
            ("source-hfm-through-thickness", 0.0),
            ("source-hfm-temperature-difference", 25.0),
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
        for t in [283.14, 313.16] {
            refuses(card, "volumetric-heat-capacity", C, &point(source, t));
        }
        for t in [288.14, 318.16] {
            refuses(card, "thermal-conductivity", K, &point(source, t));
        }
        refuses(
            card,
            "computed-volumetric-heat-capacity-sample-standard-deviation",
            C,
            &at(context(source, 0.0, false), "temperature", 298.15),
        );
        refuses(
            card,
            "reference-normalization-density",
            RHO,
            &context(source, 1.0, false)
                .with("source-density-temperature-known", 1.0)
                .unwrap(),
        );
        for name in ["density", "young-modulus", "thermal-expansion-coefficient"] {
            assert!(card.claims().claims_for(name).is_empty());
        }
    }
}

#[test]
fn g1_g3_conditioned_native_storage_and_integrated_conduction_change_product_response() {
    let mut heated = Vec::new();
    let mut powers = Vec::new();
    for source in &SOURCES {
        let pack = compiled(source);
        let card = pack.card();
        let (ckey, kkey) = (
            key("volumetric-heat-capacity", C),
            key("thermal-conductivity", K),
        );
        let (lower, upper) = (point(source, 288.15), point(source, 313.15));
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
        let c20 = mean(&source.capacity.map(|x| x[1]));
        let c30 = mean(&source.capacity.map(|x| x[2]));
        let slope = (c30 - c20) / 10.0;
        let k25 = source.k[0] + (source.k[1] - source.k[0]) / 3.0;
        let sample = path.at_temperature(298.15).unwrap();
        close(sample.volumetric_heat_capacity_j_m3_k(), 0.5 * (c20 + c30));
        close(sample.conductivity_w_m_k(), k25);
        close(sample.diffusivity_m2_s(), k25 / (0.5 * (c20 + c30)));
        let storage = 5.0 * (c20 + c30);
        close(path.storage_change_j_m3(293.15, 303.15).unwrap(), storage);
        close(path.storage_change_j_m3(303.15, 293.15).unwrap(), -storage);
        // Same authored fixed reference volume and sensible-caloric constraint.
        // The dry/wet comparison does not turn these populations into a drying model.
        let delivered = 600.0 / 0.002;
        let delta = 2.0 * delivered / (c20 + (c20 * c20 + 2.0 * slope * delivered).sqrt());
        let temperature = path
            .temperature_after_storage_change(293.15, delivered)
            .unwrap();
        close(temperature, 293.15 + delta);
        heated.push(temperature);
        let entropy = slope * 10.0 + (c20 - slope * 293.15) * (303.15_f64 / 293.15).ln();
        close(
            path.sensible_entropy_change_j_m3_k(293.15, 303.15).unwrap(),
            entropy,
        );
        let conductivity = ConductivityCurve::new(
            path.source_lines()[1]
                .knots
                .iter()
                .map(|(t, a)| [*t, a.evidence.value.value])
                .collect(),
        )
        .unwrap();
        let layers = TwoLayerConduction::new([&conductivity; 2], [0.01; 2], 0.01, 0.0).unwrap();
        let flux = layers.evaluate([293.15, 303.15], LAYER, || Ok(())).unwrap();
        // Independent integral of affine k across 20..30 C and total length .02 m.
        close(flux.power_a_to_b_w, -0.01 / 0.02 * 10.0 * k25);
        close(
            layers
                .evaluate([303.15, 293.15], LAYER, || Ok(()))
                .unwrap()
                .power_a_to_b_w,
            -flux.power_a_to_b_w,
        );
        let contact = TwoLayerConduction::new([&conductivity; 2], [0.01; 2], 0.01, 0.1)
            .unwrap()
            .evaluate([293.15, 303.15], LAYER, || Ok(()))
            .unwrap();
        assert!(contact.power_a_to_b_w.abs() < flux.power_a_to_b_w.abs());
        powers.push(flux.power_a_to_b_w.abs());
        assert!(
            PreparedVolumetricThermalTransport::try_new(
                card,
                VolumetricThermalPathQuery {
                    lower: &point(source, 283.15),
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
        assert!(path.temperature_after_storage_change(293.15, 1e9).is_err());
        for line in path.source_lines() {
            for (_, answer) in &line.knots {
                card.claims().verify_receipt(&answer.receipt).unwrap();
            }
        }
        assert_eq!(path, original);
    }
    assert!(powers[0] > powers[1] && powers[1] > powers[2]);
    // Published k rounds to the same values; measured native storage differs.
    assert_eq!(powers[2], powers[3]);
    assert!(heated[3] > heated[2]);
    assert!(heated[2] > heated[1] && heated[1] > heated[0]);
}
