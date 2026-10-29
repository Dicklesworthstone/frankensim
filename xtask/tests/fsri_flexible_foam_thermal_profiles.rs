//! Native flexible-foam heat storage and transport, with conditioning and actual source populations.
#![deny(unsafe_code)]

use fs_matdb::{
    ClaimSelection, MaterialCard, MaterialStateId, NormalizedMaterialCardPack, NormalizedPack,
    PropertyKey, PropertyValue, QueryPoint, SelectionPolicy, UncertaintyModel,
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
type KGroup = (f64, &'static [f64]);

struct Foam {
    id: &'static str,
    condition: &'static str,
    year: u16,
    stock: Option<f64>,
    capacity: &'static [[f64; 4]],
    density: &'static [f64],
    cp: Option<([f64; 4], [f64; 4])>,
    conductivity: &'static [KGroup],
    k_specimens: usize,
    selected_r1: bool,
    claim_count: usize,
    value_count: usize,
}

// Independent selected numerical facts from pinned primary SI result tables.
// The mattress June 1 R1 run is aborted: only June 2 R1/R2/R3 are included.
// Chair R2 capacity is a different 33.11 mm sample; this ~51 mm selection uses R1 only.
const FOAMS: [Foam; 4] = [
    Foam {
        id: "rebond-foam-pad",
        condition: "unconditioned",
        year: 2022,
        stock: Some(0.009),
        capacity: &[
            [141793.0, 142826.0, 144715.0, 147373.0],
            [142485.0, 144755.0, 146454.0, 149966.0],
            [158383.0, 164573.0, 160724.0, 170319.0],
        ],
        density: &[84.0, 84.0, 89.0],
        cp: Some((
            [1722.0, 1759.0, 1758.0, 1820.0],
            [110.0, 141.0, 103.0, 147.0],
        )),
        conductivity: &[
            (
                15.0,
                &[
                    0.03646, 0.03648, 0.03646, 0.03605, 0.03610, 0.03606, 0.03723, 0.03726, 0.03721,
                ],
            ),
            (
                65.0,
                &[
                    0.04360, 0.04356, 0.04352, 0.04298, 0.04305, 0.04300, 0.04415, 0.04420, 0.04415,
                ],
            ),
        ],
        k_specimens: 3,
        selected_r1: false,
        claim_count: 8,
        value_count: 23,
    },
    Foam {
        id: "rebond-foam-pad",
        condition: "dried",
        year: 2022,
        stock: Some(0.009),
        capacity: &[
            [137759.0, 138932.0, 140449.0, 144051.0],
            [137905.0, 138963.0, 140808.0, 143746.0],
            [144136.0, 145083.0, 146560.0, 150139.0],
        ],
        density: &[83.0, 83.0, 88.0],
        cp: Some(([1653.0, 1665.0, 1684.0, 1724.0], [43.0, 42.0, 41.0, 43.0])),
        conductivity: &[
            (
                15.0,
                &[0.03671, 0.03675, 0.03675, 0.03638, 0.03645, 0.03646],
            ),
            (
                45.0,
                &[0.04078, 0.04081, 0.04081, 0.04037, 0.04039, 0.04042],
            ),
        ],
        k_specimens: 2,
        selected_r1: false,
        claim_count: 8,
        value_count: 23,
    },
    Foam {
        id: "puf-mattress-pad",
        condition: "unconditioned",
        year: 2022,
        stock: None,
        capacity: &[
            [41549.0, 42108.0, 41406.0, 38709.0],
            [39886.0, 39495.0, 31250.0, 53684.0],
            [39005.0, 38924.0, 31565.0, 49401.0],
        ],
        density: &[],
        cp: None,
        conductivity: &[
            (
                15.0,
                &[
                    0.03577, 0.03589, 0.03594, 0.03591, 0.03608, 0.03607, 0.03589, 0.03612, 0.03612,
                ],
            ),
            (
                65.0,
                &[
                    0.04561, 0.04561, 0.04562, 0.04586, 0.04586, 0.04584, 0.04597, 0.04598, 0.04596,
                ],
            ),
        ],
        k_specimens: 3,
        selected_r1: false,
        claim_count: 5,
        value_count: 14,
    },
    Foam {
        id: "chair-pu-foam-r1",
        condition: "unconditioned",
        year: 2021,
        stock: Some(0.2),
        capacity: &[[54854.0, 53861.0, 58335.0, 61219.0]],
        density: &[],
        cp: None,
        conductivity: &[
            (25.0, &[0.03762]),
            (30.0, &[0.03874]),
            (35.0, &[0.03969]),
            (40.0, &[0.04055]),
            (45.0, &[0.04155]),
        ],
        k_specimens: 1,
        selected_r1: true,
        claim_count: 3,
        value_count: 14,
    },
];

fn compile(p: &Foam) -> NormalizedMaterialCardPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let slug = format!("fsri-{}-{}-{}", p.id, p.condition, p.year);
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let output = std::env::temp_dir().join(format!(
        "foam-{slug}-{nonce}-{}.fsmatpk",
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
    let pack = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(pack.content_hash(), &bytes).unwrap();
    assert_eq!(pack.to_bytes(), bytes);
    assert_eq!(pack.pack_id(), slug);
    assert_eq!(pack.claims().claims_ordered().count(), p.claim_count);
    assert_eq!(
        pack.claims()
            .claims_ordered()
            .map(|(_, c)| match &c.value {
                PropertyValue::Scalar { .. } => 1,
                PropertyValue::Curve { knots, .. } => knots.len(),
                _ => panic!("unexpected foam representation"),
            })
            .sum::<usize>(),
        p.value_count
    );
    let pack = NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: slug,
            phase: "solid".into(),
            process: format!(
                "FSRI {} foam; source-selected population and unknown formulation",
                p.condition
            ),
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
    ) && dims == K
    {
        PropertyKey::with_quantity(name, semantic(QuantityKind::ThermalConductivity))
    } else {
        PropertyKey::new(name, dims)
    }
}
fn statistic(p: &Foam) -> f64 {
    if p.selected_r1 { 0.0 } else { 1.0 }
}
fn context(p: &Foam, stat: f64, local: bool) -> QueryPoint {
    let mut q = QueryPoint::new();
    for (axis, value) in [
        (format!("source-fsri-{}", p.id), 1.0),
        (
            "source-unconditioned".into(),
            if p.condition == "unconditioned" {
                1.0
            } else {
                0.0
            },
        ),
        (
            "source-oven-dried".into(),
            if p.condition == "dried" { 1.0 } else { 0.0 },
        ),
        ("source-grade-known".into(), 0.0),
        ("source-formulation-known".into(), 0.0),
        ("source-moisture-content-known".into(), 0.0),
        ("source-contact-pressure-known".into(), 0.0),
        ("source-statistic-is-mean".into(), stat),
        (
            "source-capacity-specimen-count".into(),
            p.capacity.len() as f64,
        ),
        (
            "source-conductivity-specimen-count".into(),
            p.k_specimens as f64,
        ),
        (
            "source-density-specimen-count".into(),
            p.density.len() as f64,
        ),
        (
            "source-selected-r1-only".into(),
            if p.selected_r1 { 1.0 } else { 0.0 },
        ),
        (
            "source-metadata-stock-thickness-known".into(),
            if p.stock.is_some() { 1.0 } else { 0.0 },
        ),
        (
            "source-conductivity-reading-count".into(),
            p.conductivity[0].1.len() as f64,
        ),
        ("source-hfm-through-thickness".into(), 1.0),
        ("source-hfm-upward-heat-flow".into(), 1.0),
        ("source-density-temperature-known".into(), 0.0),
    ] {
        q = q
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    if let Some(stock) = p.stock {
        q = q
            .with_quantity(
                "source-metadata-stock-thickness",
                QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0])),
                stock,
            )
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
fn point(p: &Foam, t: f64) -> QueryPoint {
    at(context(p, statistic(p), true), "temperature", t)
}
fn mean(v: &[f64]) -> f64 {
    v.iter().sum::<f64>() / v.len() as f64
}
fn sd(v: &[f64]) -> f64 {
    let m = mean(v);
    (v.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (v.len() - 1) as f64).sqrt()
}
fn reduction(p: &Foam, value: f64) -> f64 {
    if p.cp.is_some() {
        (value * 1000.0).round_ties_even() / 1000.0
    } else {
        value
    }
}
fn capacity(p: &Foam, i: usize) -> f64 {
    mean(&p.capacity.iter().map(|r| r[i]).collect::<Vec<_>>())
}
fn k_sd_name(p: &Foam) -> &'static str {
    if p.cp.is_some() {
        "reported-apparent-thermal-conductivity-standard-deviation"
    } else {
        "computed-apparent-thermal-conductivity-sample-standard-deviation"
    }
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-12 * expected.abs().max(1e-10),
        "{actual} != {expected}"
    );
}
fn sample(card: &MaterialCard, name: &str, dims: Dims, q: &QueryPoint, expected: f64) {
    let r = card
        .claims()
        .query_typed(&key(name, dims), q, SelectionPolicy::SingleClaimOnly)
        .unwrap();
    close(r.evidence.value.value, expected);
    assert_eq!(r.evidence.value.dims, dims);
    assert_eq!(r.evidence.value.uncertainty, UncertaintyModel::Unstated);
    assert!(r.receipt.observation_backed && !r.receipt.source_hashes.is_empty());
    card.claims().verify_receipt(&r.receipt).unwrap();
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
fn native_foam_values_and_distinct_density_reductions_survive_real_compiler() {
    for p in &FOAMS {
        let pack = compile(p);
        let card = pack.card();
        for i in 0..4 {
            let t = 283.15 + 10.0 * i as f64;
            let raw: Vec<_> = p.capacity.iter().map(|r| r[i]).collect();
            sample(
                card,
                "volumetric-heat-capacity",
                C,
                &point(p, t),
                mean(&raw),
            );
            if p.capacity.len() > 1 {
                sample(
                    card,
                    "computed-volumetric-heat-capacity-sample-standard-deviation",
                    C,
                    &at(context(p, 0.0, true), "temperature", t),
                    sd(&raw),
                );
            }
            if let Some((cp, scatter)) = p.cp {
                let rho = mean(p.density);
                sample(card, "specific-heat-capacity", CP, &point(p, t), cp[i]);
                sample(
                    card,
                    "reported-specific-heat-capacity-standard-deviation",
                    CP,
                    &at(context(p, 0.0, true), "temperature", t),
                    scatter[i],
                );
                assert_eq!((mean(&raw) / rho).round_ties_even(), cp[i]);
                assert_eq!((sd(&raw) / rho).round_ties_even(), scatter[i]);
            }
        }
        if !p.density.is_empty() {
            sample(
                card,
                "reference-normalization-density",
                RHO,
                &context(p, statistic(p), true),
                mean(p.density),
            );
        } else {
            for name in [
                "reference-normalization-density",
                "specific-heat-capacity",
                "reported-specific-heat-capacity-standard-deviation",
            ] {
                assert!(
                    card.claims().claims_for(name).is_empty(),
                    "{name} cannot be invented"
                );
            }
        }
        for &(t, raw) in p.conductivity {
            let value = reduction(p, mean(raw));
            sample(
                card,
                "thermal-conductivity",
                K,
                &point(p, t + 273.15),
                value,
            );
            sample(
                card,
                "apparent-thermal-conductivity",
                K,
                &at(
                    context(p, statistic(p), false),
                    "test-mean-temperature",
                    t + 273.15,
                ),
                value,
            );
            if raw.len() > 1 {
                sample(
                    card,
                    k_sd_name(p),
                    K,
                    &at(context(p, 0.0, false), "test-mean-temperature", t + 273.15),
                    reduction(p, sd(raw)),
                );
            }
        }
        if p.selected_r1 {
            assert!(
                card.claims()
                    .claims_for("computed-volumetric-heat-capacity-sample-standard-deviation")
                    .is_empty()
            );
            assert!(card.claims().claims_for(k_sd_name(p)).is_empty());
        }
    }
}

#[test]
fn condition_specimen_unknowns_native_axes_and_finite_support_remain_enforced() {
    for p in &FOAMS {
        let pack = compile(p);
        let card = pack.card();
        for (axis, value) in [
            (
                "source-unconditioned",
                if p.condition == "unconditioned" {
                    0.0
                } else {
                    1.0
                },
            ),
            (
                "source-oven-dried",
                if p.condition == "dried" { 0.0 } else { 1.0 },
            ),
            ("source-grade-known", 1.0),
            ("source-formulation-known", 1.0),
            ("source-moisture-content-known", 1.0),
            ("source-contact-pressure-known", 1.0),
            (
                "source-capacity-specimen-count",
                p.capacity.len() as f64 + 1.0,
            ),
            (
                "source-conductivity-specimen-count",
                p.k_specimens as f64 + 1.0,
            ),
            (
                "source-density-specimen-count",
                p.density.len() as f64 + 1.0,
            ),
            (
                "source-selected-r1-only",
                if p.selected_r1 { 0.0 } else { 1.0 },
            ),
            (
                "source-metadata-stock-thickness-known",
                if p.stock.is_some() { 0.0 } else { 1.0 },
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
            (
                "source-conductivity-reading-count",
                p.conductivity[0].1.len() as f64 + 1.0,
            ),
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
            &at(context(p, statistic(p), false), "temperature", 308.15),
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
        refuses(
            card,
            "thermal-conductivity",
            K,
            &point(p, 308.15)
                .with_quantity(
                    "source-hfm-temperature-difference",
                    semantic(QuantityKind::AbsoluteTemperature),
                    20.0,
                )
                .unwrap(),
        );
        refuses(
            card,
            "volumetric-heat-capacity",
            C,
            &point(p, 308.15)
                .with_quantity(
                    "temperature",
                    semantic(QuantityKind::TemperatureDifference),
                    308.15,
                )
                .unwrap(),
        );
        for name in ["density", "young-modulus", "thermal-expansion-coefficient"] {
            assert!(card.claims().claims_for(name).is_empty());
        }
        if !p.density.is_empty() {
            let answer = card
                .claims()
                .query_typed(
                    &key("reference-normalization-density", RHO),
                    &context(p, statistic(p), true),
                    SelectionPolicy::SingleClaimOnly,
                )
                .unwrap();
            assert!(
                card.claims()
                    .claim(answer.receipt.selected)
                    .unwrap()
                    .validity
                    .bound("temperature")
                    .is_none()
            );
            refuses(
                card,
                "reference-normalization-density",
                RHO,
                &context(p, statistic(p), true)
                    .with("source-density-temperature-known", 1.0)
                    .unwrap(),
            );
        }
        if p.capacity.len() > 1 {
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
fn prepared_foam_curves_drive_volume_heating_inverse_and_transport_without_invented_density() {
    for p in &FOAMS {
        let pack = compile(p);
        let card = pack.card();
        let c_key = key("volumetric-heat-capacity", C);
        let k_key = key("thermal-conductivity", K);
        let lower = point(p, 298.15);
        let upper = point(p, 313.15);
        let query = VolumetricThermalPathQuery {
            volumetric_heat_capacity: &c_key,
            conductivity: &k_key,
            temperature_axis: "temperature",
            lower: &lower,
            upper: &upper,
            selections: [POLICY; 2],
        };
        let path = PreparedVolumetricThermalTransport::try_new(card, query).unwrap();
        assert_eq!(path.material_card_identity(), card.content_hash());
        let c25 = (capacity(p, 1) + capacity(p, 2)) / 2.0;
        let c30 = capacity(p, 2);
        let c40 = capacity(p, 3);
        let c35 = (c30 + c40) / 2.0;
        let bounds = p
            .conductivity
            .windows(2)
            .find(|pair| pair[0].0 <= 35.0 && pair[1].0 >= 35.0)
            .unwrap();
        let k0 = reduction(p, mean(bounds[0].1));
        let k1 = reduction(p, mean(bounds[1].1));
        let k35 = k0 + (k1 - k0) * (35.0 - bounds[0].0) / (bounds[1].0 - bounds[0].0);
        let sample = path.at_temperature(308.15).unwrap();
        close(sample.volumetric_heat_capacity_j_m3_k(), c35);
        close(sample.conductivity_w_m_k(), k35);
        close(sample.diffusivity_m2_s(), k35 / c35);
        close(sample.effusivity_j_m2_k_sqrt_s(), (k35 * c35).sqrt());
        let full = 5.0 * (c25 + c30) / 2.0 + 10.0 * (c30 + c40) / 2.0;
        close(path.storage_change_j_m3(298.15, 313.15).unwrap(), full);
        close(path.storage_change_j_m3(313.15, 298.15).unwrap(), -full);
        close(
            path.temperature_after_storage_change(298.15, full).unwrap(),
            313.15,
        );
        let partial = 5.0 * (c30 + c35) / 2.0;
        close(
            path.temperature_after_storage_change(303.15, partial)
                .unwrap(),
            308.15,
        );
        close(
            path.temperature_after_storage_change(308.15, -partial)
                .unwrap(),
            303.15,
        );
        assert!(
            path.temperature_after_storage_change(298.15, full * 1.001)
                .is_err()
        );
        assert!(
            PreparedVolumetricThermalTransport::try_new(
                card,
                VolumetricThermalPathQuery {
                    upper: &point(p, 313.16),
                    ..query
                }
            )
            .is_err()
        );
        let no_opt_in = at(context(p, statistic(p), false), "temperature", 298.15);
        assert!(
            PreparedVolumetricThermalTransport::try_new(
                card,
                VolumetricThermalPathQuery {
                    lower: &no_opt_in,
                    ..query
                }
            )
            .is_err()
        );
        if p.selected_r1 {
            assert!(
                PreparedVolumetricThermalTransport::try_new(
                    card,
                    VolumetricThermalPathQuery {
                        lower: &point(p, 293.15),
                        ..query
                    }
                )
                .is_err()
            );
        }
        if p.id == "puf-mattress-pad" {
            assert!(
                capacity(p, 2) < capacity(p, 1),
                "retain the measured nonmonotonic capacity"
            );
        }
    }
    // Separate conditions are separate thermal populations, not a drying kinetics model.
    assert!(capacity(&FOAMS[0], 2) > capacity(&FOAMS[1], 2));
}
