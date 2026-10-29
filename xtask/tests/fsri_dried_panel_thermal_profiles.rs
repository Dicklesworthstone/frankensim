//! G0/G1/G3: source conditioning changes the actual supported panel heat response.
#![deny(unsafe_code)]

use fs_matdb::{
    ClaimSelection, MaterialCard, MaterialStateId, NormalizedMaterialCardPack, NormalizedPack,
    PropertyKey, PropertyValue, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::thermal::transport::volumetric::{
    PreparedVolumetricThermalTransport as ThermalPath, VolumetricThermalPathQuery,
    timeline::FixedVolumeThermalLaw,
};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionSimulation, Event,
    MaterialId, StepConfig, ThermalLaw, TimedEvent, series_contact_conductance_w_k,
};
use std::{path::Path, process::Command};

const C: Dims = Dims([-1, 1, -2, -1, 0, 0]);
const CP: Dims = Dims([2, 0, -2, -1, 0, 0]);
const K: Dims = Dims([1, 1, -3, -1, 0, 0]);
const RHO: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const POLICY: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);

struct Panel {
    id: &'static str,
    year: u16,
    stock: f64,
    rho: [f64; 3],
    capacity: &'static [[f64; 4]],
    cp: [f64; 4],
    cp_sd: [f64; 4],
    k_temperature: [f64; 2],
    k_readings: [[f64; 9]; 2],
    k: [f64; 2],
    k_sd: [f64; 2],
    wet_capacity_20_30: [f64; 2],
}

// Independently transcribed native dry HFM rows and reported CSV density columns.
// A conductivity group has nine readings from THREE specimens, not nine coupons.
const PANELS: [Panel; 4] = [
    Panel {
        id: "mdf",
        year: 2021,
        stock: 0.0159,
        rho: [701.6333691, 701.4060172, 702.8080209],
        capacity: &[],
        cp: [0.0; 4],
        cp_sd: [0.0; 4],
        k_temperature: [30.0, 45.0],
        k_readings: [
            [
                0.1078, 0.1079, 0.1078, 0.1069, 0.1073, 0.1073, 0.1071, 0.1074, 0.1074,
            ],
            [
                0.1110, 0.1110, 0.1109, 0.1101, 0.1101, 0.1101, 0.1102, 0.1103, 0.1103,
            ],
        ],
        k: [0.107, 0.110],
        k_sd: [0.0, 0.0],
        wet_capacity_20_30: [0.0; 2],
    },
    Panel {
        id: "osb",
        year: 2022,
        stock: 0.0159,
        rho: [552.4637966, 575.4162888, 587.0479533],
        capacity: &[
            [694699.0, 718164.0, 754813.0, 784823.0],
            [718127.0, 739097.0, 777132.0, 806901.0],
            [694524.0, 721566.0, 751244.0, 780264.0],
        ],
        cp: [1229.0, 1271.0, 1331.0, 1383.0],
        cp_sd: [24.0, 20.0, 25.0, 25.0],
        k_temperature: [15.0, 65.0],
        k_readings: [
            [
                0.09198, 0.09223, 0.09220, 0.09877, 0.09920, 0.09931, 0.09937, 0.1001, 0.1004,
            ],
            [
                0.1018, 0.1018, 0.1017, 0.1090, 0.1089, 0.1088, 0.1097, 0.1097, 0.1095,
            ],
        ],
        k: [0.097, 0.107],
        k_sd: [0.004, 0.004],
        wet_capacity_20_30: [
            (850111.0 + 878342.0 + 858945.0) / 3.0,
            (894755.0 + 926534.0 + 895473.0) / 3.0,
        ],
    },
    Panel {
        id: "plywood",
        year: 2021,
        stock: 0.0159,
        rho: [477.6, 451.2, 460.7],
        capacity: &[
            [561361.0, 587984.0, 605953.0, 626381.0],
            [508948.0, 528002.0, 538640.0, 564851.0],
            [519698.0, 538665.0, 558273.0, 588353.0],
        ],
        cp: [1144.0, 1191.0, 1226.0, 1281.0],
        cp_sd: [60.0, 69.0, 75.0, 67.0],
        k_temperature: [15.0, 45.0],
        k_readings: [
            [
                0.1016, 0.1020, 0.1021, 0.09131, 0.09138, 0.09143, 0.09263, 0.09274, 0.09271,
            ],
            [
                0.1085, 0.1086, 0.1088, 0.09757, 0.09749, 0.09743, 0.09879, 0.09891, 0.09871,
            ],
        ],
        k: [0.095, 0.102],
        k_sd: [0.005, 0.005],
        wet_capacity_20_30: [
            (720096.0 + 654160.0 + 693532.0) / 3.0,
            (761710.0 + 687223.0 + 723167.0) / 3.0,
        ],
    },
    Panel {
        id: "particleboard",
        year: 2021,
        stock: 0.0189,
        rho: [634.0, 626.0, 621.0],
        capacity: &[
            [771864.0, 791145.0, 827181.0, 864787.0],
            [752149.0, 793034.0, 813466.0, 860715.0],
            [742773.0, 773100.0, 802549.0, 822120.0],
        ],
        cp: [1205.0, 1253.0, 1299.0, 1354.0],
        cp_sd: [24.0, 18.0, 20.0, 38.0],
        k_temperature: [15.0, 45.0],
        k_readings: [
            [
                0.1112, 0.1114, 0.1113, 0.1108, 0.1112, 0.1111, 0.1098, 0.1102, 0.1100,
            ],
            [
                0.1182, 0.1181, 0.1180, 0.1177, 0.1177, 0.1176, 0.1163, 0.1163, 0.1162,
            ],
        ],
        k: [0.111, 0.117],
        k_sd: [0.001, 0.001],
        wet_capacity_20_30: [
            (995109.0 + 971335.0 + 952548.0) / 3.0,
            (1034608.0 + 1025477.0 + 1000509.0) / 3.0,
        ],
    },
];

fn compile(p: &Panel, dried: bool) -> NormalizedMaterialCardPack {
    let condition = if dried { "dried" } else { "unconditioned" };
    let slug = format!("fsri-{}-{condition}-{}", p.id, p.year);
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let output =
        std::env::temp_dir().join(format!("{slug}-{}-{nonce}.fsmatpk", std::process::id()));
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
        "{}",
        String::from_utf8_lossy(&run.stderr)
    );
    let bytes = std::fs::read(output).unwrap();
    let raw = NormalizedPack::from_bytes(&bytes).unwrap();
    let raw = NormalizedPack::from_bytes_verified(raw.content_hash(), &bytes).unwrap();
    assert_eq!(raw.pack_id(), slug);
    if dried {
        assert_eq!(
            raw.claims().claims_ordered().count(),
            if p.capacity.is_empty() { 4 } else { 8 }
        );
        let values: usize = raw
            .claims()
            .claims_ordered()
            .map(|(_, c)| match &c.value {
                PropertyValue::Scalar { .. } => 1,
                PropertyValue::Curve { knots, .. } => knots.len(),
                _ => panic!("unexpected source payload"),
            })
            .sum();
        assert_eq!(values, if p.capacity.is_empty() { 7 } else { 23 });
    }
    let pack = NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: slug,
            phase: "solid".into(),
            process: format!("FSRI {condition} panel"),
            revision: 0,
        },
        raw,
    )
    .unwrap();
    NormalizedMaterialCardPack::from_bytes_verified(pack.content_hash(), &pack.to_bytes()).unwrap()
}
fn quantity(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn key(name: &str, dims: Dims) -> PropertyKey {
    if matches!(
        name,
        "thermal-conductivity" | "apparent-thermal-conductivity"
    ) {
        PropertyKey::with_quantity(name, quantity(QuantityKind::ThermalConductivity))
    } else {
        PropertyKey::new(name, dims)
    }
}
fn context(p: &Panel, dried: bool, statistic: f64, local: bool) -> QueryPoint {
    let mut q = QueryPoint::new();
    let count = if dried {
        p.capacity.len() as f64
    } else if p.id == "mdf" {
        1.0
    } else {
        3.0
    };
    for (axis, value) in [
        (format!("source-fsri-{}", p.id), 1.0),
        ("source-unconditioned".into(), if dried { 0.0 } else { 1.0 }),
        ("source-grade-known".into(), 0.0),
        ("source-formulation-known".into(), 0.0),
        ("source-moisture-content-known".into(), 0.0),
        ("source-contact-pressure-known".into(), 0.0),
        ("source-statistic-is-mean".into(), statistic),
        ("source-capacity-specimen-count".into(), count),
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
    if dried {
        q = q
            .with_quantity("source-dried", QuantitySpec::dimensional(Dims::NONE), 1.0)
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
        p.stock,
    )
    .unwrap()
    .with_quantity(
        "source-hfm-temperature-difference",
        quantity(QuantityKind::TemperatureDifference),
        20.0,
    )
    .unwrap()
}
fn at(q: QueryPoint, axis: &str, t: f64) -> QueryPoint {
    q.with_quantity(axis, quantity(QuantityKind::AbsoluteTemperature), t)
        .unwrap()
}
fn point(p: &Panel, dried: bool, t: f64) -> QueryPoint {
    at(context(p, dried, 1.0, true), "temperature", t)
}
fn mean(v: &[f64]) -> f64 {
    v.iter().sum::<f64>() / v.len() as f64
}
fn sd(v: &[f64]) -> f64 {
    let average = mean(v);
    (v.iter().map(|x| (x - average) * (x - average)).sum::<f64>() / (v.len() - 1) as f64).sqrt()
}
fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 3e-12 * expected.abs().max(1e-10),
        "{actual} != {expected}"
    );
}
fn sample(card: &MaterialCard, name: &str, dims: Dims, q: &QueryPoint) -> f64 {
    let answer = card
        .claims()
        .query_typed(&key(name, dims), q, SelectionPolicy::SingleClaimOnly)
        .unwrap();
    assert_eq!(answer.evidence.value.dims, dims);
    assert_eq!(
        answer.evidence.value.uncertainty,
        UncertaintyModel::Unstated
    );
    assert!(answer.receipt.observation_backed && !answer.receipt.source_hashes.is_empty());
    card.claims().verify_receipt(&answer.receipt).unwrap();
    answer.evidence.value.value
}
fn refuses(card: &MaterialCard, name: &str, dims: Dims, q: &QueryPoint) {
    assert!(
        card.claims()
            .query_typed(&key(name, dims), q, SelectionPolicy::SingleClaimOnly)
            .is_err()
    );
}
fn path(
    card: &MaterialCard,
    p: &Panel,
    dried: bool,
) -> Result<ThermalPath, fs_material::thermal::transport::ThermalTransportError> {
    ThermalPath::try_new(
        card,
        VolumetricThermalPathQuery {
            volumetric_heat_capacity: &key("volumetric-heat-capacity", C),
            conductivity: &key("thermal-conductivity", K),
            temperature_axis: "temperature",
            lower: &point(p, dried, if p.id == "mdf" { 303.15 } else { 288.15 }),
            upper: &point(p, dried, 313.15),
            selections: [POLICY; 2],
        },
    )
}

#[test]
fn g1_dried_observations_and_reductions_survive_compiler_and_card_roundtrip() {
    for p in &PANELS {
        let pack = compile(p, true);
        let card = pack.card();
        let rho = mean(&p.rho);
        close(
            sample(
                card,
                "reference-normalization-density",
                RHO,
                &context(p, true, 1.0, true),
            ),
            rho,
        );
        for i in 0..4 {
            if p.capacity.is_empty() {
                break;
            }
            let values: Vec<_> = p.capacity.iter().map(|v| v[i]).collect();
            let q = point(p, true, 283.15 + 10.0 * i as f64);
            close(
                sample(card, "volumetric-heat-capacity", C, &q),
                mean(&values),
            );
            close(sample(card, "specific-heat-capacity", CP, &q), p.cp[i]);
            assert_eq!((mean(&values) / rho).round_ties_even(), p.cp[i]);
            let q = q.with("source-statistic-is-mean", 0.0).unwrap();
            close(
                sample(
                    card,
                    "computed-volumetric-heat-capacity-sample-standard-deviation",
                    C,
                    &q,
                ),
                sd(&values),
            );
            close(
                sample(
                    card,
                    "reported-specific-heat-capacity-standard-deviation",
                    CP,
                    &q,
                ),
                p.cp_sd[i],
            );
            assert_eq!((sd(&values) / rho).round_ties_even(), p.cp_sd[i]);
        }
        for i in 0..2 {
            let t = p.k_temperature[i] + 273.15;
            close(
                sample(card, "thermal-conductivity", K, &point(p, true, t)),
                p.k[i],
            );
            close(
                sample(
                    card,
                    "apparent-thermal-conductivity",
                    K,
                    &at(context(p, true, 1.0, false), "test-mean-temperature", t),
                ),
                p.k[i],
            );
            close(
                sample(
                    card,
                    "reported-apparent-thermal-conductivity-standard-deviation",
                    K,
                    &at(context(p, true, 0.0, false), "test-mean-temperature", t),
                ),
                p.k_sd[i],
            );
            assert_eq!(
                (mean(&p.k_readings[i]) * 1000.0).round_ties_even() / 1000.0,
                p.k[i]
            );
            assert_eq!(
                (sd(&p.k_readings[i]) * 1000.0).round_ties_even() / 1000.0,
                p.k_sd[i]
            );
            assert!(sd(&p.k_readings[i]) > 0.0); // A rounded zero is not zero raw scatter.
        }
    }
}

#[test]
fn g0_conditioning_and_absent_dried_capacity_cannot_be_filled_by_defaults() {
    for p in &PANELS {
        let pack = compile(p, true);
        let card = pack.card();
        let q = point(p, true, 308.15);
        for (axis, wrong) in [
            ("source-dried", 0.0),
            ("source-unconditioned", 1.0),
            ("source-grade-known", 1.0),
            ("source-moisture-content-known", 1.0),
            ("source-contact-pressure-known", 1.0),
            ("source-conductivity-specimen-count", 9.0),
            ("source-hfm-through-thickness", 0.0),
            ("source-hfm-temperature-difference", 25.0),
            ("mean-test-temperature-as-local-model", 0.0),
        ] {
            refuses(
                card,
                "thermal-conductivity",
                K,
                &q.clone().with(axis, wrong).unwrap(),
            );
        }
        refuses(card, "thermal-conductivity", K, &point(p, false, 308.15));
        refuses(
            card,
            "thermal-conductivity",
            K,
            &at(context(p, true, 1.0, false), "temperature", 308.15),
        );
        refuses(card, "apparent-thermal-conductivity", K, &q);
        for t in [p.k_temperature[0] + 273.14, p.k_temperature[1] + 273.16] {
            refuses(card, "thermal-conductivity", K, &point(p, true, t));
        }
        let wrong_kind = q
            .clone()
            .with_quantity(
                "source-hfm-temperature-difference",
                quantity(QuantityKind::AbsoluteTemperature),
                20.0,
            )
            .unwrap();
        refuses(card, "thermal-conductivity", K, &wrong_kind);
        for name in ["density", "thermal-expansion-coefficient", "young-modulus"] {
            assert!(card.claims().claims_for(name).is_empty());
        }
        let density = card.claims().claims_for("reference-normalization-density")[0].1;
        assert!(density.validity.bound("temperature").is_none());
        refuses(
            card,
            "reference-normalization-density",
            RHO,
            &context(p, true, 1.0, true)
                .with("source-density-temperature-known", 1.0)
                .unwrap(),
        );
        if p.capacity.is_empty() {
            for name in [
                "volumetric-heat-capacity",
                "specific-heat-capacity",
                "computed-volumetric-heat-capacity-sample-standard-deviation",
            ] {
                assert!(card.claims().claims_for(name).is_empty());
            }
            assert!(path(card, p, true).is_err());
        } else {
            for t in [283.14, 313.16] {
                refuses(card, "volumetric-heat-capacity", C, &point(p, true, t));
            }
            refuses(card, "volumetric-heat-capacity", CP, &q);
            refuses(
                card,
                "computed-volumetric-heat-capacity-sample-standard-deviation",
                C,
                &q.with("source-statistic-is-mean", 0.0).unwrap(),
            );
        }
    }
}

fn heat(path: &ThermalPath, rho: f64) -> (f64, f64) {
    // Each condition is a separate, explicitly authored fixed-volume sensible
    // scenario. This comparison does not simulate drying, mass loss or damage.
    let law = FixedVolumeThermalLaw::try_new(path, rho, 293.15)
        .unwrap()
        .with_sensible_entropy();
    let cell = DepositedCell {
        id: CellId(1),
        segment_id: 1,
        body_id: 1,
        material_id: MaterialId(1),
        mass_kg: law.reference_mass_kg(0.002).unwrap(),
        specific_h_j_kg: 0.0,
        geometry: CellGeometry {
            start_m: [0.0; 3],
            end_m: [0.1, 0.0, 0.0],
            volume_m3: 0.002,
            surface_area_m2: 0.1,
        },
        boundary: Boundary {
            ambient_k: 293.15,
            htc_w_m2_k: 0.0,
            reservoir: None,
            source_w: 30.0,
        },
    };
    law.admit_cell(&cell).unwrap();
    let initial = DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(cell),
        }],
        &[(MaterialId(1), &law as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 1,
            max_substeps: 128,
        },
    )
    .unwrap();
    let (heated, energy, entropy) = initial
        .advance_to_with_entropy(
            20.0,
            StepConfig {
                max_dt_s: 1.0,
                energy_tolerance_j: 1e-8,
            },
            || Ok(()),
        )
        .unwrap();
    close(energy.external_heat_j, 600.0);
    close(energy.stored_energy_change_j, 600.0);
    assert!(energy.residual_j.abs() < 1e-9);
    assert!(entropy.has_unclassified_sources);
    assert_eq!(entropy.conduction_generation_j_k, 0.0);
    assert_eq!(entropy.reservoir_entropy_input_j_k, 0.0);
    assert!(entropy.temporal_defect_j_k < 0.0);
    (
        heated.temperature_k(CellId(1)).unwrap(),
        entropy.thermal_storage_change_j_k,
    )
}

#[test]
fn g1_g3_conditioned_sources_change_the_existing_conduction_and_heating_owners() {
    for p in &PANELS {
        let dry = compile(p, true);
        let wet = compile(p, false);
        // Declared equal geometry and original 20 K plate difference. This uses
        // the measured apparent coefficient, without inventing a storage law.
        let effective = |card: &MaterialCard, dried| {
            sample(
                card,
                "apparent-thermal-conductivity",
                K,
                &at(
                    context(p, dried, 1.0, false),
                    "test-mean-temperature",
                    303.15,
                ),
            )
        };
        let dry_k = effective(dry.card(), true);
        let wet_k = effective(wet.card(), false);
        assert!(dry_k < wet_k);
        let dry_g = series_contact_conductance_w_k(0.04, [0.01; 2], [dry_k; 2], 0.0).unwrap();
        let wet_g = series_contact_conductance_w_k(0.04, [0.01; 2], [wet_k; 2], 0.0).unwrap();
        close(dry_g * 20.0, dry_k * 0.04 * 20.0 / 0.02);
        assert!(dry_g < wet_g);
        if p.capacity.is_empty() {
            assert!(path(dry.card(), p, true).is_err());
            continue;
        }
        let dry_path = path(dry.card(), p, true).unwrap();
        let wet_path = path(wet.card(), p, false).unwrap();
        let dry_capacity =
            [1, 2].map(|i| mean(&p.capacity.iter().map(|r| r[i]).collect::<Vec<_>>()));
        let mut final_temperature = Vec::new();
        for (path, [c0, c1]) in [(&dry_path, dry_capacity), (&wet_path, p.wet_capacity_20_30)] {
            let slope = (c1 - c0) / 10.0;
            // Independent inverse of the affine capacity over 20-30 C.
            let delta = 600000.0 / (c0 + (c0 * c0 + 600000.0 * slope).sqrt());
            assert!((0.0..10.0).contains(&delta));
            let expected_entropy =
                0.002 * (slope * delta + (c0 - slope * 293.15) * (delta / 293.15).ln_1p());
            close(
                path.storage_change_j_m3(293.15, 303.15).unwrap(),
                5.0 * (c0 + c1),
            );
            for rho in [500.0, 1000.0] {
                // Authored reference basis cancels; these are not measured densities.
                let (t, s) = heat(path, rho);
                close(t, 293.15 + delta);
                close(s, expected_entropy);
                close(
                    path.temperature_after_storage_change(t, -300000.0).unwrap(),
                    293.15,
                );
                final_temperature.push(t);
            }
        }
        assert_eq!(final_temperature[0], final_temperature[1]);
        assert_eq!(final_temperature[2], final_temperature[3]);
        assert!(final_temperature[0] > final_temperature[2]);
    }
}
