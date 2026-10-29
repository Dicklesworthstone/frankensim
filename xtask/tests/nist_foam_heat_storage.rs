//! Primary foam heat capacities drive mass-based heating without k or density.
#![deny(unsafe_code)]

use fs_matdb::{
    ClaimSelection, MaterialCard, MaterialStateId, NormalizedPack, PropertyKey, QueryPoint,
    SelectionPolicy, UncertaintyModel,
};
use fs_material::thermal::transport::storage::{
    HeatStorageAmount, HeatStorageComponentQuery, PreparedHeatStorage,
};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use std::{path::Path, process::Command};

const SINGLE: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);

struct Foam {
    slug: &'static str,
    family: &'static str,
    gas_axis: &'static str,
    gas_value: f64,
    density_label: f64,
    labels_consistent: f64,
    knots: usize,
    coefficients: &'static [f64],
}

// Independent transcription of the original three NIST Cp columns. The
// numerical ordinates were evaluated separately with 70-digit decimal math.
const FOAMS: &[Foam] = &[
    Foam {
        slug: "polyurethane-air-389p25-nist-thermal-reference",
        family: "polyurethane-air",
        gas_axis: "source-cell-gas-air",
        gas_value: 1.0,
        density_label: 389.25,
        labels_consistent: 1.0,
        knots: 28,
        coefficients: &[
            4894.36, -11608.63, 10463.31, -3895.8, -40.0053, 497.4517, -147.7555, 14.19365, 0.0,
        ],
    },
    Foam {
        slug: "polystyrene-9p93-nist-thermal-reference",
        family: "polystyrene",
        gas_axis: "source-cell-gas-known",
        gas_value: 0.0,
        density_label: 9.93,
        labels_consistent: 0.0,
        knots: 34,
        coefficients: &[
            -734.172, 1163.613, -135.157, -878.514, 791.787, -308.6236, 58.6764, -4.4494, 0.0,
        ],
    },
    Foam {
        slug: "polystyrene-60p07-nist-thermal-reference",
        family: "polystyrene",
        gas_axis: "source-cell-gas-known",
        gas_value: 0.0,
        density_label: 60.07,
        labels_consistent: 1.0,
        knots: 39,
        coefficients: &[
            2139.33, -6518.015, 9650.919, -9229.889, 6104.0571, -2739.14, 784.878, -128.40103,
            9.08405,
        ],
    },
];

fn semantic(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}
fn cp_key() -> PropertyKey {
    PropertyKey::new("specific-heat-capacity", Dims([2, 0, -2, -1, 0, 0]))
}
fn point(foam: &Foam, t: f64) -> QueryPoint {
    let mut p = QueryPoint::new();
    for (name, value) in [
        ("engineering-reference-model", 1.0),
        ("source-grade-known", 0.0),
        ("source-processing-known", 0.0),
        ("source-moisture-condition-known", 0.0),
        ("source-pressure-known", 0.0),
        ("source-specimens-matched", 0.0),
        ("source-orientation-known", 0.0),
        ("source-cell-gas-composition-known", 0.0),
        (
            "source-density-unit-labels-consistent",
            foam.labels_consistent,
        ),
        (foam.gas_axis, foam.gas_value),
    ] {
        p = p
            .with_quantity(name, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    p.with_quantity(
        format!("source-family-nist-{}", foam.family),
        QuantitySpec::dimensional(Dims::NONE),
        1.0,
    )
    .unwrap()
    .with_quantity(
        "source-reference-density",
        QuantitySpec::dimensional(Dims([-3, 1, 0, 0, 0, 0])),
        foam.density_label,
    )
    .unwrap()
    .with_quantity(
        "temperature",
        semantic(QuantityKind::AbsoluteTemperature),
        t,
    )
    .unwrap()
}
fn compile(foam: &Foam) -> NormalizedPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let output = std::env::temp_dir().join(format!(
        "nist-foam-capacity-{}-{nonce}-{}.fsmatpk",
        std::process::id(),
        foam.slug
    ));
    let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(root.join(format!("data/matdb/seed-v1/{}/manifest.tsv", foam.slug)))
        .arg("--out")
        .arg(&output)
        .current_dir(root)
        .env("CARGO_WORKSPACE_DIR", root)
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "{}: {}",
        foam.slug,
        String::from_utf8_lossy(&run.stderr)
    );
    let bytes = std::fs::read(output).unwrap();
    let raw = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(raw.content_hash(), &bytes).unwrap();
    assert_eq!(pack.to_bytes(), bytes);
    assert_eq!(pack.claims().claims_ordered().count(), 2);
    pack
}
fn card(foam: &Foam, pack: &NormalizedPack) -> MaterialCard {
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: foam.slug.into(),
            phase: "solid".into(),
            process: "NIST source reference; supplied mass; no density or conductivity claim"
                .into(),
            revision: 0,
        },
        pack.claims().clone(),
        vec![],
    )
    .unwrap()
}
fn original(foam: &Foam, t: f64) -> f64 {
    let x = t.log10();
    10.0_f64.powf(foam.coefficients.iter().rev().fold(0.0, |y, a| y * x + a))
}
fn close(a: f64, b: f64, relative: f64, absolute: f64) {
    assert!(
        (a - b).abs() <= relative * b.abs() + absolute,
        "{a:.16e} != {b:.16e}"
    );
}
fn original_integral(foam: &Foam, lo: f64, hi: f64, entropy: bool) -> f64 {
    let value = |t| original(foam, t) / if entropy { t } else { 1.0 };
    let n = 2048;
    let step = (hi - lo) / f64::from(n);
    let mut sum = value(lo) + value(hi);
    for i in 1..n {
        sum += if i % 2 == 0 { 2.0 } else { 4.0 } * value(lo + f64::from(i) * step);
    }
    sum * step / 3.0
}

#[test]
fn g1_g3_literal_capacities_units_receipts_and_conflicting_density_label_survive() {
    let mut count = 0;
    for foam in FOAMS {
        let pack = compile(foam);
        let line = pack
            .claims()
            .query_line_typed(
                &cp_key(),
                "temperature",
                &point(foam, 100.0),
                &point(foam, 300.0),
                SINGLE,
            )
            .unwrap();
        assert_eq!(line.knots.len(), foam.knots);
        count += line.knots.len() + 1;
        for (t, answer) in &line.knots {
            pack.claims().verify_receipt(&answer.receipt).unwrap();
            assert_eq!(
                answer.evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
            close(answer.evidence.value.value, original(foam, *t), 2e-7, 1e-10);
        }
        for pair in line.knots.windows(2) {
            for fraction in [0.25, 0.5, 0.75] {
                let t = pair[0].0 + fraction * (pair[1].0 - pair[0].0);
                let answer = pack
                    .claims()
                    .query_typed(&cp_key(), &point(foam, t), SelectionPolicy::SingleClaimOnly)
                    .unwrap();
                close(
                    answer.evidence.value.value,
                    original(foam, t),
                    8.1e-5,
                    1e-10,
                );
            }
        }
        let density_point = point(foam, 200.0)
            .with_quantity(
                "source-density-temperature-known",
                QuantitySpec::dimensional(Dims::NONE),
                0.0,
            )
            .unwrap();
        let density = pack
            .claims()
            .query_typed(
                &PropertyKey::new("nominal-source-density", Dims([-3, 1, 0, 0, 0, 0])),
                &density_point,
                SelectionPolicy::SingleClaimOnly,
            )
            .unwrap();
        assert_eq!(density.evidence.value.value, foam.density_label);
        pack.claims().verify_receipt(&density.receipt).unwrap();
    }
    assert_eq!(count, 104);
    // The source's "9.93 kg/m3 (=0.062 lb/ft3)" labels are inconsistent.
    // Neither is silently corrected or promoted into a physical density law.
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let source = std::fs::read_to_string(
        root.join("data/matdb/seed-v1/polystyrene-9p93-nist-thermal-reference/properties.tsv"),
    )
    .unwrap();
    assert!(source.contains("9.93 kg/m3 to 0.062 lb/ft3"));
    assert!(source.contains("source-density-unit-labels-consistent=0"));
}

#[test]
fn g0_unknown_conditions_and_distinct_foam_populations_cannot_supply_missing_laws() {
    for foam in FOAMS {
        let pack = compile(foam);
        let p = point(foam, 200.0);
        for axis in [
            "engineering-reference-model",
            "source-grade-known",
            "source-processing-known",
            "source-moisture-condition-known",
            "source-pressure-known",
            "source-specimens-matched",
            "source-orientation-known",
            "source-cell-gas-composition-known",
            "source-density-unit-labels-consistent",
            foam.gas_axis,
        ] {
            let wrong = p.clone().with(axis, 1.0 - p.axes()[axis]).unwrap();
            assert!(
                pack.claims()
                    .query_typed(&cp_key(), &wrong, SelectionPolicy::SingleClaimOnly)
                    .is_err()
            );
        }
        assert!(
            pack.claims()
                .query_typed(
                    &cp_key(),
                    &p.clone()
                        .with("source-reference-density", foam.density_label + 0.01)
                        .unwrap(),
                    SelectionPolicy::SingleClaimOnly
                )
                .is_err()
        );
        for t in [100.0_f64.next_down(), 300.0_f64.next_up()] {
            assert!(
                pack.claims()
                    .query_typed(&cp_key(), &point(foam, t), SelectionPolicy::SingleClaimOnly)
                    .is_err()
            );
        }
        let wrong_temperature = p
            .with_quantity(
                "temperature",
                semantic(QuantityKind::TemperatureDifference),
                200.0,
            )
            .unwrap();
        assert!(
            pack.claims()
                .query_typed(
                    &cp_key(),
                    &wrong_temperature,
                    SelectionPolicy::SingleClaimOnly
                )
                .is_err()
        );
        for missing in [
            "density",
            "thermal-conductivity",
            "volumetric-heat-capacity",
            "linear-thermal-expansion-coefficient",
        ] {
            assert!(pack.claims().claims_for(missing).is_empty());
        }
    }
}

#[test]
fn g1_g3_foam_mass_heating_and_common_temperature_storage_need_no_density_or_k() {
    let packs: Vec<_> = FOAMS.iter().map(compile).collect();
    let cards: Vec<_> = FOAMS.iter().zip(&packs).map(|(f, p)| card(f, p)).collect();
    let bounds: Vec<_> = FOAMS
        .iter()
        .map(|f| [point(f, 100.0), point(f, 300.0)])
        .collect();
    let capacity = cp_key();
    let query = |i: usize, mass: f64| HeatStorageComponentQuery {
        id: (i + 1) as u64,
        card: &cards[i],
        capacity: &capacity,
        amount: HeatStorageAmount::MassKg(mass),
        temperature_axis: "temperature",
        lower: &bounds[i][0],
        upper: &bounds[i][1],
        selection: SINGLE,
    };
    let mut heated = Vec::new();
    for (i, foam) in FOAMS.iter().enumerate() {
        let storage = PreparedHeatStorage::try_new(&[query(i, 2.0)]).unwrap();
        assert_eq!(storage.temperature_span_k(), [100.0, 300.0]);
        close(
            storage.heat_capacity_j_k(200.0).unwrap(),
            2.0 * original(foam, 200.0),
            2e-7,
            1e-9,
        );
        let e = storage.energy_change_j(150.0, 250.0).unwrap();
        close(
            e,
            2.0 * original_integral(foam, 150.0, 250.0, false),
            8.1e-5,
            1e-7,
        );
        close(
            storage.energy_change_j(250.0, 150.0).unwrap(),
            -e,
            1e-13,
            1e-8,
        );
        close(
            storage.temperature_after_energy_change(150.0, e).unwrap(),
            250.0,
            0.0,
            2e-9,
        );
        close(
            storage.sensible_entropy_change_j_k(150.0, 250.0).unwrap(),
            2.0 * original_integral(foam, 150.0, 250.0, true),
            8.1e-5,
            1e-8,
        );
        assert!(storage.heat_capacity_j_k(99.0).is_err());
        assert!(storage.temperature_after_energy_change(200.0, 1e9).is_err());
        heated.push(
            storage
                .temperature_after_energy_change(200.0, 20000.0)
                .unwrap(),
        );
    }
    assert!(heated[0] < heated[1]);
    assert!((heated[1] - heated[2]).abs() > 0.01);

    // These are two separately sourced components held at a common lumped
    // temperature. Their heat capacities add; no homogenized k or rho is made.
    let mixture = PreparedHeatStorage::try_new(&[query(0, 1.0), query(2, 3.0)]).unwrap();
    assert_eq!(mixture.components().len(), 2);
    let total = original(&FOAMS[0], 200.0) + 3.0 * original(&FOAMS[2], 200.0);
    close(mixture.heat_capacity_j_k(200.0).unwrap(), total, 2e-7, 1e-8);
    let expected = original_integral(&FOAMS[0], 150.0, 250.0, false)
        + 3.0 * original_integral(&FOAMS[2], 150.0, 250.0, false);
    close(
        mixture.energy_change_j(150.0, 250.0).unwrap(),
        expected,
        8.1e-5,
        1e-7,
    );
    let energy = mixture.energy_change_j(175.0, 260.0).unwrap();
    close(
        mixture
            .temperature_after_energy_change(175.0, energy)
            .unwrap(),
        260.0,
        0.0,
        2e-9,
    );
    let substituted = PreparedHeatStorage::try_new(&[query(0, 1.0), query(1, 3.0)]).unwrap();
    assert_ne!(mixture.identity(), substituted.identity());
    assert!(
        (mixture
            .temperature_after_energy_change(200.0, 20000.0)
            .unwrap()
            - substituted
                .temperature_after_energy_change(200.0, 20000.0)
                .unwrap())
        .abs()
            > 0.01
    );
    let mut volume_query = query(0, 1.0);
    volume_query.amount = HeatStorageAmount::ReferenceVolumeM3(1.0);
    assert!(PreparedHeatStorage::try_new(&[volume_query]).is_err());
}
