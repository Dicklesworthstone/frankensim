//! SAC305 secondary creep measurements and versioned solder ribbon references.
//! Exercise actual pack compilation, typed queries, receipts and domain refusals.
#![deny(unsafe_code)]

use fs_matdb::{
    ClaimSelection, InterpolationPolicy, NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy,
    UncertaintyModel,
};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use std::{path::Path, process::Command};

const PA: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const FORCE: Dims = Dims([1, 1, -2, 0, 0, 0]);
const RATE: Dims = Dims([0, 0, -1, 0, 0, 0]);
const RHO: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const K: Dims = Dims([1, 1, -3, -1, 0, 0]);
const TEMP: Dims = Dims([0, 0, 0, 1, 0, 0]);
const SPECIFIC_ENERGY: Dims = Dims([2, 0, -2, 0, 0, 0]);
const CP: Dims = Dims([2, 0, -2, -1, 0, 0]);

fn compile(slug: &str, count: usize) -> NormalizedPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "frankensim-solder-expansion-{}-{nonce}-{slug}",
        std::process::id()
    ));
    std::fs::create_dir(&dir).unwrap();
    let output = dir.join("material.fsmatpk");
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
    assert_eq!(pack.claims().claims_ordered().count(), count, "{slug}");
    pack
}

fn axis(point: QueryPoint, name: &str, dims: Dims, value: f64) -> QueryPoint {
    point
        .with_quantity(name, QuantitySpec::dimensional(dims), value)
        .unwrap()
}

fn flags(values: &[(&str, f64)]) -> QueryPoint {
    values.iter().fold(QueryPoint::new(), |p, (name, value)| {
        axis(p, name, Dims::NONE, *value)
    })
}

fn temperature(point: QueryPoint, value: f64) -> QueryPoint {
    point
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            value,
        )
        .unwrap()
}

fn key(name: &str, dims: Dims) -> PropertyKey {
    let spec = match name {
        "thermal-conductivity" => QuantitySpec::semantic(SemanticType::new(
            QuantityKind::ThermalConductivity,
            ValueForm::Static,
        )),
        "solidus-temperature" | "liquidus-temperature" => QuantitySpec::semantic(
            SemanticType::new(QuantityKind::AbsoluteTemperature, ValueForm::Static),
        ),
        _ => QuantitySpec::dimensional(dims),
    };
    PropertyKey::with_quantity(name, spec)
}

fn sample(pack: &NormalizedPack, name: &str, dims: Dims, point: &QueryPoint, expected: f64) {
    let answer = pack
        .claims()
        .query_typed(&key(name, dims), point, SelectionPolicy::SingleClaimOnly)
        .unwrap_or_else(|error| panic!("{} {name}: {error:?}", pack.pack_id()));
    assert!((answer.evidence.value.value - expected).abs() <= 1e-12 * expected.abs());
    assert_eq!(answer.evidence.value.dims, dims);
    assert_eq!(
        answer.evidence.value.uncertainty,
        UncertaintyModel::Unstated
    );
    assert!(answer.receipt.observation_backed);
    assert!(!answer.receipt.source_hashes.is_empty());
    pack.claims().verify_receipt(&answer.receipt).unwrap();
}

fn refuses(pack: &NormalizedPack, name: &str, dims: Dims, point: &QueryPoint) {
    assert!(
        pack.claims()
            .query_typed(&key(name, dims), point, SelectionPolicy::SingleClaimOnly)
            .is_err(),
        "{} must refuse {name}",
        pack.pack_id()
    );
}

fn creep_point(mode: &str, specimen: usize, temp: f64, stress_mpa: f64, load: f64) -> QueryPoint {
    let p = flags(&[
        (&format!("source-cohort-sac305-hsu2013-{mode}"), 1.0),
        ("source-specimen-number", specimen as f64),
        ("secondary-creep-stage", 1.0),
        ("source-table-coordinates", 1.0),
        ("source-aging-history-known", 0.0),
        ("source-cooling-history-known", 0.0),
        ("source-pressure-known", 0.0),
    ]);
    temperature(
        axis(
            axis(p, "applied-load", FORCE, load),
            &format!("average-{mode}-stress"),
            PA,
            stress_mpa * 1e6,
        ),
        temp,
    )
}

#[test]
fn sac305_secondary_creep_preserves_exact_observations_without_inventing_a_law() {
    let temps = [393.0, 405.0, 423.0, 438.0, 393.0, 405.0, 423.0, 438.0];
    for (mode, stress, rates) in [
        (
            "tensile",
            [0.669, 1.145, 1.076, 1.266, 3.657, 3.057, 2.890, 3.057],
            [
                8.54e-8, 3.36e-7, 5.06e-7, 7.95e-7, 1.28e-6, 2.10e-6, 3.87e-6, 6.12e-6,
            ],
        ),
        (
            "shear",
            [0.977, 0.711, 1.021, 1.193, 7.424, 6.282, 6.074, 5.833],
            [
                5.50e-8, 6.29e-8, 1.54e-7, 2.95e-7, 3.52e-7, 4.85e-7, 8.69e-7, 1.43e-6,
            ],
        ),
    ] {
        let pack = compile(&format!("solder-sac305-{mode}-creep-hsu2013"), 8);
        let property = format!("secondary-{mode}-creep-strain-rate");
        for (_, claim) in pack.claims().claims_ordered() {
            assert_eq!(claim.interpolation, InterpolationPolicy::TabulatedOnly);
            assert!(claim.provenance.source.contains("10.1155/2013/838075"));
        }
        assert!(pack.claims().observation_ids().any(|id| {
            let notes = &pack.claims().observation(id).unwrap().caveats;
            notes.contains("405 K") && notes.contains("408 K") && notes.contains("7.53 N")
        }));
        for i in 0..8 {
            let load = if i < 4 { 7.35 } else { 29.4 };
            let p = creep_point(mode, i + 1, temps[i], stress[i], load);
            sample(&pack, &property, RATE, &p, rates[i]);
            refuses(
                &pack,
                &property,
                RATE,
                &temperature(p.clone(), temps[i] + 0.1),
            );
            refuses(
                &pack,
                &property,
                RATE,
                &axis(
                    p.clone(),
                    &format!("average-{mode}-stress"),
                    PA,
                    stress[i] * 1e6 + 1.0,
                ),
            );
            refuses(
                &pack,
                &property,
                RATE,
                &axis(p, "source-aging-history-known", Dims::NONE, 1.0),
            );
        }
        let first = creep_point(mode, 1, temps[0], stress[0], 7.35);
        let second = creep_point(mode, 2, temps[1], stress[1], 7.35);
        refuses(&pack, &property, RATE, &temperature(second, 408.0));
        refuses(
            &pack,
            &property,
            RATE,
            &axis(first.clone(), "applied-load", FORCE, 7.53),
        );
        refuses(&pack, &property, RATE, &QueryPoint::new());
        refuses(
            &pack,
            if mode == "tensile" {
                "secondary-shear-creep-strain-rate"
            } else {
                "secondary-tensile-creep-strain-rate"
            },
            RATE,
            &first,
        );
        // Two admissible samples do not establish continuous support between them.
        let last = creep_point(mode, 8, temps[7], stress[7], 29.4);
        assert!(
            pack.claims()
                .query_envelope_typed(
                    &key(&property, RATE),
                    &first,
                    &last,
                    ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
                )
                .is_err()
        );
    }
}

fn ribbon(alloy: &str, temperature_known: f64) -> QueryPoint {
    flags(&[
        ("engineering-reference-model", 1.0),
        ("source-specimens-matched", 0.0),
        ("source-test-temperature-known", temperature_known),
        ("source-test-method-known", 0.0),
        (&format!("source-alloy-{alloy}"), 1.0),
    ])
}

#[test]
fn ribbon_phase_density_and_strength_retain_exact_pdf_revision_and_reference_scope() {
    let psi_to_pa = 0.45359237 * 9.80665 / 0.0254_f64.powi(2);
    for (alloy, phase, density, tensile_psi) in [
        ("au80sn20", 553.15, 14510.0, 40000.0),
        ("au88ge12", 629.15, 14670.0, 26825.0),
        ("in97ag3", 416.15, 7380.0, 800.0),
    ] {
        let pack = compile(&format!("solder-{alloy}-indium-ribbon-r9"), 5);
        let p = ribbon(alloy, 0.0);
        sample(&pack, "solidus-temperature", TEMP, &p, phase);
        sample(&pack, "liquidus-temperature", TEMP, &p, phase);
        sample(&pack, "density", RHO, &p, density);
        sample(&pack, "tensile-strength", PA, &p, tensile_psi * psi_to_pa);
        assert!(pack.claims().claims_ordered().all(|(_, claim)| {
            claim
                .provenance
                .source
                .contains("Solder-Ribbon-and-Foil-PDS-97621-R9.pdf")
        }));
        refuses(&pack, "density", RHO, &ribbon(alloy, 1.0));
        refuses(
            &pack,
            "tensile-strength",
            PA,
            &axis(p.clone(), "source-specimens-matched", Dims::NONE, 1.0),
        );
        refuses(
            &pack,
            "density",
            RHO,
            &axis(p.clone(), "engineering-reference-model", Dims::NONE, 0.0),
        );
        refuses(
            &pack,
            "density",
            RHO,
            &axis(p.clone(), &format!("source-alloy-{alloy}"), Dims::NONE, 0.0),
        );
        // Eutectic endpoints cannot create an unmeasured caloric or bond model.
        refuses(&pack, "latent-heat-of-fusion", SPECIFIC_ENERGY, &p);
        refuses(&pack, "specific-heat-capacity", CP, &p);
        refuses(&pack, "bond-shear-strength", PA, &p);
    }
}

#[test]
fn ribbon_conductivity_is_qualified_only_at_the_printed_85c_temperature() {
    for (alloy, conductivity) in [("au80sn20", 57.0), ("au88ge12", 44.0), ("in97ag3", 73.0)] {
        let pack = compile(&format!("solder-{alloy}-indium-ribbon-r9"), 5);
        let p = temperature(ribbon(alloy, 1.0), 358.15);
        sample(&pack, "thermal-conductivity", K, &p, conductivity);
        refuses(&pack, "thermal-conductivity", K, &ribbon(alloy, 1.0));
        refuses(
            &pack,
            "thermal-conductivity",
            K,
            &temperature(ribbon(alloy, 0.0), 358.15),
        );
        for t in [293.15, 358.14, 358.16, 400.0] {
            refuses(&pack, "thermal-conductivity", K, &temperature(p.clone(), t));
        }
        // A typed absolute temperature is required; a same-dimension scalar is insufficient.
        refuses(
            &pack,
            "thermal-conductivity",
            K,
            &axis(ribbon(alloy, 1.0), "temperature", TEMP, 358.15),
        );
    }
}
