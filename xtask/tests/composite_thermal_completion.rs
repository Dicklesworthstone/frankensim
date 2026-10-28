//! Source thermal reductions distinguish cured adhesive volume from carbon fiber constituents.
#![deny(unsafe_code)]

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use std::{path::Path, process::Command};

const RHO: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const CP: Dims = Dims([2, 0, -2, -1, 0, 0]);
const CAPACITY: Dims = Dims([-1, 1, -2, -1, 0, 0]);
const K: Dims = Dims([1, 1, -3, -1, 0, 0]);
const DIFFUSIVITY: Dims = Dims([2, 0, -1, 0, 0, 0]);
const PA: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const METRE: Dims = Dims([1, 0, 0, 0, 0, 0]);
const SECOND: Dims = Dims([0, 0, 1, 0, 0, 0]);
const WATT: Dims = Dims([2, 1, -3, 0, 0, 0]);
const OHM_M: Dims = Dims([3, 1, -3, 0, -2, 0]);
const INV_K: Dims = Dims([0, 0, 0, -1, 0, 0]);

fn compile(slug: &str, count: usize) -> NormalizedPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let path = std::env::temp_dir().join(format!(
        "frankensim-composite-thermal-{}-{nonce}-{slug}.fsmatpk",
        std::process::id()
    ));
    let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(root.join(format!("data/matdb/seed-v1/{slug}/manifest.tsv")))
        .arg("--out")
        .arg(&path)
        .env("CARGO_WORKSPACE_DIR", root)
        .current_dir(root)
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "{slug}: {}",
        String::from_utf8_lossy(&run.stderr)
    );
    let bytes = std::fs::read(path).unwrap();
    let pack = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(pack.content_hash(), &bytes).unwrap();
    assert_eq!(pack.to_bytes(), bytes);
    assert_eq!(pack.pack_id(), slug);
    assert_eq!(pack.claims().claims_ordered().count(), count);
    pack
}

fn axis(p: QueryPoint, name: &str, dims: Dims, value: f64) -> QueryPoint {
    p.with_quantity(name, QuantitySpec::dimensional(dims), value)
        .unwrap()
}

fn flags(values: &[(&str, f64)]) -> QueryPoint {
    values.iter().fold(QueryPoint::new(), |p, (name, value)| {
        axis(p, name, Dims::NONE, *value)
    })
}

fn temp(p: QueryPoint, t: f64) -> QueryPoint {
    p.with_quantity(
        "temperature",
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        )),
        t,
    )
    .unwrap()
}

fn key(name: &str, dims: Dims) -> PropertyKey {
    let quantity = if name.starts_with("thermal-conductivity") {
        assert_eq!(dims, K);
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::ThermalConductivity,
            ValueForm::Static,
        ))
    } else {
        QuantitySpec::dimensional(dims)
    };
    PropertyKey::with_quantity(name, quantity)
}

fn sample(pack: &NormalizedPack, name: &str, dims: Dims, point: &QueryPoint, expected: f64) -> f64 {
    let answer = pack
        .claims()
        .query_typed(&key(name, dims), point, SelectionPolicy::SingleClaimOnly)
        .unwrap_or_else(|e| panic!("{} {name}: {e:?}", pack.pack_id()));
    let value = answer.evidence.value.value;
    assert!((value - expected).abs() <= 1e-12 * expected.abs().max(1e-10));
    assert_eq!(answer.evidence.value.dims, dims);
    assert_eq!(
        answer.evidence.value.uncertainty,
        UncertaintyModel::Unstated
    );
    assert!(answer.receipt.observation_backed);
    assert!(!answer.receipt.source_hashes.is_empty());
    pack.claims().verify_receipt(&answer.receipt).unwrap();
    value
}

fn refuses(pack: &NormalizedPack, name: &str, dims: Dims, p: &QueryPoint) {
    assert!(
        pack.claims()
            .query_typed(&key(name, dims), p, SelectionPolicy::SingleClaimOnly)
            .is_err(),
        "{} must refuse {name} at {p:?}",
        pack.pack_id()
    );
}

fn adhesive(grade: &str) -> QueryPoint {
    let loctite = grade == "loctite-aa-h4500";
    let p = flags(&[
        (&format!("source-grade-{grade}"), 1.0),
        ("source-cured-material", 1.0),
        ("source-room-temperature-cure", 1.0),
        ("source-cure-temperature-known", 0.0),
        ("source-mix-ratio-known", 0.0),
        ("source-moisture-condition-known", 0.0),
        ("source-pressure-known", 0.0),
        ("source-tps-isotropic-reduction", 1.0),
        ("source-independent-thermal-measurements-known", 0.0),
        ("source-repeat-count", 5.0),
        (
            "source-test-duration-known",
            if loctite { 1.0 } else { 0.0 },
        ),
    ]);
    let p = axis(
        axis(
            axis(p, "cure-time", SECOND, 86400.0),
            "source-sensor-radius",
            METRE,
            0.0032,
        ),
        "source-test-power",
        WATT,
        if loctite { 0.1 } else { 0.2 },
    );
    let p = if loctite {
        axis(p, "source-test-duration", SECOND, 5.0)
    } else {
        p
    };
    temp(p, 296.15)
}

#[test]
fn g0_g3_tps_outputs_retain_volumetric_capacity_and_common_measurement_context() {
    for (grade, k, diffusivity, capacity) in [
        ("loctite-aa-h4500", 0.200, 0.117e-6, 1.717e6),
        ("prima-solder-eg8020", 1.905, 1.849e-6, 1.030e6),
    ] {
        let pack = compile(&format!("{grade}-thermtest-23c"), 3);
        let p = adhesive(grade);
        let k = sample(&pack, "thermal-conductivity", K, &p, k);
        let reported = sample(&pack, "thermal-diffusivity", DIFFUSIVITY, &p, diffusivity);
        let capacity = sample(&pack, "volumetric-heat-capacity", CAPACITY, &p, capacity);
        // All outputs share one TPS reduction, not independent measurements.
        // Rounded published numbers remain intact, rather than enforcing closure.
        let derived = k / capacity;
        assert_ne!(derived, reported);
        assert!((derived / reported - 1.0).abs() < 0.005);
        refuses(&pack, "density", RHO, &p);
        refuses(&pack, "specific-heat-capacity", CP, &p);
        for t in [293.15, 296.14, 296.16, 300.0] {
            refuses(&pack, "thermal-conductivity", K, &temp(p.clone(), t));
        }
        for (name, value) in [
            ("source-cured-material", 0.0),
            ("source-cure-temperature-known", 1.0),
            ("source-mix-ratio-known", 1.0),
            ("source-independent-thermal-measurements-known", 1.0),
            ("source-tps-isotropic-reduction", 0.0),
            ("source-repeat-count", 1.0),
        ] {
            refuses(
                &pack,
                "volumetric-heat-capacity",
                CAPACITY,
                &axis(p.clone(), name, Dims::NONE, value),
            );
        }
        refuses(
            &pack,
            "thermal-conductivity",
            K,
            &axis(p.clone(), "cure-time", SECOND, 3600.0),
        );
        refuses(
            &pack,
            "thermal-conductivity",
            K,
            &axis(p.clone(), "source-sensor-radius", METRE, 0.0064),
        );
        refuses(
            &pack,
            "thermal-conductivity",
            K,
            &axis(p.clone(), "temperature", Dims([0, 0, 0, 1, 0, 0]), 296.15),
        );
        if grade == "prima-solder-eg8020" {
            refuses(
                &pack,
                "thermal-conductivity",
                K,
                &axis(p, "source-test-duration-known", Dims::NONE, 1.0),
            );
        } else {
            refuses(
                &pack,
                "thermal-conductivity",
                K,
                &axis(p, "source-test-duration", SECOND, 2.0),
            );
        }
    }
}

fn fiber(grade: &str) -> QueryPoint {
    flags(&[
        (&format!("source-grade-hextow-{grade}"), 1.0),
        ("material-form-fiber", 1.0),
        ("engineering-reference-model", 1.0),
        ("source-test-temperature-known", 0.0),
        ("source-moisture-condition-known", 0.0),
        ("source-pressure-known", 0.0),
        ("source-sizing-known", 0.0),
        ("source-tow-count-known", 0.0),
        ("source-specimens-matched", 0.0),
        ("thermochemical-calorie-conversion", 1.0),
        ("source-measurement-direction-known", 0.0),
        ("source-expansion-interval-known", 0.0),
        ("source-test-rate-known", 0.0),
    ])
}

#[test]
fn g0_g3_fiber_thermal_mechanical_and_electrical_values_preserve_constituent_roles() {
    for (grade, rho, cal, k, modulus, strength, diameter, resistivity, expansion) in [
        (
            "as4", 1790.0, 0.27, 6.83, 231e9, 4501e6, 7.1e-6, 1.7e-5, -0.63e-6,
        ),
        (
            "im7", 1780.0, 0.21, 5.40, 276e9, 5670e6, 5.2e-6, 1.5e-5, -0.64e-6,
        ),
    ] {
        let pack = compile(&format!("hexcel-hextow-{grade}-fiber-2023"), 9);
        let p = fiber(grade);
        sample(&pack, "density", RHO, &p, rho);
        sample(&pack, "specific-heat-capacity", CP, &p, cal * 4.184 / 0.001);
        sample(&pack, "fiber-tensile-modulus", PA, &p, modulus);
        sample(&pack, "fiber-tensile-strength", PA, &p, strength);
        sample(&pack, "fiber-tensile-elongation", Dims::NONE, &p, 0.018);
        sample(&pack, "filament-diameter", METRE, &p, diameter);
        sample(
            &pack,
            "thermal-conductivity-direction-unspecified",
            K,
            &p,
            k,
        );
        sample(
            &pack,
            "electrical-resistivity-direction-unspecified",
            OHM_M,
            &p,
            resistivity,
        );
        sample(
            &pack,
            "reported-linear-expansion-coefficient-direction-unspecified",
            INV_K,
            &p,
            expansion,
        );
        // The fiber and laminate tables on the same sheet are separate materials.
        for (name, dims) in [
            ("young-modulus", PA),
            ("poisson-ratio", Dims::NONE),
            ("thermal-conductivity", K),
            ("electrical-resistivity", OHM_M),
            ("linear-thermal-expansion", INV_K),
        ] {
            refuses(&pack, name, dims, &p);
        }
        for (name, value) in [
            ("material-form-fiber", 0.0),
            ("engineering-reference-model", 0.0),
            ("source-test-temperature-known", 1.0),
            ("source-specimens-matched", 1.0),
            ("source-sizing-known", 1.0),
        ] {
            refuses(
                &pack,
                "density",
                RHO,
                &axis(p.clone(), name, Dims::NONE, value),
            );
        }
        refuses(
            &pack,
            "specific-heat-capacity",
            CP,
            &axis(
                p.clone(),
                "thermochemical-calorie-conversion",
                Dims::NONE,
                0.0,
            ),
        );
        refuses(
            &pack,
            "thermal-conductivity-direction-unspecified",
            K,
            &axis(p, "source-measurement-direction-known", Dims::NONE, 1.0),
        );
        // There is deliberately no authored or measured temperature trajectory.
        assert!(
            pack.claims()
                .claims_ordered()
                .all(|(_, claim)| !claim.validity.bounds().contains_key("temperature"))
        );
    }
}

#[test]
fn g3_primary_source_discrepancies_and_identity_limits_remain_in_the_pack() {
    let loctite = compile("loctite-aa-h4500-thermtest-23c", 3);
    let silver = compile("prima-solder-eg8020-thermtest-23c", 3);
    // Compiler retains source citations on every immutable claim, not just an index.
    for pack in [&loctite, &silver] {
        assert!(pack.claims().claims_ordered().all(|(_, claim)| {
            claim
                .provenance
                .source
                .contains("thermtest.com/application/thermal-conductivity-of-epoxy")
        }));
    }
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let original = |slug: &str| {
        std::fs::read_to_string(root.join(format!("data/matdb/seed-v1/{slug}/properties.tsv")))
            .unwrap()
    };
    assert!(original("loctite-aa-h4500-thermtest-23c").contains("methacrylate/acrylic"));
    let text = original("prima-solder-eg8020-thermtest-23c");
    assert!(text.contains("2 s measurement") && text.contains("5 s at 200 mW"));
    assert_ne!(loctite.content_hash(), silver.content_hash());
    refuses(
        &silver,
        "thermal-conductivity",
        K,
        &adhesive("loctite-aa-h4500"),
    );
    refuses(
        &loctite,
        "thermal-conductivity",
        K,
        &adhesive("prima-solder-eg8020"),
    );
}
