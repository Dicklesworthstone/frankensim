//! G0/G3: sourced construction properties retain conditions and measurement roles.
//! The explicit local-k references are engineering models, not experimental claims.
#![deny(unsafe_code)]

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use std::{path::Path, process::Command};

const RHO: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const PA: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const CP: Dims = Dims([2, 0, -2, -1, 0, 0]);
const K: Dims = Dims([1, 1, -3, -1, 0, 0]);
const FORCE: Dims = Dims([1, 1, -2, 0, 0, 0]);

fn compile(slug: &str, count: usize) -> NormalizedPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "frankensim-construction-thermal-{}-{nonce}-{slug}",
        std::process::id()
    ));
    std::fs::create_dir(&dir).unwrap();
    let path = dir.join("material.fsmatpk");
    let output = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(root.join(format!("data/matdb/seed-v1/{slug}/manifest.tsv")))
        .arg("--out")
        .arg(&path)
        .env("CARGO_WORKSPACE_DIR", root)
        .current_dir(root)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{slug}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let bytes = std::fs::read(path).unwrap();
    let pack = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(pack.content_hash(), &bytes).unwrap();
    assert_eq!(pack.to_bytes(), bytes);
    assert_eq!(pack.pack_id(), slug);
    assert_eq!(pack.claims().claims_ordered().count(), count);
    pack
}

fn flags(p: QueryPoint, values: &[(&str, f64)]) -> QueryPoint {
    values.iter().fold(p, |p, (axis, value)| {
        p.with_quantity(*axis, QuantitySpec::dimensional(Dims::NONE), *value)
            .unwrap()
    })
}

fn base(slug: &str) -> QueryPoint {
    flags(
        QueryPoint::new(),
        &[
            (&format!("source-grade-{slug}"), 1.0),
            ("source-moisture-condition-known", 0.0),
            ("source-gas-pressure-known", 0.0),
            ("source-specimens-matched", 0.0),
        ],
    )
}

fn unknown(slug: &str) -> QueryPoint {
    flags(base(slug), &[("source-temperature-known", 0.0)])
}

fn temperature(p: QueryPoint, axis: &str, kelvin: f64) -> QueryPoint {
    p.with_quantity(
        axis,
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        )),
        kelvin,
    )
    .unwrap()
}

fn key(name: &str, dims: Dims) -> PropertyKey {
    let spec = if matches!(
        name,
        "thermal-conductivity" | "apparent-thermal-conductivity" | "declared-thermal-conductivity"
    ) {
        assert_eq!(dims, K);
        QuantitySpec::semantic(SemanticType::new(
            QuantityKind::ThermalConductivity,
            ValueForm::Static,
        ))
    } else {
        QuantitySpec::dimensional(dims)
    };
    PropertyKey::with_quantity(name, spec)
}

fn sample(pack: &NormalizedPack, name: &str, dims: Dims, p: &QueryPoint, expected: f64) {
    let result = pack
        .claims()
        .query_typed(&key(name, dims), p, SelectionPolicy::SingleClaimOnly)
        .unwrap_or_else(|e| panic!("{} {name}: {e:?}", pack.pack_id()));
    assert!((result.evidence.value.value - expected).abs() <= 1e-11 * expected.abs().max(1e-6));
    assert_eq!(result.evidence.value.dims, dims);
    assert!(matches!(
        result.evidence.value.uncertainty,
        UncertaintyModel::Unstated
    ));
    assert!(result.receipt.observation_backed);
    assert!(!result.receipt.source_hashes.is_empty());
    pack.claims().verify_receipt(&result.receipt).unwrap();
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

fn transport(slug: &str) -> QueryPoint {
    flags(
        base(slug),
        &[
            ("source-thermal-direction-known", 0.0),
            ("source-test-temperature-difference-known", 0.0),
            ("source-test-thickness-known", 0.0),
        ],
    )
}

#[test]
fn g0_g3_foamglas_preserves_mean_temperature_polynomial_samples_and_strength_limits() {
    let slug = "foamglas-one-astm-2022";
    let pack = compile(slug, 9);
    let p = unknown(slug);
    sample(&pack, "density", RHO, &p, 115.0);
    sample(
        &pack,
        "density-relative-manufacturing-tolerance",
        Dims::NONE,
        &p,
        0.10,
    );
    let mechanical = flags(
        p.clone(),
        &[
            ("source-measurement-direction-known", 0.0),
            ("source-test-rate-known", 0.0),
        ],
    );
    sample(
        &pack,
        "compression-strength-average",
        PA,
        &mechanical,
        620_000.0,
    );
    sample(
        &pack,
        "compression-strength-lower-specification-limit",
        PA,
        &mechanical,
        414_000.0,
    );
    sample(
        &pack,
        "flexural-strength-lower-specification-limit",
        PA,
        &mechanical,
        283_000.0,
    );
    let modulus = flags(
        p.clone(),
        &[
            ("source-measurement-direction-known", 0.0),
            ("source-assumed-poisson-ratio", 0.25),
            ("source-frequency-known", 0.0),
        ],
    );
    sample(&pack, "elastic-modulus-approximation", PA, &modulus, 900e6);
    refuses(
        &pack,
        "elastic-modulus-approximation",
        PA,
        &flags(modulus.clone(), &[("source-assumed-poisson-ratio", 0.30)]),
    );
    refuses(&pack, "young-modulus", PA, &modulus);
    refuses(&pack, "poisson-ratio", Dims::NONE, &modulus);
    sample(
        &pack,
        "specific-heat-capacity",
        CP,
        &temperature(base(slug), "temperature", 298.15),
        770.0,
    );
    refuses(
        &pack,
        "specific-heat-capacity",
        CP,
        &temperature(base(slug), "temperature", 297.15),
    );

    // Product polynomial evaluations, not the neighboring ASTM C552 limit row.
    for (celsius, conductivity) in [
        (-165.0, 0.020),
        (-157.0, 0.021),
        (-129.0, 0.023),
        (-101.0, 0.026),
        (-73.0, 0.029),
        (-46.0, 0.032),
        (-18.0, 0.036),
        (10.0, 0.040),
        (24.0, 0.042),
        (38.0, 0.044),
        (93.0, 0.054),
        (149.0, 0.066),
        (204.0, 0.078),
    ] {
        sample(
            &pack,
            "apparent-thermal-conductivity",
            K,
            &temperature(transport(slug), "mean-test-temperature", celsius + 273.15),
            conductivity,
        );
    }
    refuses(
        &pack,
        "apparent-thermal-conductivity",
        K,
        &temperature(transport(slug), "mean-test-temperature", 304.15),
    );
    let local = temperature(
        flags(
            transport(slug),
            &[
                ("engineering-reference-model", 1.0),
                ("mean-test-temperature-as-local-model", 1.0),
            ],
        ),
        "temperature",
        304.15,
    );
    sample(&pack, "thermal-conductivity", K, &local, 0.043);
    refuses(
        &pack,
        "thermal-conductivity",
        K,
        &flags(
            local.clone(),
            &[("mean-test-temperature-as-local-model", 0.0)],
        ),
    );
    refuses(
        &pack,
        "thermal-conductivity",
        K,
        &temperature(local, "temperature", 478.15),
    );
    refuses(
        &pack,
        "density",
        RHO,
        &flags(p, &[("source-moisture-condition-known", 1.0)]),
    );
}

#[test]
fn g0_g3_calcium_silicate_boards_keep_heat_capacity_and_directional_domains() {
    let l = "promatect-l-2017";
    let pack = compile(l, 7);
    let p = unknown(l);
    sample(&pack, "density", RHO, &p, 450.0);
    sample(&pack, "specific-heat-capacity", CP, &p, 950.0);
    sample(
        &pack,
        "water-vapor-diffusion-resistance-factor",
        Dims::NONE,
        &p,
        3.2,
    );
    let mechanical = flags(
        p.clone(),
        &[
            ("source-measurement-direction-known", 0.0),
            ("source-test-rate-known", 0.0),
        ],
    );
    sample(&pack, "cold-compression-strength", PA, &mechanical, 2.4e6);
    sample(&pack, "bending-strength", PA, &mechanical, 3.1e6);
    for (celsius, k) in [(20.0, 0.08), (200.0, 0.10), (400.0, 0.11)] {
        sample(
            &pack,
            "apparent-thermal-conductivity",
            K,
            &temperature(transport(l), "mean-test-temperature", celsius + 273.15),
            k,
        );
    }
    refuses(
        &pack,
        "apparent-thermal-conductivity",
        K,
        &temperature(transport(l), "mean-test-temperature", 383.15),
    );
    let local = temperature(
        flags(
            transport(l),
            &[
                ("engineering-reference-model", 1.0),
                ("mean-test-temperature-as-local-model", 1.0),
            ],
        ),
        "temperature",
        383.15,
    );
    sample(&pack, "thermal-conductivity", K, &local, 0.09);
    refuses(
        &pack,
        "thermal-conductivity",
        K,
        &flags(local, &[("engineering-reference-model", 0.0)]),
    );

    let h = "promatect-h-2017";
    let pack = compile(h, 12);
    let p = unknown(h);
    sample(&pack, "density", RHO, &p, 870.0);
    sample(
        &pack,
        "water-vapor-diffusion-resistance-factor",
        Dims::NONE,
        &p,
        20.0,
    );
    sample(
        &pack,
        "specific-heat-capacity",
        CP,
        &temperature(base(h), "temperature", 673.15),
        920.0,
    );
    refuses(
        &pack,
        "specific-heat-capacity",
        CP,
        &temperature(base(h), "temperature", 293.15),
    );
    for (longitudinal, modulus, bending, tension) in
        [(1.0, 4200e6, 7.6e6, 4.8e6), (0.0, 2900e6, 4.8e6, 2.6e6)]
    {
        let mech = flags(
            p.clone(),
            &[
                ("load-along-board-longitudinal", longitudinal),
                ("source-test-rate-known", 0.0),
            ],
        );
        sample(&pack, "board-elastic-modulus", PA, &mech, modulus);
        sample(&pack, "bending-strength", PA, &mech, bending);
        sample(&pack, "tensile-strength", PA, &mech, tension);
        refuses(&pack, "young-modulus", PA, &mech);
    }
    sample(
        &pack,
        "cold-compression-strength",
        PA,
        &flags(
            p,
            &[
                ("source-measurement-direction-known", 0.0),
                ("source-test-rate-known", 0.0),
            ],
        ),
        9.3e6,
    );
    let raw = flags(
        transport(h),
        &[("source-temperature-definition-known", 0.0)],
    );
    for (celsius, k) in [(20.0, 0.17), (100.0, 0.19), (200.0, 0.21)] {
        sample(
            &pack,
            "apparent-thermal-conductivity",
            K,
            &temperature(raw.clone(), "source-test-temperature", celsius + 273.15),
            k,
        );
    }
    let local = temperature(
        flags(
            raw,
            &[
                ("engineering-reference-model", 1.0),
                ("source-test-temperature-as-local-model", 1.0),
            ],
        ),
        "temperature",
        333.15,
    );
    sample(&pack, "thermal-conductivity", K, &local, 0.18);
    refuses(
        &pack,
        "thermal-conductivity",
        K,
        &temperature(local.clone(), "temperature", 673.15),
    );
    refuses(
        &pack,
        "thermal-conductivity",
        K,
        &flags(local, &[("source-temperature-definition-known", 1.0)]),
    );
}

#[test]
fn g0_g3_wood_fibre_batt_preserves_declared_values_and_acoustic_resistivity_units() {
    let slug = "steico-flex036-v5-2026";
    let pack = compile(slug, 5);
    let p = unknown(slug);
    sample(&pack, "density", RHO, &p, 50.0);
    sample(&pack, "specific-heat-capacity", CP, &p, 2100.0);
    sample(
        &pack,
        "water-vapor-diffusion-resistance-factor",
        Dims::NONE,
        &p,
        2.0,
    );
    let declared = flags(
        p.clone(),
        &[
            ("source-is-declared-value", 1.0),
            ("source-measurement-direction-known", 0.0),
        ],
    );
    sample(&pack, "declared-thermal-conductivity", K, &declared, 0.036);
    sample(
        &pack,
        "airflow-resistivity-lower-limit",
        Dims([-3, 1, -1, 0, 0, 0]),
        &declared,
        5000.0,
    );
    refuses(&pack, "thermal-conductivity", K, &declared);
    refuses(
        &pack,
        "declared-thermal-conductivity",
        K,
        &flags(declared, &[("source-is-declared-value", 0.0)]),
    );
    refuses(
        &pack,
        "specific-heat-capacity",
        CP,
        &flags(p.clone(), &[("source-temperature-known", 1.0)]),
    );
    refuses(
        &pack,
        "density",
        RHO,
        &flags(p, &[(&format!("source-grade-{slug}"), 0.0)]),
    );
}

#[test]
fn g0_g3_acetylated_wood_keeps_conditioning_finite_shrinkage_and_janka_force() {
    let slug = "accoya-radiata-pine-datasheet-2025";
    let pack = compile(slug, 13);
    let b = flags(
        QueryPoint::new(),
        &[
            (&format!("source-grade-{slug}"), 1.0),
            ("source-strength-grade-known", 0.0),
            ("source-gas-pressure-known", 0.0),
            ("source-specimens-matched", 0.0),
        ],
    );
    let conditioned = temperature(
        flags(b.clone(), &[("conditioning-relative-humidity", 0.65)]),
        "temperature",
        293.15,
    );
    for (name, dims, value) in [
        ("density", RHO, 515.0),
        ("density-reported-range-low", RHO, 400.0),
        ("density-reported-range-high", RHO, 600.0),
        ("equilibrium-moisture-content-low", Dims::NONE, 0.03),
        ("equilibrium-moisture-content-high", Dims::NONE, 0.05),
    ] {
        sample(&pack, name, dims, &conditioned, value);
    }
    refuses(
        &pack,
        "density",
        RHO,
        &flags(
            conditioned.clone(),
            &[("conditioning-relative-humidity", 0.05)],
        ),
    );
    refuses(
        &pack,
        "density",
        RHO,
        &temperature(conditioned.clone(), "temperature", 298.15),
    );
    refuses(
        &pack,
        "density",
        RHO,
        &flags(conditioned.clone(), &[("source-strength-grade-known", 1.0)]),
    );
    for (radial, conditioned_shrinkage, dry_shrinkage) in [(1.0, 0.004, 0.007), (0.0, 0.008, 0.015)]
    {
        let wet_to_conditioned = flags(
            conditioned.clone(),
            &[
                ("radial-direction", radial),
                ("initial-state-wet", 1.0),
                ("final-state-oven-dry", 0.0),
            ],
        );
        sample(
            &pack,
            "linear-shrinkage",
            Dims::NONE,
            &wet_to_conditioned,
            conditioned_shrinkage,
        );
        let wet_to_dry = flags(
            b.clone(),
            &[
                ("source-temperature-known", 0.0),
                ("radial-direction", radial),
                ("initial-state-wet", 1.0),
                ("final-state-oven-dry", 1.0),
            ],
        );
        sample(
            &pack,
            "linear-shrinkage",
            Dims::NONE,
            &wet_to_dry,
            dry_shrinkage,
        );
        refuses(
            &pack,
            "linear-shrinkage",
            Dims::NONE,
            &flags(wet_to_dry, &[("initial-state-wet", 0.0)]),
        );
    }
    let unconditioned = flags(
        b,
        &[
            ("source-temperature-known", 0.0),
            ("source-moisture-condition-known", 0.0),
            ("source-measurement-direction-known", 0.0),
            ("source-test-rate-known", 0.0),
        ],
    );
    sample(
        &pack,
        "apparent-thermal-conductivity",
        K,
        &unconditioned,
        0.12,
    );
    sample(&pack, "bending-moe", PA, &unconditioned, 8.8e9);
    for (end, force) in [(0.0, 4100.0), (1.0, 6600.0)] {
        let p = flags(unconditioned.clone(), &[("load-on-end-grain", end)]);
        sample(&pack, "janka-hardness-force", FORCE, &p, force);
        refuses(&pack, "janka-hardness-force", PA, &p);
    }
    refuses(&pack, "young-modulus", PA, &unconditioned);
    refuses(&pack, "specific-heat-capacity", CP, &unconditioned);
}
