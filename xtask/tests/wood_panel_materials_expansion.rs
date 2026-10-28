//! Producer panel declarations retain lay-up, loading mode and statistical role.
#![deny(unsafe_code)]

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use std::{path::Path, process::Command};

const RHO: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const PA: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const K: Dims = Dims([1, 1, -3, -1, 0, 0]);
const CP: Dims = Dims([2, 0, -2, -1, 0, 0]);
const WISA: &str = "wisa-spruce-plywood-18mm-upm001-2026";
const LVL_R: &str = "steico-lvl-r-39mm-design-guide-2026";
const LVL_X: &str = "steico-lvl-x-39mm-design-guide-2026";

fn compile(slug: &str, count: usize) -> NormalizedPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "frankensim-wood-panel-{}-{nonce}-{slug}",
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
    let decoded = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap();
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

fn length(p: QueryPoint, axis: &str, value: f64) -> QueryPoint {
    p.with_quantity(
        axis,
        QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0])),
        value,
    )
    .unwrap()
}

fn base(slug: &str) -> QueryPoint {
    let wisa = slug == WISA;
    let p = flags(
        length(
            QueryPoint::new(),
            "nominal-thickness",
            if wisa { 0.018 } else { 0.039 },
        ),
        &[
            (&format!("source-grade-{slug}"), 1.0),
            ("ply-count", if wisa { 7.0 } else { 13.0 }),
            ("source-temperature-known", 0.0),
            ("source-moisture-condition-known", 0.0),
            ("source-pressure-known", 0.0),
            ("source-is-declared-value", 1.0),
            ("source-specimens-matched", 0.0),
        ],
    );
    if wisa {
        flags(p, &[("surface-is-coated", 0.0)])
    } else {
        flags(
            p,
            &[("cross-ply-count", if slug == LVL_X { 3.0 } else { 0.0 })],
        )
    }
}

fn loading(p: QueryPoint, parallel: f64, edge: Option<f64>) -> QueryPoint {
    let p = flags(p, &[("load-parallel-to-face-grain", parallel)]);
    edge.map_or_else(|| p.clone(), |v| flags(p.clone(), &[("load-edgewise", v)]))
}

fn key(name: &str, dims: Dims) -> PropertyKey {
    let spec = if matches!(name, "thermal-conductivity" | "design-thermal-conductivity") {
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

fn sample(pack: &NormalizedPack, name: &str, dims: Dims, p: &QueryPoint, expected: f64) -> f64 {
    let result = pack
        .claims()
        .query_typed(&key(name, dims), p, SelectionPolicy::SingleClaimOnly)
        .unwrap_or_else(|e| panic!("{} {name}: {e:?}", pack.pack_id()));
    let value = result.evidence.value.value;
    assert!((value - expected).abs() <= 1e-11 * expected.abs().max(1e-6));
    assert_eq!(result.evidence.value.dims, dims);
    assert!(matches!(
        result.evidence.value.uncertainty,
        UncertaintyModel::Unstated
    ));
    assert!(result.receipt.observation_backed);
    assert!(!result.receipt.source_hashes.is_empty());
    pack.claims().verify_receipt(&result.receipt).unwrap();
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

fn check_identity_and_unstated_conditions(pack: &NormalizedPack, p: &QueryPoint) {
    for (axis, value) in [
        ("source-temperature-known", 1.0),
        ("source-moisture-condition-known", 1.0),
        ("source-is-declared-value", 0.0),
        ("source-specimens-matched", 1.0),
        ("ply-count", 8.0),
    ] {
        refuses(pack, "density", RHO, &flags(p.clone(), &[(axis, value)]));
    }
    refuses(
        pack,
        "density",
        RHO,
        &length(p.clone(), "nominal-thickness", 0.027),
    );
    refuses(
        pack,
        "density",
        RHO,
        &flags(
            p.clone(),
            &[(&format!("source-grade-{}", pack.pack_id()), 0.0)],
        ),
    );
    refuses(pack, "specific-heat-capacity", CP, p);
    refuses(pack, "poisson-ratio", Dims::NONE, p);
}

fn check_transport(pack: &NormalizedPack, p: &QueryPoint, thermal_name: &str, wet: f64, dry: f64) {
    let thermal = flags(p.clone(), &[("source-conductivity-direction-known", 0.0)]);
    sample(pack, thermal_name, K, &thermal, 0.13);
    refuses(
        pack,
        thermal_name,
        K,
        &flags(thermal, &[("source-conductivity-direction-known", 1.0)]),
    );
    let vapor = flags(
        p.clone(),
        &[
            ("source-vapor-humidity-endpoints-known", 0.0),
            ("source-vapor-direction-known", 0.0),
        ],
    );
    for (state, value) in [(1.0, wet), (0.0, dry)] {
        sample(
            pack,
            "water-vapor-diffusion-resistance-factor",
            Dims::NONE,
            &flags(vapor.clone(), &[("source-vapor-method-wet", state)]),
            value,
        );
    }
    refuses(
        pack,
        "water-vapor-diffusion-resistance-factor",
        Dims::NONE,
        &vapor,
    );
    refuses(
        pack,
        "water-vapor-diffusion-resistance-factor",
        Dims::NONE,
        &flags(vapor, &[("source-vapor-method-wet", 0.5)]),
    );
}

#[test]
fn g0_g3_wisa_spruce_18mm_keeps_layup_and_axial_bending_distinctions() {
    let pack = compile(WISA, 15);
    let p = base(WISA);
    let density = sample(&pack, "density", RHO, &p, 460.0);
    sample(&pack, "characteristic-density-for-racking", RHO, &p, 420.0);
    // Areal mass remains a calculation from nominal geometry, not a new measured claim.
    assert!((density * 0.018 - 8.28).abs() < 1e-12);
    for (parallel, bending, axial, fb, fc, ft) in [
        (1.0, 8170e6, 6682e6, 20.4e6, 16.7e6, 10e6),
        (0.0, 3830e6, 5318e6, 13e6, 13.3e6, 8e6),
    ] {
        let q = loading(p.clone(), parallel, None);
        for (name, value) in [
            ("bending-modulus", bending),
            ("axial-modulus", axial),
            ("characteristic-bending-strength", fb),
            ("characteristic-compression-strength", fc),
            ("characteristic-tensile-strength", ft),
        ] {
            sample(&pack, name, PA, &q, value);
        }
    }
    refuses(&pack, "bending-modulus", PA, &p);
    refuses(&pack, "bending-modulus", PA, &loading(p.clone(), 0.5, None));
    refuses(&pack, "young-modulus", PA, &loading(p.clone(), 1.0, None));
    refuses(
        &pack,
        "density",
        RHO,
        &flags(p.clone(), &[("surface-is-coated", 1.0)]),
    );
    check_transport(&pack, &p, "thermal-conductivity", 66.0, 190.0);
    check_identity_and_unstated_conditions(&pack, &p);
}

fn check_lvl(pack: &NormalizedPack, slug: &str, is_x: bool) {
    let p = base(slug);
    sample(pack, "density", RHO, &p, if is_x { 530.0 } else { 550.0 });
    sample(
        pack,
        "characteristic-density-fifth-percentile",
        RHO,
        &p,
        480.0,
    );
    for edge in [1.0, 0.0] {
        let q = loading(p.clone(), 1.0, Some(edge));
        sample(
            pack,
            "bending-modulus",
            PA,
            &q,
            if is_x { 10600e6 } else { 14000e6 },
        );
        sample(
            pack,
            "characteristic-bending-modulus",
            PA,
            &q,
            if is_x { 9000e6 } else { 12000e6 },
        );
        let strength = if edge == 1.0 {
            refuses(pack, "characteristic-bending-strength", PA, &q);
            refuses(
                pack,
                "characteristic-bending-strength",
                PA,
                &length(q.clone(), "reference-member-depth", 0.6),
            );
            length(q.clone(), "reference-member-depth", 0.3)
        } else {
            q.clone()
        };
        let fb = if edge == 1.0 {
            if is_x { 32e6 } else { 44e6 }
        } else if is_x {
            36e6
        } else {
            50e6
        };
        sample(pack, "characteristic-bending-strength", PA, &strength, fb);
        let g = if edge == 1.0 {
            600e6
        } else if is_x {
            150e6
        } else {
            560e6
        };
        sample(pack, "shear-modulus", PA, &q, g);
        sample(
            pack,
            "characteristic-shear-modulus",
            PA,
            &q,
            if edge == 0.0 && is_x { 130e6 } else { 400e6 },
        );
        sample(
            pack,
            "characteristic-shear-strength",
            PA,
            &q,
            if edge == 1.0 {
                4.6e6
            } else if is_x {
                1.1e6
            } else {
                2.6e6
            },
        );
        sample(
            pack,
            "characteristic-compression-strength",
            PA,
            &loading(p.clone(), 0.0, Some(edge)),
            if edge == 1.0 {
                if is_x { 9e6 } else { 7.5e6 }
            } else if is_x {
                4e6
            } else {
                3.6e6
            },
        );
    }
    let axial = loading(p.clone(), 1.0, None);
    sample(
        pack,
        "characteristic-compression-strength",
        PA,
        &axial,
        if is_x { 30e6 } else { 40e6 },
    );
    refuses(pack, "characteristic-tensile-strength", PA, &axial);
    sample(
        pack,
        "characteristic-tensile-strength",
        PA,
        &length(axial.clone(), "reference-member-length", 3.0),
        if is_x { 18e6 } else { 36e6 },
    );
    refuses(
        pack,
        "characteristic-tensile-strength",
        PA,
        &length(axial, "reference-member-length", 6.0),
    );
    sample(
        pack,
        "characteristic-tensile-strength",
        PA,
        &loading(p.clone(), 0.0, Some(1.0)),
        if is_x { 5e6 } else { 0.9e6 },
    );
    refuses(pack, "bending-modulus", PA, &loading(p.clone(), 1.0, None));
    refuses(
        pack,
        "bending-modulus",
        PA,
        &loading(p.clone(), 1.0, Some(0.5)),
    );
    refuses(
        pack,
        "density",
        RHO,
        &flags(
            p.clone(),
            &[("cross-ply-count", if is_x { 0.0 } else { 3.0 })],
        ),
    );
    check_transport(pack, &p, "design-thermal-conductivity", 75.0, 205.0);
    refuses(
        pack,
        "thermal-conductivity",
        K,
        &flags(p.clone(), &[("source-conductivity-direction-known", 0.0)]),
    );
    check_identity_and_unstated_conditions(pack, &p);

    let moisture = flags(
        p,
        &[
            ("source-below-fiber-saturation", 1.0),
            ("source-fiber-saturation-value-known", 0.0),
        ],
    );
    let name = "swelling-shrinkage-strain-per-wood-moisture-fraction-change";
    for (direction, coefficient) in [
        (1.0, 0.01),
        (2.0, if is_x { 0.03 } else { 0.32 }),
        (3.0, 0.32),
    ] {
        let q = flags(
            moisture.clone(),
            &[("dimensional-change-direction", direction)],
        );
        sample(pack, name, Dims::NONE, &q, coefficient);
        refuses(
            pack,
            name,
            Dims::NONE,
            &flags(q, &[("source-below-fiber-saturation", 0.0)]),
        );
    }
    refuses(pack, name, Dims::NONE, &moisture);
    refuses(
        pack,
        name,
        Dims::NONE,
        &flags(moisture, &[("dimensional-change-direction", 1.5)]),
    );
}

#[test]
fn g0_g3_lvl_r_keeps_transverse_compression_modulus_out_of_bending() {
    let pack = compile(LVL_R, 27);
    check_lvl(&pack, LVL_R, false);
    let p = loading(base(LVL_R), 0.0, Some(1.0));
    sample(&pack, "compression-modulus", PA, &p, 430e6);
    sample(&pack, "characteristic-compression-modulus", PA, &p, 350e6);
    refuses(&pack, "bending-modulus", PA, &p);
    refuses(
        &pack,
        "characteristic-bending-strength",
        PA,
        &length(p.clone(), "reference-member-depth", 0.3),
    );
    refuses(&pack, "shear-modulus", PA, &p);
    refuses(&pack, "characteristic-shear-strength", PA, &p);
}

#[test]
fn g0_g3_lvl_x_uses_thick_cross_laminated_column_and_discrete_load_cases() {
    let pack = compile(LVL_X, 35);
    check_lvl(&pack, LVL_X, true);
    for (edge, mean, characteristic) in [(1.0, 3000e6, 2300e6), (0.0, 2500e6, 1800e6)] {
        let p = loading(base(LVL_X), 0.0, Some(edge));
        sample(&pack, "bending-modulus", PA, &p, mean);
        sample(
            &pack,
            "characteristic-bending-modulus",
            PA,
            &p,
            characteristic,
        );
        let strength = if edge == 1.0 {
            length(p.clone(), "reference-member-depth", 0.3)
        } else {
            p.clone()
        };
        sample(&pack, "characteristic-bending-strength", PA, &strength, 8e6);
        sample(
            &pack,
            "characteristic-shear-strength",
            PA,
            &p,
            if edge == 1.0 { 4.6e6 } else { 1.1e6 },
        );
        if edge == 0.0 {
            sample(&pack, "shear-modulus", PA, &p, 150e6);
            sample(&pack, "characteristic-shear-modulus", PA, &p, 130e6);
        } else {
            refuses(&pack, "shear-modulus", PA, &p);
        }
        refuses(&pack, "compression-modulus", PA, &p);
    }
}
