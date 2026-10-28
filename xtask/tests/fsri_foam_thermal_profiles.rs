//! Primary HFM foam curves retain their measurement basis and explicit local-model opt-in.
#![deny(unsafe_code)]

use fs_matdb::{
    MaterialAnswer, MaterialCard, MaterialStateId, NormalizedMaterialCardPack, NormalizedPack,
    PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use std::{path::Path, process::Command};

const C: Dims = Dims([-1, 1, -2, -1, 0, 0]);
const CP: Dims = Dims([2, 0, -2, -1, 0, 0]);
const RHO: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const K: Dims = Dims([1, 1, -3, -1, 0, 0]);

struct Foam {
    id: &'static str,
    stock: f64,
    densities: [f64; 3],
    capacity: [[f64; 4]; 3],
    cp: [f64; 4],
    cp_sd: [f64; 4],
    k: [f64; 2],
}

const FOAMS: [Foam; 4] = [
    Foam {
        id: "pe-foam-pipe",
        stock: 0.0135,
        densities: [35.0, 35.0, 36.0],
        capacity: [
            [76013.0, 79483.0, 84135.0, 89786.0],
            [76310.0, 80500.0, 84687.0, 91536.0],
            [76709.0, 80564.0, 86028.0, 92244.0],
        ],
        cp: [2161.0, 2269.0, 2404.0, 2581.0],
        cp_sd: [10.0, 17.0, 28.0, 36.0],
        k: [0.049, 0.058],
    },
    Foam {
        id: "rubber-foam-pipe",
        stock: 0.012,
        densities: [52.0, 52.0, 51.0],
        capacity: [
            [68116.0, 69732.0, 72042.0, 74765.0],
            [69293.0, 70005.0, 73986.0, 76291.0],
            [69043.0, 70991.0, 73602.0, 76075.0],
        ],
        cp: [1332.0, 1360.0, 1417.0, 1465.0],
        cp_sd: [12.0, 13.0, 20.0, 16.0],
        k: [0.034, 0.038],
    },
    Foam {
        id: "memory-foam-pad",
        stock: 0.0128,
        densities: [136.0, 137.0, 130.0],
        capacity: [
            [240014.0, 247027.0, 251688.0, 254275.0],
            [240566.0, 247201.0, 252307.0, 254538.0],
            [231341.0, 239049.0, 241341.0, 245492.0],
        ],
        cp: [1767.0, 1820.0, 1849.0, 1872.0],
        cp_sd: [39.0, 35.0, 46.0, 38.0],
        k: [0.039, 0.042],
    },
    Foam {
        id: "xps-board",
        stock: 0.0254,
        densities: [27.0, 28.0, 28.0],
        capacity: [
            [34276.0, 35606.0, 36667.0, 37297.0],
            [36018.0, 37169.0, 38656.0, 39199.0],
            [35710.0, 36349.0, 38272.0, 38080.0],
        ],
        cp: [1277.0, 1315.0, 1369.0, 1380.0],
        cp_sd: [34.0, 28.0, 38.0, 35.0],
        k: [0.028, 0.032],
    },
];

fn compile(f: &Foam) -> NormalizedMaterialCardPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let slug = format!("fsri-{}-unconditioned-2022", f.id);
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let path = std::env::temp_dir().join(format!(
        "fsri-foam-{}-{nonce}-{}.fsmatpk",
        std::process::id(),
        f.id
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
    assert_eq!(pack.claims().claims_ordered().count(), 8);
    let card = NormalizedMaterialCardPack::new(
        MaterialStateId {
            chemistry: slug,
            phase: "porous solid".into(),
            process: "FSRI unconditioned January 2022 specimen population".into(),
            revision: 0,
        },
        pack,
    )
    .unwrap();
    let bytes = card.to_bytes();
    let restored =
        NormalizedMaterialCardPack::from_bytes_verified(card.content_hash(), &bytes).unwrap();
    assert_eq!(restored.to_bytes(), bytes);
    restored
}

fn semantic(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

fn axis(p: QueryPoint, name: &str, value: f64) -> QueryPoint {
    p.with_quantity(name, QuantitySpec::dimensional(Dims::NONE), value)
        .unwrap()
}

fn temperature(p: QueryPoint, name: &str, kelvin: f64) -> QueryPoint {
    p.with_quantity(name, semantic(QuantityKind::AbsoluteTemperature), kelvin)
        .unwrap()
}

fn point(f: &Foam, kelvin: f64) -> QueryPoint {
    let mut p = QueryPoint::new();
    for (name, value) in [
        (format!("source-fsri-{}", f.id), 1.0),
        ("source-unconditioned".into(), 1.0),
        ("source-grade-known".into(), 0.0),
        ("source-formulation-known".into(), 0.0),
        ("source-moisture-content-known".into(), 0.0),
        ("source-contact-pressure-known".into(), 0.0),
        ("source-statistic-is-mean".into(), 1.0),
        ("source-specimen-count".into(), 3.0),
        ("source-hfm-through-thickness".into(), 1.0),
        ("source-hfm-upward-heat-flow".into(), 1.0),
        ("source-density-temperature-known".into(), 0.0),
        ("mean-test-temperature-as-local-model".into(), 1.0),
    ] {
        p = axis(p, &name, value);
    }
    p = p
        .with_quantity(
            "source-stock-thickness",
            QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0])),
            f.stock,
        )
        .unwrap();
    p = p
        .with_quantity(
            "source-hfm-temperature-difference",
            semantic(QuantityKind::TemperatureDifference),
            20.0,
        )
        .unwrap();
    temperature(p, "temperature", kelvin)
}

fn key(name: &str, dims: Dims) -> PropertyKey {
    let quantity = if matches!(
        name,
        "thermal-conductivity" | "apparent-thermal-conductivity"
    ) {
        semantic(QuantityKind::ThermalConductivity)
    } else {
        QuantitySpec::dimensional(dims)
    };
    PropertyKey::with_quantity(name, quantity)
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-12 * expected.abs().max(1e-10),
        "{actual} != {expected}"
    );
}

fn sample(
    card: &MaterialCard,
    name: &str,
    dims: Dims,
    p: &QueryPoint,
    expected: f64,
) -> MaterialAnswer {
    let answer = card
        .claims()
        .query_typed(&key(name, dims), p, SelectionPolicy::SingleClaimOnly)
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
    answer
}

fn refuses(card: &MaterialCard, name: &str, dims: Dims, p: &QueryPoint) {
    assert!(
        card.claims()
            .query_typed(&key(name, dims), p, SelectionPolicy::SingleClaimOnly)
            .is_err(),
        "must refuse {name} at {p:?}"
    );
}

#[test]
fn g0_g3_foam_native_volume_capacity_and_published_cp_keep_their_normalization() {
    for f in &FOAMS {
        let pack = compile(f);
        let card = pack.card();
        let density = f.densities.iter().sum::<f64>() / 3.0;
        sample(
            card,
            "reference-normalization-density",
            RHO,
            &point(f, 303.15),
            density,
        );
        for j in 0..4 {
            let p = point(f, 283.15 + 10.0 * j as f64);
            let mean = f.capacity.iter().map(|row| row[j]).sum::<f64>() / 3.0;
            sample(card, "volumetric-heat-capacity", C, &p, mean);
            sample(card, "specific-heat-capacity", CP, &p, f.cp[j]);
            assert_eq!((mean / density).round(), f.cp[j]);
            assert_ne!(mean, density * f.cp[j], "published rounding must survive");
        }
        let p = point(f, 298.15);
        let expected_c = f
            .capacity
            .iter()
            .map(|row| (row[1] + row[2]) / 2.0)
            .sum::<f64>()
            / 3.0;
        sample(card, "volumetric-heat-capacity", C, &p, expected_c);
        sample(
            card,
            "specific-heat-capacity",
            CP,
            &p,
            (f.cp[1] + f.cp[2]) / 2.0,
        );
        sample(
            card,
            "thermal-conductivity",
            K,
            &p,
            f.k[0] + (f.k[1] - f.k[0]) / 3.0,
        );
        refuses(card, "density", RHO, &p);
        refuses(card, "thermal-diffusivity", Dims([2, 0, -1, 0, 0, 0]), &p);
    }
}

#[test]
fn g0_g3_reported_scatter_remains_a_discrete_statistic_not_measurement_confidence() {
    for f in &FOAMS {
        let pack = compile(f);
        let card = pack.card();
        let density = f.densities.iter().sum::<f64>() / 3.0;
        for j in 0..4 {
            let p = axis(
                point(f, 283.15 + 10.0 * j as f64),
                "source-statistic-is-mean",
                0.0,
            );
            let mean = f.capacity.iter().map(|row| row[j]).sum::<f64>() / 3.0;
            let sd = (f
                .capacity
                .iter()
                .map(|row| (row[j] - mean).powi(2))
                .sum::<f64>()
                / 2.0)
                .sqrt();
            sample(
                card,
                "computed-volumetric-heat-capacity-sample-standard-deviation",
                C,
                &p,
                sd,
            );
            sample(
                card,
                "reported-specific-heat-capacity-standard-deviation",
                CP,
                &p,
                f.cp_sd[j],
            );
            assert_eq!((sd / density).round(), f.cp_sd[j]);
        }
        for t in [288.15, 318.15] {
            let p = axis(
                temperature(point(f, t), "test-mean-temperature", t),
                "source-statistic-is-mean",
                0.0,
            );
            sample(
                card,
                "reported-apparent-thermal-conductivity-standard-deviation",
                K,
                &p,
                0.0,
            );
        }
        let p = axis(
            temperature(point(f, 298.15), "test-mean-temperature", 298.15),
            "source-statistic-is-mean",
            0.0,
        );
        for (name, dims) in [
            (
                "computed-volumetric-heat-capacity-sample-standard-deviation",
                C,
            ),
            ("reported-specific-heat-capacity-standard-deviation", CP),
            (
                "reported-apparent-thermal-conductivity-standard-deviation",
                K,
            ),
        ] {
            refuses(card, name, dims, &p);
        }
    }
}

#[test]
fn g0_g3_foam_state_direction_temperature_and_local_model_boundaries_refuse() {
    for f in &FOAMS {
        let pack = compile(f);
        let card = pack.card();
        let p = point(f, 303.15);
        refuses(card, "apparent-thermal-conductivity", K, &p);
        let native = temperature(p.clone(), "test-mean-temperature", 303.15);
        let a = sample(
            card,
            "apparent-thermal-conductivity",
            K,
            &native,
            (f.k[0] + f.k[1]) / 2.0,
        );
        let b = sample(card, "thermal-conductivity", K, &p, (f.k[0] + f.k[1]) / 2.0);
        assert_ne!(a.receipt.selected, b.receipt.selected);
        for (name, value) in [
            ("source-unconditioned", 0.0),
            ("source-moisture-content-known", 1.0),
            ("source-grade-known", 1.0),
            ("source-formulation-known", 1.0),
            ("source-contact-pressure-known", 1.0),
            ("source-specimen-count", 9.0),
        ] {
            let changed = axis(p.clone(), name, value);
            refuses(card, "volumetric-heat-capacity", C, &changed);
            refuses(card, "thermal-conductivity", K, &changed);
        }
        for name in [
            "source-hfm-through-thickness",
            "source-hfm-upward-heat-flow",
            "mean-test-temperature-as-local-model",
        ] {
            refuses(card, "thermal-conductivity", K, &axis(p.clone(), name, 0.0));
        }
        for t in [283.14, 313.16] {
            refuses(card, "volumetric-heat-capacity", C, &point(f, t));
            refuses(card, "specific-heat-capacity", CP, &point(f, t));
        }
        for t in [288.14, 318.16] {
            refuses(card, "thermal-conductivity", K, &point(f, t));
        }
        let changed = p
            .clone()
            .with_quantity(
                "source-hfm-temperature-difference",
                semantic(QuantityKind::TemperatureDifference),
                25.0,
            )
            .unwrap();
        refuses(card, "thermal-conductivity", K, &changed);
        let changed = p
            .with_quantity(
                "source-stock-thickness",
                QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0])),
                f.stock * 2.0,
            )
            .unwrap();
        refuses(card, "volumetric-heat-capacity", C, &changed);
    }
}
