//! G0/G3: actual compiler -> immutable varnish packs -> typed source queries.
//! These checks establish data transport, not historical acoustic validation.
#![deny(unsafe_code)]

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{FrequencyConvention, QuantityKind, SemanticType, ValueForm},
};
use std::{path::Path, process::Command};

const PA: Dims = Dims([-1, 1, -2, 0, 0, 0]);
const RHO: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const SECOND: Dims = Dims([0, 0, 1, 0, 0, 0]);
const METRE: Dims = Dims([1, 0, 0, 0, 0, 0]);

fn compile(slug: &str, count: usize) -> NormalizedPack {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "frankensim-varnish-{}-{nonce}-{slug}",
        std::process::id()
    ));
    std::fs::create_dir(&dir).unwrap();
    let out = dir.join("material.fsmatpk");
    let result = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(root.join(format!("data/matdb/seed-v1/{slug}/manifest.tsv")))
        .arg("--out")
        .arg(&out)
        .env("CARGO_WORKSPACE_DIR", root)
        .current_dir(root)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{slug}: {}",
        String::from_utf8_lossy(&result.stderr)
    );
    let bytes = std::fs::read(out).unwrap();
    let decoded = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap();
    assert_eq!(pack.to_bytes(), bytes);
    assert_eq!(pack.pack_id(), slug);
    assert_eq!(pack.claims().claims_ordered().count(), count);
    pack
}

fn flag(point: QueryPoint, name: &str, value: f64) -> QueryPoint {
    point
        .with_quantity(name, QuantitySpec::dimensional(Dims::NONE), value)
        .unwrap()
}

fn semantic(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

fn point(kelvin: f64, rh: f64, hz: f64, matched: f64) -> QueryPoint {
    let p = QueryPoint::new()
        .with_quantity(
            "temperature",
            semantic(QuantityKind::AbsoluteTemperature),
            kelvin,
        )
        .unwrap()
        .with_quantity(
            "frequency",
            semantic(QuantityKind::Frequency(FrequencyConvention::Cyclic)),
            hz,
        )
        .unwrap();
    flag(
        flag(
            flag(p, "conditioning_relative_humidity", rh),
            "source_fitted_response",
            1.,
        ),
        "source_density_sample_matched",
        matched,
    )
}

fn value(pack: &NormalizedPack, name: &str, dims: Dims, p: &QueryPoint) -> f64 {
    let answer = pack
        .claims()
        .query_typed(
            &PropertyKey::with_quantity(name, QuantitySpec::dimensional(dims)),
            p,
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap_or_else(|e| panic!("{} {name}: {e}", pack.pack_id()));
    assert_eq!(answer.evidence.value.dims, dims);
    assert!(matches!(
        answer.evidence.value.uncertainty,
        UncertaintyModel::Unstated
    ));
    assert!(answer.receipt.observation_backed);
    assert!(!answer.receipt.source_hashes.is_empty());
    pack.claims().verify_receipt(&answer.receipt).unwrap();
    answer.evidence.value.value
}

fn refuses(pack: &NormalizedPack, name: &str, dims: Dims, p: &QueryPoint) {
    assert!(
        pack.claims()
            .query_typed(
                &PropertyKey::with_quantity(name, QuantitySpec::dimensional(dims)),
                p,
                SelectionPolicy::SingleClaimOnly,
            )
            .is_err(),
        "{} unexpectedly admits {name}",
        pack.pack_id()
    );
}

#[test]
fn g0_six_varnish_samples_preserve_published_coefficients_and_density_statistics() {
    for (family, sample, e0, einf, alpha, beta, tau, rho, sd, n) in [
        (
            "oil-copal",
            1,
            20e6,
            1.37e9,
            0.391,
            0.360,
            0.127,
            1097.,
            26.,
            5.,
        ),
        (
            "oil-copal",
            2,
            20e6,
            1.41e9,
            0.402,
            0.370,
            0.075,
            1097.,
            26.,
            5.,
        ),
        (
            "spirit-shellac",
            1,
            70e6,
            1.89e9,
            0.330,
            0.314,
            0.125,
            984.,
            30.,
            7.,
        ),
        (
            "spirit-shellac",
            2,
            70e6,
            1.83e9,
            0.335,
            0.313,
            0.179,
            984.,
            30.,
            7.,
        ),
        (
            "mineral-ground",
            1,
            500e6,
            2.08e9,
            0.309,
            0.288,
            0.268,
            1155.,
            30.,
            4.,
        ),
        (
            "mineral-ground",
            2,
            500e6,
            2.12e9,
            0.346,
            0.318,
            0.132,
            1155.,
            30.,
            4.,
        ),
    ] {
        let pack = compile(&format!("varnish-{family}-lammlein2021-s{sample}"), 8);
        for hz in [20., 1000., 20000.] {
            let p = point(293.15, 0.35, hz, 0.);
            for (name, dims, expected) in [
                ("varnish_e0", PA, e0),
                ("varnish_einf", PA, einf),
                ("varnish_alpha", Dims::NONE, alpha),
                ("varnish_beta", Dims::NONE, beta),
                ("varnish_tau", SECOND, tau),
                ("density", RHO, rho),
                ("density_standard_deviation", RHO, sd),
                ("density_sample_count", Dims::NONE, n),
            ] {
                assert_eq!(
                    value(&pack, name, dims, &p),
                    expected,
                    "{family} sample{sample} {name}"
                );
            }
        }
    }
}

#[test]
fn g3_varnish_refuses_extrapolation_convention_changes_and_invented_specimen_matching() {
    let pack = compile("varnish-oil-copal-lammlein2021-s1", 8);
    for p in [
        point(294.15, 0.35, 1000., 0.),
        point(293.15, 0.65, 1000., 0.),
        point(293.15, 0.35, 19.99, 0.),
        point(293.15, 0.35, 20000.01, 0.),
        QueryPoint::new(),
    ] {
        refuses(&pack, "varnish_einf", PA, &p);
    }
    refuses(&pack, "density", RHO, &point(293.15, 0.35, 1000., 1.));
    let angular = point(293.15, 0.35, 1000., 0.)
        .with_quantity(
            "frequency",
            semantic(QuantityKind::Frequency(FrequencyConvention::Angular)),
            1000.,
        )
        .unwrap();
    refuses(&pack, "varnish_einf", PA, &angular);
    refuses(
        &pack,
        "poisson_ratio",
        Dims::NONE,
        &point(293.15, 0.35, 1000., 0.),
    );
}

#[test]
fn g3_historical_layer_ranges_are_categorical_and_do_not_supply_mechanical_properties() {
    let pack = compile("varnish-stradivari-echard2010-stratigraphy", 4);
    let p = |code| {
        flag(
            flag(QueryPoint::new(), "specimen_code", code),
            "source_stratigraphy_only",
            1.,
        )
    };
    for code in 1..=5 {
        assert_eq!(
            value(
                &pack,
                "external_coating_thickness_observed_lower",
                METRE,
                &p(f64::from(code))
            ),
            10e-6
        );
        assert_eq!(
            value(
                &pack,
                "external_coating_thickness_observed_upper",
                METRE,
                &p(f64::from(code))
            ),
            30e-6
        );
    }
    for (code, lower, upper) in [
        (1., 30e-6, 100e-6),
        (2., 10e-6, 30e-6),
        (3., 10e-6, 30e-6),
        (4., 30e-6, 100e-6),
    ] {
        assert_eq!(
            value(
                &pack,
                "penetrating_oil_depth_observed_lower",
                METRE,
                &p(code)
            ),
            lower
        );
        assert_eq!(
            value(
                &pack,
                "penetrating_oil_depth_observed_upper",
                METRE,
                &p(code)
            ),
            upper
        );
    }
    refuses(
        &pack,
        "external_coating_thickness_observed_lower",
        METRE,
        &p(2.5),
    );
    refuses(&pack, "penetrating_oil_depth_observed_upper", METRE, &p(5.));
    refuses(&pack, "density", RHO, &p(2.));
    refuses(&pack, "varnish_einf", PA, &p(2.));
}
