//! Original reported-carbon-mass adsorption data, without an assumed dry basis.
//! A native material compiler/query test, not a kinetic or plant-calibration test.
#![deny(unsafe_code)]

use std::{
    path::{Path, PathBuf},
    process::Command,
    sync::OnceLock,
};

use fs_matdb::{
    MaterialAnswer, NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

const PACK: &str = "calgon-dg11-53um-gold-cyanide-martinez2023-reference";
const MASS_CONCENTRATION: Dims = Dims([-3, 1, 0, 0, 0, 0]);
const AFFINITY: Dims = Dims([3, -1, 0, 0, 0, 0]);
const MASS: Dims = Dims([0, 1, 0, 0, 0, 0]);
const LENGTH: Dims = Dims([1, 0, 0, 0, 0, 0]);

fn workspace() -> PathBuf {
    std::env::var_os("CARGO_WORKSPACE_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| {
            Path::new(env!("CARGO_MANIFEST_DIR"))
                .parent()
                .unwrap()
                .into()
        })
}

fn pack() -> &'static NormalizedPack {
    static PACKED: OnceLock<NormalizedPack> = OnceLock::new();
    PACKED.get_or_init(|| {
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let out =
            std::env::temp_dir().join(format!("frankensim-calgon-{}-{nonce}", std::process::id()));
        std::fs::create_dir(&out).unwrap();
        let target = out.join(format!("{PACK}.fsmatpk"));
        let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
            .args(["matdb-pack", "--manifest"])
            .arg(workspace().join(format!("data/matdb/seed-v1/{PACK}/manifest.tsv")))
            .arg("--out")
            .arg(&target)
            .env("CARGO_WORKSPACE_DIR", workspace())
            .current_dir(workspace())
            .output()
            .unwrap();
        assert!(
            run.status.success(),
            "{}",
            String::from_utf8_lossy(&run.stderr)
        );
        let bytes = std::fs::read(target).unwrap();
        let decoded = NormalizedPack::from_bytes(&bytes).unwrap();
        let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap();
        assert_eq!(pack.pack_id(), PACK);
        assert_eq!(pack.to_bytes(), bytes);
        assert_eq!(pack.claims().claim_count(), 12);
        pack
    })
}

fn point(fit: bool, carbon_mass_kg: Option<f64>) -> QueryPoint {
    let mut at = QueryPoint::new()
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            298.15,
        )
        .unwrap();
    for (name, dims, value) in [
        ("source-reported-pH-lower", Dims::NONE, 10.5),
        ("source-reported-pH-upper", Dims::NONE, 11.0),
        ("source-exact-pH-known", Dims::NONE, 0.0),
        ("source-carbon-sieve-size", LENGTH, 53e-6),
        (
            "source-nacn-preparation-mass-concentration",
            MASS_CONCENTRATION,
            1.0,
        ),
        ("source-calgon-dg11-hb184an4", Dims::NONE, 1.0),
        ("source-dry-carbon-mass-basis-established", Dims::NONE, 0.0),
        ("source-synthetic-gold-cyanide-medium", Dims::NONE, 1.0),
        ("engineering-reference-model", Dims::NONE, 1.0),
        (
            "source-published-langmuir-fit",
            Dims::NONE,
            if fit { 1.0 } else { 0.0 },
        ),
    ] {
        at = at
            .with_quantity(name, QuantitySpec::dimensional(dims), value)
            .unwrap();
    }
    if let Some(mass) = carbon_mass_kg {
        at = at
            .with_quantity("source-carbon-mass", QuantitySpec::dimensional(MASS), mass)
            .unwrap();
    }
    at
}

fn answer(name: &str, dims: Dims, at: &QueryPoint) -> MaterialAnswer {
    let answer = pack()
        .claims()
        .query_typed(
            &PropertyKey::new(name, dims),
            at,
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
    assert_eq!(
        answer.evidence.value.uncertainty,
        UncertaintyModel::Unstated
    );
    assert!(answer.receipt.observation_backed && !answer.receipt.source_hashes.is_empty());
    pack().claims().verify_receipt(&answer.receipt).unwrap();
    answer
}

fn value(name: &str, dims: Dims, at: &QueryPoint) -> f64 {
    answer(name, dims, at).evidence.value.value
}

fn relative(actual: f64, expected: f64) {
    assert!(actual.is_finite() && actual > 0.0 && expected > 0.0);
    assert!(
        (actual / expected - 1.0).abs() < 2e-13,
        "{actual} != {expected}"
    );
}

#[test]
fn g0_published_fit_keeps_reported_mass_units_and_receipts() {
    let at = point(true, None);
    relative(
        value("langmuir-saturation-mass-loading", Dims::NONE, &at),
        0.039_525_7,
    );
    relative(value("langmuir-affinity-mass-basis", AFFINITY, &at), 3953.1);
    relative(
        value(
            "langmuir-half-saturation-mass-concentration",
            MASS_CONCENTRATION,
            &at,
        ),
        0.000_252_966_026_662_619_2,
    );
    relative(
        value("langmuir-linearized-fit-r-squared", Dims::NONE, &at),
        0.9737,
    );
}

#[test]
fn g3_four_observed_pairs_close_original_gold_mass_accounts_without_interpolation() {
    for (carbon_g, ce_mg_l, qe_mg_g) in [
        (0.010, 9.258, 37.100),
        (0.050, 5.766, 42.340),
        (0.125, 1.837, 32.652),
        (0.250, 0.244, 19.512),
    ] {
        let carbon_kg = carbon_g * 0.001;
        let at = point(false, Some(carbon_kg));
        let concentration = value(
            "observed-equilibrium-gold-mass-concentration",
            MASS_CONCENTRATION,
            &at,
        );
        let loading = value("observed-gold-mass-loading", Dims::NONE, &at);
        relative(concentration, ce_mg_l * 0.001);
        relative(loading, qe_mg_g * 0.001);
        // Table 3 uses the nominal 0.5 L bath and 10 mg/L initial Au.
        // This checks the printed summary, not a reconstructed sampled-volume history.
        let final_gold_kg = concentration * 0.0005 + loading * carbon_kg;
        relative(final_gold_kg, 5e-6);
    }
    for at in [point(false, None), point(false, Some(0.0001))] {
        assert!(
            pack()
                .claims()
                .query_typed(
                    &PropertyKey::new("observed-gold-mass-loading", Dims::NONE),
                    &at,
                    SelectionPolicy::SingleClaimOnly,
                )
                .is_err()
        );
    }
}

#[test]
fn g3_half_saturation_is_algebraic_and_fitted_capacity_is_not_a_measured_hard_cap() {
    let at = point(true, None);
    let qmax = value("langmuir-saturation-mass-loading", Dims::NONE, &at);
    let affinity = value("langmuir-affinity-mass-basis", AFFINITY, &at);
    let half = value(
        "langmuir-half-saturation-mass-concentration",
        MASS_CONCENTRATION,
        &at,
    );
    relative(affinity * half, 1.0);
    relative(qmax * affinity * half / (1.0 + affinity * half), qmax / 2.0);
    let observed = point(false, Some(0.00005));
    let loading = value("observed-gold-mass-loading", Dims::NONE, &observed);
    let concentration = value(
        "observed-equilibrium-gold-mass-concentration",
        MASS_CONCENTRATION,
        &observed,
    );
    assert!(loading > qmax);
    let predicted = qmax * affinity * concentration / (1.0 + affinity * concentration);
    assert!(predicted < qmax && (predicted - loading).abs() > 0.001);
    // The fit's R² is not a claim that its upper asymptote bounds every datum.
}

#[test]
fn g0_wrong_chemistry_particle_temperature_and_dry_basis_refuse() {
    let original = point(true, None);
    let bad = [
        QueryPoint::new(),
        original
            .clone()
            .with("source-dry-carbon-mass-basis-established", 1.0)
            .unwrap(),
        original
            .clone()
            .with("source-calgon-dg11-hb184an4", 0.0)
            .unwrap(),
        original
            .clone()
            .with("source-synthetic-gold-cyanide-medium", 0.0)
            .unwrap(),
        original
            .clone()
            .with("source-published-langmuir-fit", 0.0)
            .unwrap(),
        original
            .clone()
            .with("source-reported-pH-lower", 9.0)
            .unwrap(),
        original
            .clone()
            .with_quantity(
                "source-carbon-sieve-size",
                QuantitySpec::dimensional(LENGTH),
                75e-6,
            )
            .unwrap(),
        original
            .clone()
            .with_quantity(
                "source-nacn-preparation-mass-concentration",
                QuantitySpec::dimensional(MASS_CONCENTRATION),
                0.5,
            )
            .unwrap(),
        original
            .clone()
            .with_quantity(
                "temperature",
                QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::AbsoluteTemperature,
                    ValueForm::Static,
                )),
                310.0,
            )
            .unwrap(),
    ];
    for at in bad {
        assert!(
            pack()
                .claims()
                .query_typed(
                    &PropertyKey::new("langmuir-saturation-mass-loading", Dims::NONE),
                    &at,
                    SelectionPolicy::SingleClaimOnly,
                )
                .is_err()
        );
    }
    for (key, dims) in [
        ("saturation-loading", Dims([0, -1, 0, 0, 0, 1])),
        ("half-saturation-concentration", Dims([-3, 0, 0, 0, 0, 1])),
        (
            "linear-driving-force-coefficient",
            Dims([0, 0, -1, 0, 0, 0]),
        ),
    ] {
        assert!(
            pack()
                .claims()
                .query_typed(
                    &PropertyKey::new(key, dims),
                    &original,
                    SelectionPolicy::SingleClaimOnly,
                )
                .is_err(),
            "invented dry-molar or kinetic property {key}"
        );
    }
}
