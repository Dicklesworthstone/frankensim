//! Original USGS hydroxide data close water autoprotolysis without mixing NBS.
#![deny(unsafe_code)]

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::OnceLock;

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_qty::{
    Amount, Dimensionless, Dims, MolarMass, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    SpeciesId,
    stoichiometry::{
        minerals::MineralReaction,
        thermodynamics::{
            FormationGibbsSource,
            acid_base::{IdealAcidBaseSystem, SpeciationOptions, WaterAutoprotolysis},
        },
    },
};

const REFERENCE: &str = "source-elemental-reference-usgs-b2131";
const PACKS: [&str; 3] = [
    "water-liquid-usgs-b2131",
    "hydrogen-ion-aqueous-usgs-b2131",
    "hydroxide-ion-aqueous-usgs-b2131",
];
const SPECIES: [&str; 3] = ["H2O(l)", "H+(aq)", "OH[-](aq)"];
const MOLAR_ENERGY: Dims = Dims([2, 1, -2, 0, 0, -1]);
const MOLAR_ENTROPY: Dims = Dims([2, 1, -2, -1, 0, -1]);

fn workspace() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .to_path_buf()
}

fn materials() -> &'static [NormalizedPack; 3] {
    static MATERIALS: OnceLock<[NormalizedPack; 3]> = OnceLock::new();
    MATERIALS.get_or_init(|| {
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let out = std::env::temp_dir().join(format!(
            "frankensim-hydroxide-{}-{nonce}",
            std::process::id()
        ));
        std::fs::create_dir(&out).unwrap();
        PACKS.map(|slug| {
            let target = out.join(format!("{slug}.fsmatpk"));
            let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
                .args(["matdb-pack", "--manifest"])
                .arg(workspace().join(format!("data/matdb/seed-v1/{slug}/manifest.tsv")))
                .arg("--out")
                .arg(&target)
                .env("CARGO_WORKSPACE_DIR", workspace())
                .current_dir(workspace())
                .output()
                .unwrap();
            assert!(
                run.status.success(),
                "{slug}: {}",
                String::from_utf8_lossy(&run.stderr)
            );
            let bytes = std::fs::read(target).unwrap();
            let decoded = NormalizedPack::from_bytes(&bytes).unwrap();
            let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap();
            assert_eq!(pack.pack_id(), slug);
            assert_eq!(pack.to_bytes(), bytes);
            pack
        })
    })
}

fn point(index: usize) -> QueryPoint {
    let mut at = QueryPoint::new()
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            298.15,
        )
        .unwrap()
        .with_quantity(
            "pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            100_000.0,
        )
        .unwrap();
    let mut flags = vec![
        (REFERENCE, 1.0),
        ("source-standard-state", 1.0),
        ("engineering-reference-model", 1.0),
        ("source-specimen-qualified", 0.0),
    ];
    if index == 0 {
        flags.extend([("source-phase-liquid", 1.0), ("source-pure-species", 1.0)]);
    } else {
        flags.extend([
            ("source-phase-aqueous", 1.0),
            ("source-aqueous-unit-activity-infinite-dilution", 1.0),
            ("source-ionic-reference-hydrogen-zero", 1.0),
            ("source-individual-ion", 1.0),
        ]);
    }
    for (axis, value) in flags {
        at = at
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    at
}

#[test]
fn g0_g3_literal_hydroxide_source_retains_conventions_and_reported_magnitudes() {
    let pack = &materials()[2];
    let at = point(2);
    assert_eq!(pack.claims().claim_count(), 6);
    for (name, dims, expected) in [
        ("molar-mass", MolarMass::DIMS, 0.017007),
        (
            "standard-molar-formation-enthalpy",
            MOLAR_ENERGY,
            -230_000.0,
        ),
        (
            "standard-molar-formation-gibbs-energy",
            MOLAR_ENERGY,
            -157_300.0,
        ),
        ("conventional-standard-molar-entropy", MOLAR_ENTROPY, -10.7),
        (
            "reported-formation-enthalpy-uncertainty-magnitude",
            MOLAR_ENERGY,
            100.0,
        ),
        (
            "reported-conventional-entropy-uncertainty-magnitude",
            MOLAR_ENTROPY,
            0.2,
        ),
    ] {
        let key = PropertyKey::new(name, dims);
        let answer = pack
            .claims()
            .query_typed(&key, &at, SelectionPolicy::SingleClaimOnly)
            .unwrap();
        assert!(
            (answer.evidence.value.value - expected).abs() <= 1e-12,
            "{name}"
        );
        assert_eq!(
            answer.evidence.value.uncertainty,
            UncertaintyModel::Unstated
        );
        pack.claims().verify_receipt(&answer.receipt).unwrap();
        assert!(answer.receipt.observation_backed);
        for (axis, value) in [
            ("temperature", 300.0),
            ("pressure", 101_325.0),
            (REFERENCE, 0.0),
            ("source-phase-aqueous", 0.0),
            ("source-ionic-reference-hydrogen-zero", 0.0),
        ] {
            assert!(
                pack.claims()
                    .query_typed(
                        &key,
                        &at.clone().with(axis, value).unwrap(),
                        SelectionPolicy::SingleClaimOnly
                    )
                    .is_err()
            );
        }
    }
}

#[test]
fn g1_usgs_water_autoprotolysis_closes_pure_water_ph_with_original_receipts() {
    let packs = materials();
    let points = std::array::from_fn::<_, 3, _>(point);
    let ids = SPECIES.map(|s| SpeciesId::new(s).unwrap());
    let mut sources = std::array::from_fn::<_, 3, _>(|i| FormationGibbsSource {
        species: &ids[i],
        pack: &packs[i],
        at: &points[i],
        reference_axis: REFERENCE,
    });
    let build = |sources: &[FormationGibbsSource<'_>]| {
        WaterAutoprotolysis::from_sources(
            MineralReaction::WaterDissociation.build().unwrap(),
            ids[0].clone(),
            ids[1].clone(),
            ids[2].clone(),
            sources,
        )
    };
    let water = build(&sources).unwrap();
    assert_eq!(water.standard().molar().value(), 79_800.0);
    assert_eq!(water.standard().reference_axis(), REFERENCE);
    for term in water.standard().terms() {
        packs
            .iter()
            .find(|pack| pack.pack_id() == term.pack_id)
            .unwrap()
            .claims()
            .verify_receipt(&term.answer.receipt)
            .unwrap();
    }
    let system = IdealAcidBaseSystem::new(water, vec![]).unwrap();
    let answer = system
        .speciate(
            &[],
            &[],
            Amount::new(1.0),
            SpeciationOptions {
                ph_tolerance: Dimensionless::new(1e-11),
                absolute_charge_tolerance: Amount::new(1e-14),
                relative_charge_tolerance: Dimensionless::new(1e-12),
                ..SpeciationOptions::default()
            },
            || false,
        )
        .unwrap();
    // Independent 60-digit Decimal calculation using the literal USGS reaction
    // sum and R=8.31446261815324. Rounded USGS data are not forced to pH exactly 7.
    assert!((answer.ph.value() - 6.990_177_571_056_616).abs() < 1e-10);
    for activity in [answer.proton_activity, answer.hydroxide_activity] {
        assert!((activity.value() - 1.022_874_681_473_317_7e-7).abs() < 3e-17);
    }
    assert!(
        (answer.proton_activity.value() * answer.hydroxide_activity.value()
            - 1.046_272_613_999_141_1e-14)
            .abs()
            < 1e-24
    );
    assert!(answer.charge_balance_residual.value().abs() <= 1e-14);
    assert!(answer.families.is_empty());

    sources[2].reference_axis = "source-elemental-reference-nbs1989";
    assert!(
        build(&sources).is_err(),
        "a different reference cannot relabel this source"
    );
}
