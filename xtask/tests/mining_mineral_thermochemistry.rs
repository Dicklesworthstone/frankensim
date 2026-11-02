//! Primary mineral references through the real compiler and typed query seam.
#![deny(unsafe_code)]

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use fs_matdb::{
    ClaimSelection, NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_qty::{
    Dims, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};

const CP: Dims = Dims([2, 0, -2, -1, 0, 0]);
const MOLAR_MASS: Dims = Dims([0, 1, 0, 0, 0, -1]);
const MOLAR_ENERGY: Dims = Dims([2, 1, -2, 0, 0, -1]);
const MOLAR_ENTROPY: Dims = Dims([2, 1, -2, -1, 0, -1]);
const JANAF: &str = "source-elemental-reference-janaf";
const USGS: &str = "source-elemental-reference-usgs-b2131";

struct Mineral {
    pack: &'static str,
    polymorph: &'static str,
    reference: &'static str,
    mass_kg_per_mol: f64,
    formation_j_per_mol: f64,
    entropy_j_per_mol_k: f64,
    // Independently transcribed primary table Cp at 298.15, 300, 400, 500 K.
    cp_molar: Option<[f64; 4]>,
    source_heat_300_400_j_per_mol: f64,
    heating_relative_tolerance: f64,
}

const MINERALS: &[Mineral] = &[
    Mineral {
        pack: "calcite-solid-usgs-b2131-thermochemistry",
        polymorph: "source-phase-calcite",
        reference: USGS,
        mass_kg_per_mol: 0.100087,
        formation_j_per_mol: -1207400.0,
        entropy_j_per_mol_k: 91.71,
        cp_molar: Some([83.47, 83.82, 97.00, 104.55]),
        // USGS B2131 p268 prints (H(T)-H(298.15))/T in J/mol/K:
        // https://digital.library.unt.edu/ark:/67531/metadc1114954/m1/274/
        source_heat_300_400_j_per_mol: 400.0 * 23.18 - 300.0 * 0.52,
        heating_relative_tolerance: 0.009,
    },
    Mineral {
        pack: "pyrite-solid-janaf-thermochemistry",
        polymorph: "source-phase-pyrite",
        reference: JANAF,
        // The selected JANAF printed table's weight, not the WebBook update.
        mass_kg_per_mol: 0.119967,
        formation_j_per_mol: -171544.0,
        entropy_j_per_mol_k: 52.916,
        cp_molar: Some([62.124, 62.289, 68.848, 72.040]),
        // https://janaf.nist.gov/tables/Fe-028.html prints H-H(298.15) in kJ/mol.
        source_heat_300_400_j_per_mol: 1000.0 * (6.705 - 0.115),
        heating_relative_tolerance: 0.006,
    },
    Mineral {
        pack: "quartz-alpha-janaf-standard-reference",
        polymorph: "source-phase-quartz-alpha",
        reference: JANAF,
        mass_kg_per_mol: 0.0600843,
        // https://janaf.nist.gov/tables/O-037.html reference row.
        formation_j_per_mol: -910857.0,
        entropy_j_per_mol_k: 41.463,
        // The existing quartz-alpha-nist-shomate-caloric pack already has Cp.
        cp_molar: None,
        source_heat_300_400_j_per_mol: 0.0,
        heating_relative_tolerance: 0.0,
    },
];

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

fn compile(mineral: &Mineral) -> NormalizedPack {
    static NEXT: AtomicU64 = AtomicU64::new(0);
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let out = std::env::temp_dir().join(format!(
        "mineral-thermochemistry-{}-{nonce}-{}.fsmatpk",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(workspace().join(format!("data/matdb/seed-v1/{}/manifest.tsv", mineral.pack)))
        .arg("--out")
        .arg(&out)
        .current_dir(workspace())
        .env("CARGO_WORKSPACE_DIR", workspace())
        .output()
        .unwrap();
    assert!(
        run.status.success(),
        "{}: {}",
        mineral.pack,
        String::from_utf8_lossy(&run.stderr)
    );
    let bytes = std::fs::read(out).unwrap();
    let raw = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(raw.content_hash(), &bytes).unwrap();
    assert_eq!(pack.to_bytes(), bytes);
    assert_eq!(pack.pack_id(), mineral.pack);
    assert_eq!(
        pack.claims().claims_ordered().count(),
        if mineral.cp_molar.is_some() { 4 } else { 3 }
    );
    pack
}

fn semantic(kind: QuantityKind) -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static))
}

fn point(mineral: &Mineral, temperature: f64) -> QueryPoint {
    let mut point = QueryPoint::new()
        .with_quantity(
            "temperature",
            semantic(QuantityKind::AbsoluteTemperature),
            temperature,
        )
        .unwrap()
        .with_quantity("pressure", semantic(QuantityKind::Pressure), 100000.0)
        .unwrap();
    for (axis, value) in [
        ("source-standard-state", 1.0),
        (mineral.reference, 1.0),
        ("source-phase-solid", 1.0),
        (mineral.polymorph, 1.0),
        ("source-pure-species", 1.0),
        ("engineering-reference-model", 1.0),
        ("source-specimen-qualified", 0.0),
    ] {
        point = point
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    point
}

fn close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

#[test]
fn mineral_reference_values_compile_and_preserve_units_and_sources() {
    for mineral in MINERALS {
        let pack = compile(mineral);
        for (name, dims, expected) in [
            ("molar-mass", MOLAR_MASS, mineral.mass_kg_per_mol),
            (
                "standard-molar-formation-enthalpy",
                MOLAR_ENERGY,
                mineral.formation_j_per_mol,
            ),
            (
                "standard-molar-entropy",
                MOLAR_ENTROPY,
                mineral.entropy_j_per_mol_k,
            ),
        ] {
            let answer = pack
                .claims()
                .query_typed(
                    &PropertyKey::new(name, dims),
                    &point(mineral, 298.15),
                    SelectionPolicy::SingleClaimOnly,
                )
                .unwrap();
            close(answer.evidence.value.value, expected, 2e-14);
            assert_eq!(
                answer.evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
            assert!(answer.receipt.source_hashes.len() >= 3);
            pack.claims().verify_receipt(&answer.receipt).unwrap();
            let selected = pack.claims().claim(answer.receipt.selected).unwrap();
            assert_eq!(selected.validity.bound(mineral.reference), Some((1.0, 1.0)));
        }
    }
}

#[test]
fn mineral_cp_curves_agree_with_printed_points_and_independent_heat_content() {
    let cp = PropertyKey::new("specific-heat-capacity", CP);
    for mineral in MINERALS {
        let Some(source_cp) = mineral.cp_molar else {
            continue;
        };
        let pack = compile(mineral);
        let line = pack
            .claims()
            .query_line_typed(
                &cp,
                "temperature",
                &point(mineral, 298.15),
                &point(mineral, 500.0),
                ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
            )
            .unwrap();
        assert_eq!(line.knots.len(), 4);
        for ((temperature, answer), (source_t, source_cp)) in line
            .knots
            .iter()
            .zip([298.15, 300.0, 400.0, 500.0].into_iter().zip(source_cp))
        {
            assert_eq!(*temperature, source_t);
            close(
                answer.evidence.value.value * mineral.mass_kg_per_mol,
                source_cp,
                2e-13,
            );
            assert_eq!(
                answer.evidence.value.uncertainty,
                UncertaintyModel::Unstated
            );
            pack.claims().verify_receipt(&answer.receipt).unwrap();
        }
        let heating = pack
            .claims()
            .query_line_typed(
                &cp,
                "temperature",
                &point(mineral, 300.0),
                &point(mineral, 400.0),
                ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
            )
            .unwrap();
        // Authored 1.5 kg constant-pressure sensible-heating example. The
        // reference is the separately printed heat-content column, not a
        // trapezoidal calculation repeated from the Cp entries.
        let mass_kg = 1.5;
        let heat_j: f64 = mass_kg
            * heating
                .knots
                .windows(2)
                .map(|w| {
                    0.5 * (w[1].0 - w[0].0)
                        * (w[0].1.evidence.value.value + w[1].1.evidence.value.value)
                })
                .sum::<f64>();
        let expected = mass_kg / mineral.mass_kg_per_mol * mineral.source_heat_300_400_j_per_mol;
        close(heat_j, expected, mineral.heating_relative_tolerance);
        assert!(heat_j > 0.0);
        // The finite comparison tolerates the selected coarse linearization;
        // it is not a continuous error certificate or source uncertainty model.
    }
}

#[test]
fn mineral_references_refuse_other_polymorphs_conventions_and_temperatures() {
    let hf = PropertyKey::new("standard-molar-formation-enthalpy", MOLAR_ENERGY);
    let cp = PropertyKey::new("specific-heat-capacity", CP);
    for mineral in MINERALS {
        let pack = compile(mineral);
        let other_reference = if mineral.reference == JANAF {
            USGS
        } else {
            JANAF
        };
        for bad in [
            QueryPoint::new().with("temperature", 298.15).unwrap(),
            point(mineral, 300.0),
            point(mineral, 298.15).with(mineral.polymorph, 0.0).unwrap(),
            point(mineral, 298.15)
                .with("source-phase-solid", 0.0)
                .unwrap(),
            point(mineral, 298.15)
                .with("source-pure-species", 0.0)
                .unwrap(),
            point(mineral, 298.15)
                .with("source-specimen-qualified", 1.0)
                .unwrap(),
            point(mineral, 298.15)
                .with_quantity("pressure", semantic(QuantityKind::Pressure), 101325.0)
                .unwrap(),
            point(mineral, 298.15)
                .with(mineral.reference, 0.0)
                .unwrap()
                .with(other_reference, 1.0)
                .unwrap(),
        ] {
            assert!(
                pack.claims()
                    .query_typed(&hf, &bad, SelectionPolicy::SingleClaimOnly)
                    .is_err(),
                "{}",
                mineral.pack
            );
        }
        for temperature in [298.0, 501.0] {
            assert!(
                pack.claims()
                    .query_typed(
                        &cp,
                        &point(mineral, temperature),
                        SelectionPolicy::SingleClaimOnly
                    )
                    .is_err()
            );
        }
        if mineral.cp_molar.is_none() {
            assert!(
                pack.claims()
                    .query_typed(
                        &cp,
                        &point(mineral, 298.15),
                        SelectionPolicy::SingleClaimOnly
                    )
                    .is_err()
            );
        }
    }
}
