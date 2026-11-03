//! Actual USGS sources close copper electrowinning heat, Gibbs and Faraday bases.
//! Thermodynamic driving force does not supply cell voltage losses or kinetics.
#![deny(unsafe_code)]

use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicU64, Ordering};

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy, UncertaintyModel};
use fs_qty::{
    Amount, Dimensionless, Dims, ElectricCharge, MolarMass, QuantitySpec, Voltage,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    SpeciesId,
    equilibrium::EquilibriumConstantStatusV1,
    stoichiometry::{
        FARADAY_C_PER_MOL, FormationEnthalpySource, electrolysis::ElectrolysisReaction,
        minerals::MineralReaction, thermodynamics::FormationGibbsSource,
    },
};

const REFERENCE: &str = "source-elemental-reference-usgs-b2131";
const HF: &str = "standard-molar-formation-enthalpy";
const GF: &str = "standard-molar-formation-gibbs-energy";
const MOLAR_ENERGY: Dims = Dims([2, 1, -2, 0, 0, -1]);
const MOLAR_ENTROPY: Dims = Dims([2, 1, -2, -1, 0, -1]);
const SLUGS: [&str; 5] = [
    "copper-ii-ion-aqueous-usgs-b2131",
    "water-liquid-usgs-b2131",
    "copper-solid-usgs-b2131",
    "oxygen-gas-usgs-b2131",
    "hydrogen-ion-aqueous-usgs-b2131",
];
const SPECIES: [&str; 5] = ["Cu[2+](aq)", "H2O(l)", "Cu(s)", "O2(g)", "H+(aq)"];
const PHASES: [&str; 5] = ["aqueous", "liquid", "solid", "gas", "aqueous"];

fn workspace() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .to_path_buf()
}

fn compile(slug: &str) -> NormalizedPack {
    static NEXT: AtomicU64 = AtomicU64::new(0);
    let nonce = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let directory = std::env::temp_dir().join(format!(
        "frankensim-copper-gibbs-{}-{nonce}-{}",
        std::process::id(),
        NEXT.fetch_add(1, Ordering::Relaxed)
    ));
    std::fs::create_dir(&directory).unwrap();
    let output = directory.join("material.fsmatpk");
    let run = Command::new(env!("CARGO_BIN_EXE_xtask"))
        .args(["matdb-pack", "--manifest"])
        .arg(workspace().join(format!("data/matdb/seed-v1/{slug}/manifest.tsv")))
        .arg("--out")
        .arg(&output)
        .env("CARGO_WORKSPACE_DIR", workspace())
        .current_dir(workspace())
        .output()
        .expect("run the actual material compiler");
    assert!(
        run.status.success(),
        "{slug}: {}",
        String::from_utf8_lossy(&run.stderr)
    );
    let bytes = std::fs::read(output).unwrap();
    let decoded = NormalizedPack::from_bytes(&bytes).unwrap();
    let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes).unwrap();
    assert_eq!(pack.to_bytes(), bytes, "{slug}: canonical round trip");
    assert_eq!(pack.pack_id(), slug);
    pack
}

fn point(phase: &str) -> QueryPoint {
    let mut point = QueryPoint::new()
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
        ("source-standard-state", 1.0),
        (REFERENCE, 1.0),
        ("engineering-reference-model", 1.0),
        ("source-specimen-qualified", 0.0),
    ];
    if phase == "aqueous" {
        flags.extend([
            ("source-phase-aqueous", 1.0),
            ("source-aqueous-unit-activity-infinite-dilution", 1.0),
            ("source-ionic-reference-hydrogen-zero", 1.0),
            ("source-individual-ion", 1.0),
        ]);
    } else {
        let axis = match phase {
            "solid" => "source-phase-solid",
            "liquid" => "source-phase-liquid",
            "gas" => "source-phase-gas",
            _ => panic!("unknown test phase"),
        };
        flags.extend([(axis, 1.0), ("source-pure-species", 1.0)]);
    }
    for (axis, value) in flags {
        point = point
            .with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)
            .unwrap();
    }
    point
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 1e-12 * expected.abs().max(1.0),
        "{actual} != {expected}"
    );
}

fn assert_value(pack: &NormalizedPack, at: &QueryPoint, name: &str, dims: Dims, expected: f64) {
    let key = PropertyKey::new(name, dims);
    let answer = pack
        .claims()
        .query_typed(&key, at, SelectionPolicy::SingleClaimOnly)
        .unwrap();
    close(answer.evidence.value.value, expected);
    assert_eq!(answer.evidence.value.dims, dims);
    assert_eq!(
        answer.evidence.value.uncertainty,
        UncertaintyModel::Unstated
    );
    assert!(answer.receipt.observation_backed);
    assert!(!answer.receipt.source_hashes.is_empty());
    pack.claims().verify_receipt(&answer.receipt).unwrap();
    let claim = pack.claims().claim(answer.receipt.selected).unwrap();
    for (axis, value) in at.axes() {
        assert_eq!(claim.validity.bound(axis), Some((*value, *value)));
        assert_eq!(
            claim.validity.axis_quantities().get(axis),
            at.axis_quantities().get(axis)
        );
    }
    for wrong in [
        QueryPoint::new(),
        at.clone().with("temperature", 300.0).unwrap(),
        at.clone().with("pressure", 101_325.0).unwrap(),
        at.clone().with(REFERENCE, 0.0).unwrap(),
        at.clone().with("source-specimen-qualified", 1.0).unwrap(),
    ] {
        assert!(
            pack.claims()
                .query_typed(&key, &wrong, SelectionPolicy::SingleClaimOnly)
                .is_err()
        );
    }
}

#[test]
fn g0_g3_copper_sources_preserve_original_columns_and_reference_conventions() {
    // Original USGS B2131 p6: Cu2+ Hf=64.9, Gf=65.1 kJ/mol; Cu reference S=33.14.
    let ion = compile(SLUGS[0]);
    let ion_at = point("aqueous");
    assert_eq!(ion.claims().claim_count(), 6);
    for (name, dims, expected) in [
        ("molar-mass", MolarMass::DIMS, 0.063546),
        (HF, MOLAR_ENERGY, 64_900.0),
        (GF, MOLAR_ENERGY, 65_100.0),
        ("conventional-standard-molar-entropy", MOLAR_ENTROPY, -98.0),
        (
            "reported-formation-enthalpy-uncertainty-magnitude",
            MOLAR_ENERGY,
            1000.0,
        ),
        (
            "reported-conventional-entropy-uncertainty-magnitude",
            MOLAR_ENTROPY,
            4.0,
        ),
    ] {
        assert_value(&ion, &ion_at, name, dims, expected);
    }
    assert!(
        ion.claims()
            .query_typed(
                &PropertyKey::new("standard-molar-entropy", MOLAR_ENTROPY),
                &ion_at,
                SelectionPolicy::SingleClaimOnly
            )
            .is_err()
    );
    let solid = compile(SLUGS[2]);
    let solid_at = point("solid");
    assert_eq!(solid.claims().claim_count(), 5);
    for (name, dims, expected) in [
        ("molar-mass", MolarMass::DIMS, 0.063546),
        (HF, MOLAR_ENERGY, 0.0),
        (GF, MOLAR_ENERGY, 0.0),
        ("standard-molar-entropy", MOLAR_ENTROPY, 33.14),
        (
            "reported-standard-entropy-uncertainty-magnitude",
            MOLAR_ENTROPY,
            0.03,
        ),
    ] {
        assert_value(&solid, &solid_at, name, dims, expected);
    }
    for (pack, at, axis) in [
        (&ion, &ion_at, "source-ionic-reference-hydrogen-zero"),
        (&solid, &solid_at, "source-phase-solid"),
    ] {
        assert!(
            pack.claims()
                .query_typed(
                    &PropertyKey::new(GF, MOLAR_ENERGY),
                    &at.clone().with(axis, 0.0).unwrap(),
                    SelectionPolicy::SingleClaimOnly
                )
                .is_err()
        );
    }
}

#[test]
fn g1_real_sources_close_copper_electrowinning_heat_and_gibbs() {
    let packs = SLUGS.map(compile);
    let ids = SPECIES.map(|id| SpeciesId::new(id).unwrap());
    let points = PHASES.map(point);
    let enthalpy: [_; 5] = std::array::from_fn(|i| FormationEnthalpySource {
        species: &ids[i],
        pack: &packs[i],
        at: &points[i],
        reference_axis: REFERENCE,
    });
    let gibbs: [_; 5] = std::array::from_fn(|i| FormationGibbsSource {
        species: &ids[i],
        pack: &packs[i],
        at: &points[i],
        reference_axis: REFERENCE,
    });
    let reaction = MineralReaction::CopperElectrowinning.build().unwrap();
    let heat = reaction.standard_enthalpy_change(&enthalpy).unwrap();
    let free = reaction.standard_gibbs_change(&gibbs).unwrap();
    // 2 Cu2+ + 2 H2O(l) -> 2 Cu(s) + O2(g) + 4 H+.
    // Separate printed source columns: -2*64900+2*285800; -2*65100+2*237100.
    assert_eq!(heat.molar().value(), 441_800.0);
    assert_eq!(free.molar().value(), 344_000.0);
    close(free.log_equilibrium_constant().value(), -138.76803667515603);
    assert_eq!(free.status(), EquilibriumConstantStatusV1::FinitePositive);
    let k = free.equilibrium_constant().unwrap().value();
    assert!((k / 5.417605875520852e-61 - 1.0).abs() < 1e-10);
    assert_eq!(free.terms().len(), 5);
    for term in free.terms() {
        let pack = packs
            .iter()
            .find(|pack| pack.pack_id() == term.pack_id)
            .unwrap();
        pack.claims().verify_receipt(&term.answer.receipt).unwrap();
        assert!(!term.answer.receipt.source_hashes.is_empty());
    }
    assert!(
        reaction.standard_gibbs_change(&gibbs[..4]).is_err(),
        "zero-reference H+ remains required"
    );
    let without_copper = [0, 1, 3, 4].map(|i| FormationGibbsSource {
        species: &ids[i],
        pack: &packs[i],
        at: &points[i],
        reference_axis: REFERENCE,
    });
    assert!(
        reaction.standard_gibbs_change(&without_copper).is_err(),
        "zero-reference copper remains required"
    );

    let mut activities: [_; 5] = std::array::from_fn(|i| (ids[i].clone(), Dimensionless::new(1.0)));
    let standard_activity = free.at_activities(&activities).unwrap();
    assert_eq!(standard_activity.log_quotient.value(), 0.0);
    close(standard_activity.molar.value(), 344_000.0);
    // Only H+ activity changes to 0.01: Q=(0.01)^4. This is not a concentration model.
    activities[4].1 = Dimensionless::new(0.01);
    let state = free.at_activities(&activities).unwrap();
    close(state.log_quotient.value(), -18.420680743952364);
    close(state.molar.value(), 298_335.92397971795);
    assert!(
        state.molar.value() > 0.0,
        "productive direction still needs external work"
    );
    assert!(free.at_activities(&activities[..4]).is_err());
}

#[test]
fn g1_source_copper_mass_connects_four_electron_equation_to_finite_inventory() {
    let copper = compile(SLUGS[2]);
    let mass = copper
        .claims()
        .query_typed(
            &PropertyKey::new("molar-mass", MolarMass::DIMS),
            &point("solid"),
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
    let ids = SPECIES.map(|id| SpeciesId::new(id).unwrap());
    let reaction = MineralReaction::CopperElectrowinning.build().unwrap();
    // Cu(II) -> Cu(0) takes two electron-moles per mole Cu; the equation forms two Cu.
    // This independent electrode association is declared, not inferred from atom balance.
    let cell = ElectrolysisReaction::new(
        reaction,
        ids[2].clone(),
        MolarMass::new(mass.evidence.value.value),
        2,
    )
    .unwrap();
    assert_eq!(cell.electrons_per_product(), 2);
    assert_eq!(cell.electrons_per_equation(), 4);
    let initial: [_; 5] =
        std::array::from_fn(|i| (ids[i].clone(), Amount::new([2.0, 2.0, 0.0, 0.0, 0.0][i])));
    let result = cell
        .advance_charge(
            &initial,
            ElectricCharge::new(4.0 * FARADAY_C_PER_MOL),
            Dimensionless::new(1.0),
        )
        .unwrap();
    assert_eq!(result.extent.value(), 1.0);
    assert_eq!(result.product_yield.amount().value(), 2.0);
    close(result.product_yield.mass().value(), 0.127092);
    for (i, expected) in [0.0, 0.0, 2.0, 1.0, 4.0].into_iter().enumerate() {
        close(
            result
                .modeled_final_inventory
                .iter()
                .find(|(id, _)| id == &ids[i])
                .unwrap()
                .1
                .value(),
            expected,
        );
    }
    assert_eq!(result.unallocated_charge.value(), 0.0);
    assert!(
        cell.advance_charge(
            &initial,
            ElectricCharge::new(4.001 * FARADAY_C_PER_MOL),
            Dimensionless::new(1.0)
        )
        .is_err()
    );
    copper.claims().verify_receipt(&mass.receipt).unwrap();
}

#[test]
fn g1_source_bound_reversible_potential_retains_sign_activities_and_receipts() {
    let packs = SLUGS.map(compile);
    let ids = SPECIES.map(|id| SpeciesId::new(id).unwrap());
    let points = PHASES.map(point);
    let sources: [_; 5] = std::array::from_fn(|i| FormationGibbsSource {
        species: &ids[i],
        pack: &packs[i],
        at: &points[i],
        reference_axis: REFERENCE,
    });
    let copper_mass = packs[2]
        .claims()
        .query_typed(
            &PropertyKey::new("molar-mass", MolarMass::DIMS),
            &points[2],
            SelectionPolicy::SingleClaimOnly,
        )
        .unwrap();
    let model = ElectrolysisReaction::new(
        MineralReaction::CopperElectrowinning.build().unwrap(),
        ids[2].clone(),
        MolarMass::new(copper_mass.evidence.value.value),
        2,
    )
    .unwrap();
    let standard = model.standard_potential(&sources).unwrap();
    let signed: Voltage = standard.potential();
    // Independently evaluated 344000 J/mol for FOUR electron-moles per equation.
    // The opposing +0.891327 V magnitude is not a predicted operating voltage.
    close(signed.value(), -0.891327190438547);
    assert!(signed.value() < 0.0);
    assert_eq!(standard.model().electrons_per_equation(), 4);
    assert_eq!(standard.model().reaction().id(), model.reaction().id());
    assert_eq!(standard.gibbs().reaction(), model.reaction().id());
    assert_eq!(standard.gibbs().molar().value(), 344_000.0);
    assert_eq!(standard.gibbs().terms().len(), 5);
    for term in standard.gibbs().terms() {
        let pack = packs
            .iter()
            .find(|pack| pack.pack_id() == term.pack_id)
            .unwrap();
        pack.claims().verify_receipt(&term.answer.receipt).unwrap();
        assert!(!term.answer.receipt.source_hashes.is_empty());
    }
    let mut activities: [_; 5] = std::array::from_fn(|i| (ids[i].clone(), Dimensionless::new(1.0)));
    close(
        standard
            .at_activities(&activities)
            .unwrap()
            .potential()
            .value(),
        signed.value(),
    );
    // a(H+)=0.01, others unity: the same complete Gibbs state gives the Nernst shift.
    activities[4].1 = Dimensionless::new(0.01);
    let state = standard.at_activities(&activities).unwrap();
    close(state.potential().value(), -0.7730084910689824);
    close(state.gibbs().molar.value(), 298_335.92397971795);
    assert!(state.potential().value() > signed.value());
    assert!(std::ptr::eq(state.standard(), &standard));
    assert!(std::ptr::eq(state.gibbs().standard, standard.gibbs()));
    assert_eq!(state.gibbs().activities.len(), 5);
    assert!(standard.at_activities(&activities[..4]).is_err());
    for bad in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        activities[4].1 = Dimensionless::new(bad);
        assert!(standard.at_activities(&activities).is_err());
    }
    assert!(model.standard_potential(&sources[..4]).is_err());
    let without_copper = [0, 1, 3, 4].map(|i| FormationGibbsSource {
        species: &ids[i],
        pack: &packs[i],
        at: &points[i],
        reference_axis: REFERENCE,
    });
    assert!(model.standard_potential(&without_copper).is_err());
}
