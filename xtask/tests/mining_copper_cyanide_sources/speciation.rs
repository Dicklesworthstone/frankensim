//! Actual matched NBS rows through simultaneous Cu(I), ligand and pH closure.
use super::*;
use fs_thermochem::stoichiometry::{
    minerals::MineralReaction,
    thermodynamics::acid_base::{
        AcidDissociation, AcidFamily, IdealAcidBaseSystem, SpectatorIon, WaterAutoprotolysis,
        complexation::{
            IdealLigandSystem, LigandAssociation, LigandSpeciation, LigandSpeciationOptions,
            MetalLigandFamily,
        },
    },
};

const CYANIDE: &str = "uncomplexed-cyanide";
const COPPER: &str = "dissolved-copper-i";

fn acid_sources() -> &'static [Source; 4] {
    static SOURCES: OnceLock<[Source; 4]> = OnceLock::new();
    SOURCES.get_or_init(|| {
        let nonce = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let out = std::env::temp_dir().join(format!(
            "frankensim-copper-cyanide-acid-{}-{nonce}",
            std::process::id()
        ));
        std::fs::create_dir(&out).unwrap();
        let rows = [
            (
                "hydrogen-cyanide-aqueous-nbs1989",
                "HCN(aq)",
                "source-phase-aqueous",
                false,
            ),
            (
                "hydrogen-ion-aqueous-nbs1989",
                "H+(aq)",
                "source-phase-aqueous",
                true,
            ),
            (
                "hydroxide-ion-aqueous-nbs1989",
                "OH[-](aq)",
                "source-phase-aqueous",
                true,
            ),
            (
                "water-liquid-nbs1989-thermochemistry",
                "H2O(l)",
                "source-phase-liquid",
                false,
            ),
        ];
        std::array::from_fn(|i| {
            let (pack, id, phase, ionic) = rows[i];
            compile_source(pack, id, reference_point(phase, ionic), &out)
        })
    })
}

fn original(id: &str) -> &'static Source {
    sources()
        .iter()
        .chain(acid_sources())
        .find(|s| s.id.as_str() == id)
        .unwrap()
}

fn acid_gibbs(reaction: &ConservedReaction) -> Vec<FormationGibbsSource<'static>> {
    reaction
        .terms()
        .iter()
        .map(|term| {
            let s = original(term.species.as_str());
            FormationGibbsSource {
                species: &s.id,
                pack: &s.pack,
                at: &s.at,
                reference_axis: REFERENCE,
            }
        })
        .collect()
}

fn model() -> IdealLigandSystem {
    let water = MineralReaction::WaterDissociation.build().unwrap();
    let water_sources = acid_gibbs(&water);
    let water = WaterAutoprotolysis::from_sources(
        water,
        original("H2O(l)").id.clone(),
        original("H+(aq)").id.clone(),
        original("OH[-](aq)").id.clone(),
        &water_sources,
    )
    .unwrap();
    let acid = MineralReaction::HydrogenCyanideDissociation
        .build()
        .unwrap();
    let acid_sources = acid_gibbs(&acid);
    let acid = AcidDissociation::from_sources(
        acid,
        original("HCN(aq)").id.clone(),
        sources()[6].id.clone(),
        original("H+(aq)").id.clone(),
        &acid_sources,
    )
    .unwrap();
    let aqueous =
        IdealAcidBaseSystem::new(water, vec![AcidFamily::new(CYANIDE, vec![acid]).unwrap()])
            .unwrap();
    let steps = [
        (formation(2), 0, 3),
        (
            reaction(
                "copper-i-tricyanide-association",
                &[(3, -1), (6, -1), (4, 1)],
            )
            .unwrap(),
            3,
            4,
        ),
        (
            reaction(
                "copper-i-tetracyanide-association",
                &[(4, -1), (6, -1), (5, 1)],
            )
            .unwrap(),
            4,
            5,
        ),
    ]
    .into_iter()
    .map(|(reaction, parent, complex)| {
        let selected = gibbs(&reaction);
        LigandAssociation::from_sources(
            reaction,
            sources()[parent].id.clone(),
            sources()[6].id.clone(),
            sources()[complex].id.clone(),
            &selected,
        )
        .unwrap()
    })
    .collect();
    IdealLigandSystem::new(
        aqueous,
        CYANIDE,
        vec![MetalLigandFamily::new(COPPER, steps).unwrap()],
    )
    .unwrap()
}

fn spectators(ligand: f64, copper: f64) -> Vec<SpectatorIon> {
    // Explicit NaCN and dissolved Cu(I)/chloride component stocks; no mineral
    // accessibility, Cu(II) reduction or phase-stability claim is implied.
    vec![
        SpectatorIon {
            species: SpeciesId::new("Na[+](aq)").unwrap(),
            charge: 1,
            amount: Amount::new(ligand),
        },
        SpectatorIon {
            species: SpeciesId::new("Cl[-](aq)").unwrap(),
            charge: -1,
            amount: Amount::new(copper),
        },
    ]
}

fn solve(model: &IdealLigandSystem, ligand: f64, copper: f64, basis: f64) -> LigandSpeciation {
    model
        .speciate(
            Amount::new(ligand),
            &[(COPPER.into(), Amount::new(copper))],
            &[],
            &spectators(ligand, copper),
            Amount::new(basis),
            // Near ligand exhaustion pH is sensitive to the remaining ligand
            // balance error; tighten that control for the high-precision oracle.
            LigandSpeciationOptions {
                log_ligand_tolerance: Dimensionless::new(1e-12),
                absolute_ligand_tolerance: Amount::new(1e-16),
                relative_ligand_tolerance: Dimensionless::new(1e-13),
                ..LigandSpeciationOptions::default()
            },
            || false,
        )
        .unwrap()
}

fn relative(actual: f64, expected: f64) {
    // A relative check also protects the physically small Cu+ and tetra tails;
    // an absolute tolerance with a unit floor would accept their accidental loss.
    assert!(actual.is_finite() && actual > 0.0 && expected > 0.0);
    assert!(
        (actual / expected - 1.0).abs() < 2e-8,
        "{actual} != {expected}"
    );
}

fn balances(solution: &LigandSpeciation, ligand: f64, copper: f64) {
    let family = &solution.metals[0];
    assert_eq!(family.id, COPPER);
    close(family.total.value(), copper);
    close(
        family.members.iter().map(|m| m.amount.value()).sum(),
        copper,
    );
    let bound: f64 = family
        .members
        .iter()
        .map(|m| m.ligand_count as f64 * m.amount.value())
        .sum();
    assert!((solution.uncomplexed_ligand.value() + bound - ligand).abs() < 3e-12);
    close(bound, solution.bound_ligand.value());
    let charge = solution.aqueous.proton_amount.value()
        - solution.aqueous.hydroxide_amount.value()
        - solution.free_ligand.value()
        + ligand
        - copper
        + family
            .members
            .iter()
            .map(|m| m.charge as f64 * m.amount.value())
            .sum::<f64>();
    assert!(charge.abs() < 3e-12, "charge residual {charge}");
    assert!(solution.aqueous.charge_balance_residual.value().abs() < 3e-12);
    assert!(solution.ligand_balance_residual.value().abs() < 3e-12);
    assert!(solution.max_log_mass_action_residual.value() < 1e-8);
}

#[test]
fn g3_actual_nbs_copper_speciation_matches_independent_finite_balance_reference() {
    let solution = solve(&model(), 0.01, 0.002, 10.0);
    // Independent 70-digit Decimal reference from the literal Gibbs rows:
    // Ka=exp(-52700/RT), Kw=exp(-79885/RT), beta_n=exp(-DeltaG_n/RT).
    // For the stated counterions, exact balances reduce to
    // h^2=Kw/(1+c/Ka), L/n_ref=c+c*h/Ka+(Cu/n_ref)*sum(n*beta_n*c^n)/sum(beta_n*c^n).
    // Solve only that scalar equation; no production solver supplies this oracle.
    relative(solution.free_ligand.value(), 0.003_341_286_961_580_165);
    relative(solution.bound_ligand.value(), 0.005_898_993_082_565_126);
    assert!((solution.aqueous.ph.value() - 9.875_900_792_902_809).abs() < 2e-9);
    let hcn = solution.aqueous.families[0]
        .forms
        .iter()
        .find(|s| s.species.as_str() == "HCN(aq)")
        .unwrap();
    relative(hcn.amount.value(), 0.000_759_719_955_854_709_2);
    for (member, (id, count, expected)) in solution.metals[0].members.iter().zip([
        ("Cu[+](aq)", 0, 1.175_205_313_948_620_7e-21),
        ("Cu(CN)2[-](aq)", 2, 0.000_130_556_661_157_821_64),
        ("Cu(CN)3[2-](aq)", 3, 0.001_839_893_595_119_230_6),
        ("Cu(CN)4[3-](aq)", 4, 0.000_029_549_743_722_947_72),
    ]) {
        assert_eq!(member.species.as_str(), id);
        assert_eq!(member.ligand_count, count);
        assert_eq!(member.charge, 1 - count as i32);
        relative(member.amount.value(), expected);
        relative(member.activity.value(), expected / 10.0);
    }
    // The returned coupled solution keeps original source receipts, not just Ks.
    for step in solution.model().families()[0].steps() {
        for term in step.standard().terms() {
            original(term.term.species.as_str())
                .pack
                .claims()
                .verify_receipt(&term.answer.receipt)
                .unwrap();
        }
    }
    balances(&solution, 0.01, 0.002);
}

#[test]
fn g0_zero_copper_reduces_to_original_hcn_water_and_zero_ligand_retains_copper() {
    let model = model();
    let empty_copper = solve(&model, 0.01, 0.0, 10.0);
    let original = model
        .aqueous_model()
        .speciate(
            &[(CYANIDE.into(), Amount::new(0.01))],
            &spectators(0.01, 0.0),
            Amount::new(10.0),
            LigandSpeciationOptions::default().aqueous,
            || false,
        )
        .unwrap();
    close(empty_copper.aqueous.ph.value(), original.ph.value());
    relative(empty_copper.free_ligand.value(), 0.008_769_228_118_621_324);
    assert_eq!(empty_copper.bound_ligand.value(), 0.0);
    assert!(
        empty_copper.metals[0]
            .members
            .iter()
            .all(|m| m.amount.value() == 0.0)
    );
    let competing = solve(&model, 0.01, 0.002, 10.0);
    assert!(competing.free_ligand.value() < 0.4 * empty_copper.free_ligand.value());
    assert!(competing.aqueous.ph.value() < empty_copper.aqueous.ph.value());
    let empty_ligand = solve(&model, 0.0, 0.002, 10.0);
    assert_eq!(empty_ligand.free_ligand.value(), 0.0);
    assert_eq!(empty_ligand.bound_ligand.value(), 0.0);
    close(empty_ligand.metals[0].members[0].amount.value(), 0.002);
    assert!(
        empty_ligand.metals[0].members[1..]
            .iter()
            .all(|m| m.amount.value() == 0.0)
    );
    balances(&empty_copper, 0.01, 0.0);
    balances(&empty_ligand, 0.0, 0.002);
}

#[test]
fn g0_ligand_scarcity_preserves_free_copper_and_positive_small_complex_tails() {
    let solution = solve(&model(), 0.01, 0.006, 10.0);
    relative(solution.free_ligand.value(), 2.241_586_736_042_033_4e-11);
    relative(
        solution.metals[0].members[0].amount.value(),
        0.001_000_002_169_388_901_6,
    );
    relative(
        solution.metals[0].members[2].amount.value(),
        4.727_213_525_954_469e-10,
    );
    relative(
        solution.metals[0].members[3].amount.value(),
        5.093_402_567_587_196e-20,
    );
    assert!(
        (solution.aqueous.ph.value() - 6.998_453_329_150_441).abs() < 2e-9,
        "pH={} charge={} ligand={} inner={:?}",
        solution.aqueous.ph.value(),
        solution.aqueous.charge_balance_residual.value(),
        solution.ligand_balance_residual.value(),
        solution.aqueous.ph_bracket,
    );
    assert!(solution.bound_ligand.value() < 0.01);
    assert!(solution.metals[0].members[1].amount.value() < 0.005);
    balances(&solution, 0.01, 0.006);
}

#[test]
fn g3_common_molality_basis_scales_amounts_and_changes_coordination_when_changed_alone() {
    let model = model();
    let first = solve(&model, 0.01, 0.002, 10.0);
    let scaled = solve(&model, 0.07, 0.014, 70.0);
    assert!((scaled.aqueous.ph.value() - first.aqueous.ph.value()).abs() < 2e-9);
    relative(scaled.free_ligand.value(), 7.0 * first.free_ligand.value());
    for (a, b) in first.metals[0]
        .members
        .iter()
        .zip(&scaled.metals[0].members)
    {
        relative(b.amount.value(), 7.0 * a.amount.value());
        relative(b.activity.value(), a.activity.value());
    }
    // Holding amounts fixed while changing solvent activity amount is a
    // different chemical state, even though a separate process volume may agree.
    let concentrated = solve(&model, 0.01, 0.002, 1.0);
    relative(
        concentrated.free_ligand.value(),
        0.003_481_336_222_621_698_6,
    );
    relative(
        concentrated.metals[0].members[3].amount.value(),
        0.000_285_036_446_771_160_77,
    );
    assert!((concentrated.aqueous.ph.value() - 10.384_816_545_483_561).abs() < 2e-9);
    assert!(
        concentrated.metals[0].members[3].amount.value()
            > 9.0 * first.metals[0].members[3].amount.value()
    );
    balances(&scaled, 0.07, 0.014);
    balances(&concentrated, 0.01, 0.002);
}
