//! Source-backed calcium-sulfate endpoints, finite saturation and reagent accounts.
//!
//! Each equilibrium call follows ONE selected precipitation coordinate with
//! declared ideal free-ion activities. It does not solve sulfate protonation,
//! calcium-sulfate ion pairing, pH, rates or competing-phase stability.
#![deny(unsafe_code)]

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy};
use fs_qty::{
    Amount, Dimensionless, Dims, Mass, MolarMass, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    SpeciesId,
    stoichiometry::{
        ConservedReaction, FormationEnthalpySource,
        minerals::MineralReaction,
        thermodynamics::{
            FormationGibbsSource,
            equilibrium::{EquilibriumInventory, EquilibriumSolveOptions, IdealActivityModel},
        },
    },
};
use std::{
    error::Error,
    path::{Path, PathBuf},
};

type Result<T> = std::result::Result<T, Box<dyn Error>>;
const REFERENCE: &str = "source-elemental-reference-nbs1989";
const USAGE: &str = "usage: calcium_sulfate --materials PACK_DIR [--phase gypsum|anhydrite] [--calcium-mol N] [--sulfate-mol N] [--solid-mol N] [--water-mol N] [--activity-basis-mol N] [--sodium-mol N] [--chloride-mol N]";
const PACKS: [(&str, &str); 7] = [
    (
        "gypsum-selenite-solid-nbs1989-thermochemistry",
        "CaSO4.2H2O(s)",
    ),
    ("anhydrite-solid-nbs1989-thermochemistry", "CaSO4(s)"),
    ("calcium-ion-aqueous-nbs1989", "Ca[2+](aq)"),
    ("sulfate-ion-aqueous-nbs1989", "SO4[2-](aq)"),
    ("water-liquid-nbs1989-thermochemistry", "H2O(l)"),
    ("hydrogen-ion-aqueous-nbs1989", "H+(aq)"),
    (
        "calcium-hydroxide-solid-nbs1989-thermochemistry",
        "Ca(OH)2(s)",
    ),
];

fn point(index: usize) -> Result<QueryPoint> {
    let mut at = QueryPoint::new()
        .with_quantity(
            "temperature",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            298.15,
        )?
        .with_quantity(
            "pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            100_000.0,
        )?;
    let mut flags = vec![
        (REFERENCE, 1.0),
        ("source-standard-state", 1.0),
        ("engineering-reference-model", 1.0),
        ("source-specimen-qualified", 0.0),
    ];
    match index {
        0 | 1 | 6 => flags.extend([("source-phase-solid", 1.0), ("source-pure-species", 1.0)]),
        4 => flags.extend([("source-phase-liquid", 1.0), ("source-pure-species", 1.0)]),
        _ => flags.extend([
            ("source-phase-aqueous", 1.0),
            ("source-individual-ion", 1.0),
            ("source-aqueous-unit-activity-infinite-dilution", 1.0),
            ("source-ionic-reference-hydrogen-zero", 1.0),
        ]),
    }
    if index == 0 {
        flags.extend([
            ("source-phase-gypsum", 1.0),
            ("source-variety-selenite", 1.0),
        ]);
    } else if index == 1 {
        flags.push(("source-phase-anhydrite", 1.0));
    }
    for (axis, value) in flags {
        at = at.with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)?;
    }
    Ok(at)
}
struct Materials {
    packs: Vec<NormalizedPack>,
    points: Vec<QueryPoint>,
    ids: Vec<SpeciesId>,
}
impl Materials {
    fn load(directory: &Path) -> Result<Self> {
        let packs = PACKS
            .iter()
            .map(|(name, _)| -> Result<NormalizedPack> {
                let bytes = std::fs::read(directory.join(format!("{name}.fsmatpk")))?;
                let decoded = NormalizedPack::from_bytes(&bytes)?;
                let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes)?;
                if pack.pack_id() != *name || pack.to_bytes() != bytes {
                    return Err(format!("expected canonical material profile {name}").into());
                }
                println!("material={} hash={}", pack.pack_id(), pack.content_hash());
                Ok(pack)
            })
            .collect::<Result<Vec<_>>>()?;
        Ok(Self {
            packs,
            points: (0..PACKS.len()).map(point).collect::<Result<Vec<_>>>()?,
            ids: PACKS
                .iter()
                .map(|(_, id)| SpeciesId::new(*id))
                .collect::<std::result::Result<Vec<_>, _>>()?,
        })
    }
    fn index(&self, species: &SpeciesId) -> usize {
        self.ids
            .iter()
            .position(|s| s == species)
            .expect("declared reaction participant")
    }
    fn gibbs(&self, reaction: &ConservedReaction) -> Vec<FormationGibbsSource<'_>> {
        reaction
            .terms()
            .iter()
            .map(|term| {
                let i = self.index(&term.species);
                FormationGibbsSource {
                    species: &self.ids[i],
                    pack: &self.packs[i],
                    at: &self.points[i],
                    reference_axis: REFERENCE,
                }
            })
            .collect()
    }
    fn enthalpies(&self, reaction: &ConservedReaction) -> Vec<FormationEnthalpySource<'_>> {
        reaction
            .terms()
            .iter()
            .map(|term| {
                let i = self.index(&term.species);
                FormationEnthalpySource {
                    species: &self.ids[i],
                    pack: &self.packs[i],
                    at: &self.points[i],
                    reference_axis: REFERENCE,
                }
            })
            .collect()
    }
    fn molar_mass(&self, i: usize) -> Result<MolarMass> {
        let answer = self.packs[i].claims().query_typed(
            &PropertyKey::new("molar-mass", MolarMass::DIMS),
            &self.points[i],
            SelectionPolicy::SingleClaimOnly,
        )?;
        self.packs[i].claims().verify_receipt(&answer.receipt)?;
        if !answer.evidence.value.value.is_finite() || answer.evidence.value.value <= 0.0 {
            return Err("source molar mass must be positive and finite".into());
        }
        Ok(MolarMass::new(answer.evidence.value.value))
    }
}
#[derive(Clone)]
struct Scenario {
    phase: String,
    calcium: f64,
    sulfate: f64,
    solid: f64,
    water: f64,
    basis: f64,
    sodium: f64,
    chloride: f64,
}
impl Default for Scenario {
    fn default() -> Self {
        Self {
            phase: "gypsum".into(),
            calcium: 0.0,
            sulfate: 0.0,
            solid: 0.02,
            water: 55.5,
            basis: 1.0,
            sodium: 0.0,
            chloride: 0.0,
        }
    }
}
fn options() -> EquilibriumSolveOptions {
    EquilibriumSolveOptions {
        extent_tolerance: Amount::new(1e-12),
        log_residual_tolerance: Dimensionless::new(1e-9),
        max_iterations: 256,
    }
}
fn print_reaction(m: &Materials, recipe: MineralReaction) -> Result<()> {
    let reaction = recipe.build()?;
    let direct = reaction.standard_gibbs_change(&m.gibbs(&reaction))?;
    let caloric = reaction.standard_equilibrium(&m.enthalpies(&reaction))?;
    println!(
        "reaction={} delta_h_j_per_mol_extent={:.9e} direct_delta_g_j_per_mol_extent={:.9e} caloric_delta_g_j_per_mol_extent={:.9e} delta_s_j_per_mol_k={:.9e}",
        reaction.id().as_str(),
        caloric.enthalpy().molar().value(),
        direct.molar().value(),
        caloric.molar_gibbs().value(),
        caloric.molar_entropy().value(),
    );
    for term in direct.terms() {
        m.packs[m.index(&term.term.species)]
            .claims()
            .verify_receipt(&term.answer.receipt)?;
        println!(
            "gibbs_source species={} pack={} receipt={} uncertainty={:?}",
            term.term.species.as_str(),
            term.pack_id,
            term.answer.receipt.content_hash(),
            term.answer.evidence.value.uncertainty,
        );
    }
    Ok(())
}
fn run_case(m: &Materials, label: &str, s: &Scenario) -> Result<()> {
    let (recipe, solid_index) = match s.phase.as_str() {
        "gypsum" => (MineralReaction::GypsumPrecipitation, 0),
        "anhydrite" => (MineralReaction::AnhydritePrecipitation, 1),
        _ => return Err("phase must be gypsum or anhydrite".into()),
    };
    // Counterions are explicit fixed stocks outside this one reaction coordinate.
    // They provide no ion-pairing law or activity correction.
    let charge = 2.0 * s.calcium - 2.0 * s.sulfate + s.sodium - s.chloride;
    let charge_scale = 2.0 * s.calcium + 2.0 * s.sulfate + s.sodium + s.chloride;
    if !charge.is_finite()
        || !charge_scale.is_finite()
        || charge.abs() > 128.0 * f64::EPSILON * charge_scale
    {
        return Err(
            "initial free ions plus explicit fixed Na+/Cl- spectators must balance charge".into(),
        );
    }
    let reaction = recipe.build()?;
    let mut stocks = vec![(2, s.calcium), (3, s.sulfate), (solid_index, s.solid)];
    if solid_index == 0 {
        stocks.push((4, s.water));
    }
    let initial = stocks
        .iter()
        .map(|&(i, n)| EquilibriumInventory {
            species: m.ids[i].clone(),
            amount: Amount::new(n),
            activity: if matches!(i, 0 | 1 | 4) {
                IdealActivityModel::PureCondensed
            } else {
                IdealActivityModel::IdealDilute {
                    amount_at_unit_activity: Amount::new(s.basis),
                }
            },
        })
        .collect::<Vec<_>>();
    let answer = reaction.equilibrate_ideal(&m.gibbs(&reaction), &initial, options())?;
    let amount = |i: usize| -> f64 {
        answer
            .final_inventory
            .iter()
            .find(|(id, _)| id == &m.ids[i])
            .expect("complete final reaction inventory")
            .1
            .value()
    };
    let calcium = amount(2);
    let sulfate = amount(3);
    let solid = amount(solid_index);
    let calcium_residual = calcium + solid - s.calcium - s.solid;
    let sulfur_residual = sulfate + solid - s.sulfate - s.solid;
    let charge_residual = 2.0 * calcium - 2.0 * sulfate + s.sodium - s.chloride;
    let mut mass_residual = 0.0;
    for &(i, before) in &stocks {
        mass_residual += (amount(i) - before) * m.molar_mass(i)?.value();
    }
    let enthalpy = reaction.standard_enthalpy_change(&m.enthalpies(&reaction))?;
    let heat = enthalpy
        .at_extent(Amount::new(answer.extent.value().abs()))?
        .value()
        * answer.extent.value().signum();
    println!(
        "case={label} phase={} status={:?} extent_mol={:.12e} ln_q_over_k={:?} iterations={}",
        s.phase,
        answer.status,
        answer.extent.value(),
        answer.log_quotient_over_equilibrium.map(|r| r.value()),
        answer.iterations,
    );
    println!(
        "free_calcium_mol={calcium:.12e} free_sulfate_mol={sulfate:.12e} solid_mol={solid:.12e} solid_mass_kg={:.12e} activity_basis_mol={:.12e}",
        solid * m.molar_mass(solid_index)?.value(),
        s.basis,
    );
    println!(
        "fixed_sodium_mol={:.12e} fixed_chloride_mol={:.12e} calcium_balance_mol={calcium_residual:.12e} sulfate_balance_mol={sulfur_residual:.12e} charge_balance_mol_equiv={charge_residual:.12e} source_mass_balance_kg={mass_residual:.12e} signed_standard_heat_j={heat:.12e}",
        s.sodium, s.chloride,
    );
    if solid_index == 0 {
        println!(
            "water_mol={:.12e} water_plus_hydrate_balance_mol={:.12e}",
            amount(4),
            amount(4) + 2.0 * solid - s.water - 2.0 * s.solid,
        );
    }
    Ok(())
}
fn reagent_example(m: &Materials) -> Result<()> {
    let reaction = MineralReaction::GypsumPrecipitationWithHydratedLime.build()?;
    let reacted_lime = Mass::new(0.00370474);
    let extent = reaction.extent_from_reacted_mass(&m.ids[6], reacted_lime, m.molar_mass(6)?)?;
    let reactants = vec![
        (m.ids[6].clone(), Amount::new(0.1)),
        (m.ids[5].clone(), Amount::new(0.12)),
        (m.ids[3].clone(), Amount::new(0.1)),
    ];
    let mut initial = reactants.clone();
    initial.push((m.ids[0].clone(), Amount::new(0.0)));
    let final_inventory = reaction.advance_inventory(&initial, extent)?;
    let heat = reaction
        .standard_enthalpy_change(&m.enthalpies(&reaction))?
        .at_extent(extent)?;
    println!(
        "case=declared-lime-neutralization reacted_lime_kg={:.12e} extent_mol={:.12e} inventory_ceiling_mol={:.12e} standard_heat_j={:.12e}",
        reacted_lime.value(),
        extent.value(),
        reaction.limiting_extent(&reactants)?.value(),
        heat.value(),
    );
    for (species, amount) in final_inventory {
        println!(
            "remaining_species={} amount_mol={:.12e}",
            species.as_str(),
            amount.value()
        );
    }
    println!("neutralization_scope=declared_extent_not_a_pH_or_complete_sulfate_equilibrium");
    Ok(())
}
fn main() -> Result<()> {
    let mut args = std::env::args().skip(1);
    let mut directory: Option<PathBuf> = None;
    let mut custom = false;
    let mut s = Scenario::default();
    while let Some(key) = args.next() {
        if key == "--help" {
            println!("{USAGE}");
            println!("All values are nonnegative finite mol; the activity basis must be positive.");
            println!("Free-ion charge requires the explicit fixed Na+/Cl- spectator stocks.");
            println!("Required compiled material packs:");
            for (name, _) in PACKS {
                println!("  {name}.fsmatpk");
            }
            return Ok(());
        }
        let v = args
            .next()
            .ok_or_else(|| format!("missing value for {key}; {USAGE}"))?;
        if key == "--materials" {
            if directory.is_some() {
                return Err("duplicate --materials".into());
            }
            directory = Some(v.into());
            continue;
        }
        custom = true;
        if key == "--phase" {
            s.phase = v;
            continue;
        }
        let n: f64 = v.parse()?;
        if !n.is_finite() || n < 0.0 {
            return Err(format!("invalid nonnegative finite amount for {key}").into());
        }
        match key.as_str() {
            "--calcium-mol" => s.calcium = n,
            "--sulfate-mol" => s.sulfate = n,
            "--solid-mol" => s.solid = n,
            "--water-mol" => s.water = n,
            "--activity-basis-mol" if n > 0.0 => s.basis = n,
            "--activity-basis-mol" => return Err("activity basis must be positive".into()),
            "--sodium-mol" => s.sodium = n,
            "--chloride-mol" => s.chloride = n,
            _ => return Err(format!("unknown option {key}; {USAGE}").into()),
        }
    }
    let directory = directory.ok_or(USAGE)?;
    let m = Materials::load(&directory)?;
    println!(
        "model=one-selected-phase-ideal-free-ion-coordinate temperature_K=298.15 reference_pressure_Pa=100000"
    );
    println!(
        "scope=no_pH_no_sulfate_protonation_no_ion_pairing_no_kinetics_no_competing_phase_selection"
    );
    for recipe in [
        MineralReaction::GypsumPrecipitation,
        MineralReaction::AnhydritePrecipitation,
        MineralReaction::AnhydriteHydration,
        MineralReaction::GypsumPrecipitationWithHydratedLime,
    ] {
        print_reaction(&m, recipe)?;
    }
    if custom {
        run_case(&m, "custom", &s)?;
    } else {
        run_case(&m, "gypsum-dissolution", &s)?;
        run_case(
            &m,
            "gypsum-precipitation",
            &Scenario {
                calcium: 0.05,
                sulfate: 0.05,
                solid: 0.0,
                ..s.clone()
            },
        )?;
        run_case(
            &m,
            "finite-solid-exhaustion",
            &Scenario {
                solid: 0.001,
                ..s.clone()
            },
        )?;
        run_case(
            &m,
            "sulfate-common-ion",
            &Scenario {
                sulfate: 0.05,
                sodium: 0.1,
                ..s.clone()
            },
        )?;
        run_case(
            &m,
            "selected-anhydrite",
            &Scenario {
                phase: "anhydrite".into(),
                ..s
            },
        )?;
    }
    reagent_example(&m)?;
    Ok(())
}
