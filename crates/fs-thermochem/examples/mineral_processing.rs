//! Evaluate mineral reaction amounts and heat from actual material packs.
//!
//! Usage: `cargo run -p fs-thermochem --example mineral_processing -- \
//! cao.fsmatpk water.fsmatpk caoh2.fsmatpk 56.077`
//! Select another endpoint by prefixing `--reaction magnesite-calcination` or
//! `--reaction magnetite-oxidation`; pack arguments must match that equation's
//! explicit source identities. Copper oxide reduction and calcite bicarbonate
//! neutralization modes take four packs; dolomite bicarbonate neutralization
//! takes five. The carbonate endpoints use compatible USGS aqueous-ion
//! reference data, with bicarbonate retained in solution.
//! `--help` lists every supported reaction and required pack order.
//! Prefix a source-pack invocation with `--equilibrium` to also resolve standard
//! entropy, Gibbs energy and log equilibrium constants at the same source state.
//! Or explicit malachite inventory accounting, with all four inputs in moles:
//! `cargo run -p fs-thermochem --example mineral_processing -- \
//! --malachite-inventory 10 12 100 3`
//!
//! In the source-pack modes, the final argument is kilograms of the first pure
//! reactant actually reacted. It is not commercial feed mass; feed purity and
//! conversion require separate evidence.
//! This is standard-state bookkeeping at 298.15 K and 1 bar, not a plant
//! temperature, reaction-rate, water-dosing or recovery prediction.

use std::{error::Error, path::Path};

use fs_matdb::{MatDbError, NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Amount, Dims, Mass, MolarMass, QuantitySpec};
use fs_thermochem::{
    MassAmountBasis, SpeciesId,
    stoichiometry::{FormationEnthalpySource, minerals::MineralReaction},
};

const REFERENCE: &str = "source-elemental-reference-janaf";
const USGS_REFERENCE: &str = "source-elemental-reference-usgs-b2131";

struct ReactionCase {
    reaction: MineralReaction,
    reference: &'static str,
    packs: &'static [&'static str],
    species: &'static [&'static str],
    phase_axes: &'static [&'static str],
    equation: &'static str,
}

const REACTIONS: [&str; 13] = [
    "lime-hydration",
    "magnesite-calcination",
    "magnetite-oxidation",
    "cuprous-oxide-oxidation",
    "cupric-oxide-reduction",
    "cuprous-oxide-reduction",
    "pyrite-roasting",
    "pyrite-oxidation-ferrous",
    "sulfur-dioxide-oxidation",
    "sulfur-trioxide-hydration",
    "pyrite-to-acid",
    "calcite-bicarbonate-neutralization",
    "dolomite-bicarbonate-neutralization",
];
const USAGE: &str = "usage: mineral_processing [--equilibrium] [--reaction NAME] PACK... REACTED_PURE_REACTANT_KG\n       mineral_processing --malachite-inventory MALACHITE_MOL PROTON_MOL WATER_MOL EXTENT_MOL";

fn case(name: &str) -> Result<ReactionCase, Box<dyn Error>> {
    Ok(match name {
        "lime-hydration" => ReactionCase {
            reaction: MineralReaction::LimeHydration,
            reference: REFERENCE,
            packs: &[
                "calcium-oxide-solid-nist-mining",
                "water-liquid-janaf-mining",
                "calcium-hydroxide-solid-nist-mining",
            ],
            species: &["CaO(s)", "H2O(l)", "Ca(OH)2(s)"],
            phase_axes: &[
                "source-phase-solid",
                "source-phase-liquid",
                "source-phase-solid",
            ],
            equation: "CaO(s) + H2O(l) -> Ca(OH)2(s)",
        },
        "magnesite-calcination" => ReactionCase {
            reaction: MineralReaction::MagnesiteCalcination,
            reference: REFERENCE,
            packs: &[
                "magnesite-solid-janaf-mining",
                "magnesium-oxide-solid-janaf-mining",
                "carbon-dioxide-gas-janaf-mining",
            ],
            species: &["MgCO3(s)", "MgO(s)", "CO2(g)"],
            phase_axes: &[
                "source-phase-solid",
                "source-phase-solid",
                "source-phase-gas",
            ],
            equation: "MgCO3(s) -> MgO(s) + CO2(g)",
        },
        "magnetite-oxidation" => ReactionCase {
            reaction: MineralReaction::MagnetiteOxidationToHematite,
            reference: REFERENCE,
            packs: &[
                "magnetite-solid-janaf-mining",
                "oxygen-gas-janaf-mining",
                "hematite-solid-janaf-mining",
            ],
            species: &["Fe3O4(s)", "O2(g)", "Fe2O3(s)"],
            phase_axes: &[
                "source-phase-solid",
                "source-phase-gas",
                "source-phase-solid",
            ],
            equation: "4 Fe3O4(s) + O2(g) -> 6 Fe2O3(s), hematite endpoint",
        },
        "cuprous-oxide-oxidation" => ReactionCase {
            reaction: MineralReaction::CuprousOxideOxidation,
            reference: REFERENCE,
            packs: &[
                "copper-i-oxide-solid-janaf-mining",
                "oxygen-gas-janaf-mining",
                "copper-ii-oxide-solid-janaf-mining",
            ],
            species: &["Cu2O(s)", "O2(g)", "CuO(s)"],
            phase_axes: &[
                "source-phase-solid",
                "source-phase-gas",
                "source-phase-solid",
            ],
            equation: "2 Cu2O(s) + O2(g) -> 4 CuO(s)",
        },
        "cupric-oxide-reduction" => ReactionCase {
            reaction: MineralReaction::CupricOxideCarbonMonoxideReduction,
            reference: REFERENCE,
            packs: &[
                "copper-ii-oxide-solid-janaf-mining",
                "carbon-monoxide-gas-janaf-mining",
                "copper-solid-janaf-mining",
                "carbon-dioxide-gas-janaf-mining",
            ],
            species: &["CuO(s)", "CO(g)", "Cu(s)", "CO2(g)"],
            phase_axes: &[
                "source-phase-solid",
                "source-phase-gas",
                "source-phase-solid",
                "source-phase-gas",
            ],
            equation: "CuO(s) + CO(g) -> Cu(s) + CO2(g)",
        },
        "cuprous-oxide-reduction" => ReactionCase {
            reaction: MineralReaction::CuprousOxideCarbonMonoxideReduction,
            reference: REFERENCE,
            packs: &[
                "copper-i-oxide-solid-janaf-mining",
                "carbon-monoxide-gas-janaf-mining",
                "copper-solid-janaf-mining",
                "carbon-dioxide-gas-janaf-mining",
            ],
            species: &["Cu2O(s)", "CO(g)", "Cu(s)", "CO2(g)"],
            phase_axes: &[
                "source-phase-solid",
                "source-phase-gas",
                "source-phase-solid",
                "source-phase-gas",
            ],
            equation: "Cu2O(s) + CO(g) -> 2 Cu(s) + CO2(g)",
        },
        "pyrite-roasting" => ReactionCase {
            reaction: MineralReaction::PyriteRoastingToHematite,
            reference: REFERENCE,
            packs: &[
                "pyrite-solid-janaf-thermochemistry",
                "oxygen-gas-janaf-mining",
                "hematite-solid-janaf-mining",
                "sulfur-dioxide-gas-janaf-mining",
            ],
            species: &["FeS2(s)", "O2(g)", "Fe2O3(s)", "SO2(g)"],
            phase_axes: &[
                "source-phase-solid",
                "source-phase-gas",
                "source-phase-solid",
                "source-phase-gas",
            ],
            equation: "4 FeS2(s) + 11 O2(g) -> 2 Fe2O3(s) + 8 SO2(g), hematite endpoint",
        },
        "pyrite-oxidation-ferrous" => ReactionCase {
            reaction: MineralReaction::PyriteOxidationFerrous,
            reference: USGS_REFERENCE,
            packs: &[
                "pyrite-solid-usgs-b2131-thermochemistry",
                "oxygen-gas-usgs-b2131",
                "water-liquid-usgs-b2131",
                "iron-ii-ion-aqueous-usgs-b2131",
                "sulfate-ion-aqueous-usgs-b2131",
                "hydrogen-ion-aqueous-usgs-b2131",
            ],
            species: &[
                "FeS2(s)",
                "O2(g)",
                "H2O(l)",
                "Fe[2+](aq)",
                "SO4[2-](aq)",
                "H+(aq)",
            ],
            phase_axes: &[
                "source-phase-solid",
                "source-phase-gas",
                "source-phase-liquid",
                "source-phase-aqueous",
                "source-phase-aqueous",
                "source-phase-aqueous",
            ],
            equation: "2 FeS2(s) + 7 O2(g) + 2 H2O(l) -> 2 Fe[2+](aq) + 4 SO4[2-](aq) + 4 H+(aq), ferrous acid-generation endpoint",
        },
        "sulfur-dioxide-oxidation" => ReactionCase {
            reaction: MineralReaction::SulfurDioxideOxidation,
            reference: REFERENCE,
            packs: &[
                "sulfur-dioxide-gas-janaf-mining",
                "oxygen-gas-janaf-mining",
                "sulfur-trioxide-gas-janaf-mining",
            ],
            species: &["SO2(g)", "O2(g)", "SO3(g)"],
            phase_axes: &["source-phase-gas", "source-phase-gas", "source-phase-gas"],
            equation: "2 SO2(g) + O2(g) -> 2 SO3(g)",
        },
        "sulfur-trioxide-hydration" => ReactionCase {
            reaction: MineralReaction::SulfurTrioxideHydration,
            reference: REFERENCE,
            packs: &[
                "sulfur-trioxide-gas-janaf-mining",
                "water-liquid-janaf-mining",
                "sulfuric-acid-liquid-janaf-mining",
            ],
            species: &["SO3(g)", "H2O(l)", "H2SO4(l)"],
            phase_axes: &[
                "source-phase-gas",
                "source-phase-liquid",
                "source-phase-liquid",
            ],
            equation: "SO3(g) + H2O(l) -> H2SO4(l), pure liquid acid endpoint",
        },
        "pyrite-to-acid" => ReactionCase {
            reaction: MineralReaction::PyriteToHematiteAndSulfuricAcid,
            reference: REFERENCE,
            packs: &[
                "pyrite-solid-janaf-thermochemistry",
                "oxygen-gas-janaf-mining",
                "water-liquid-janaf-mining",
                "hematite-solid-janaf-mining",
                "sulfuric-acid-liquid-janaf-mining",
            ],
            species: &["FeS2(s)", "O2(g)", "H2O(l)", "Fe2O3(s)", "H2SO4(l)"],
            phase_axes: &[
                "source-phase-solid",
                "source-phase-gas",
                "source-phase-liquid",
                "source-phase-solid",
                "source-phase-liquid",
            ],
            equation: "4 FeS2(s) + 15 O2(g) + 8 H2O(l) -> 2 Fe2O3(s) + 8 H2SO4(l), all sulfur assigned to pure liquid acid",
        },
        "calcite-bicarbonate-neutralization" => ReactionCase {
            reaction: MineralReaction::CalciteBicarbonateNeutralization,
            reference: USGS_REFERENCE,
            packs: &[
                "calcite-solid-usgs-b2131-thermochemistry",
                "hydrogen-ion-aqueous-usgs-b2131",
                "calcium-ion-aqueous-usgs-b2131",
                "bicarbonate-ion-aqueous-usgs-b2131",
            ],
            species: &["CaCO3(s)", "H+(aq)", "Ca[2+](aq)", "HCO3[-](aq)"],
            phase_axes: &[
                "source-phase-solid",
                "source-phase-aqueous",
                "source-phase-aqueous",
                "source-phase-aqueous",
            ],
            equation: "CaCO3(calcite) + H+(aq) -> Ca[2+](aq) + HCO3[-](aq), bicarbonate endpoint",
        },
        "dolomite-bicarbonate-neutralization" => ReactionCase {
            reaction: MineralReaction::DolomiteBicarbonateNeutralization,
            reference: USGS_REFERENCE,
            packs: &[
                "dolomite-solid-usgs-b2131-thermochemistry",
                "hydrogen-ion-aqueous-usgs-b2131",
                "calcium-ion-aqueous-usgs-b2131",
                "magnesium-ion-aqueous-usgs-b2131",
                "bicarbonate-ion-aqueous-usgs-b2131",
            ],
            species: &[
                "CaMg(CO3)2(s)",
                "H+(aq)",
                "Ca[2+](aq)",
                "Mg[2+](aq)",
                "HCO3[-](aq)",
            ],
            phase_axes: &[
                "source-phase-solid",
                "source-phase-aqueous",
                "source-phase-aqueous",
                "source-phase-aqueous",
                "source-phase-aqueous",
            ],
            equation: "CaMg(CO3)2(dolomite) + 2 H+(aq) -> Ca[2+](aq) + Mg[2+](aq) + 2 HCO3[-](aq), bicarbonate endpoint",
        },
        _ => return Err(format!("unknown reaction {name:?}; use --help for source order").into()),
    })
}

fn load(path: &Path, expected_id: &str) -> Result<NormalizedPack, Box<dyn Error>> {
    let bytes = std::fs::read(path)?;
    let decoded = NormalizedPack::from_bytes(&bytes)?;
    let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes)?;
    if pack.to_bytes() != bytes {
        return Err(format!("{} is not a canonical material pack", path.display()).into());
    }
    // This example explicitly associates known source profiles with phases;
    // arbitrary property packs do not authenticate a chemical formula.
    if pack.pack_id() != expected_id {
        return Err(format!("expected pack {expected_id}, received {}", pack.pack_id()).into());
    }
    Ok(pack)
}

fn reference_point(
    reference: &str,
    phase_axis: &str,
    pack_id: &str,
) -> Result<QueryPoint, MatDbError> {
    let mut point = QueryPoint::new()
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
    for (axis, value) in [
        ("source-standard-state", 1.0),
        (reference, 1.0),
        ("engineering-reference-model", 1.0),
        ("source-specimen-qualified", 0.0),
        (phase_axis, 1.0),
    ] {
        point = point.with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)?;
    }
    // Aqueous ions use the source's hypothetical unit-activity, infinite-
    // dilution standard state. They are not pure isolated-ion solutions.
    let state_axes: &[&str] = if phase_axis == "source-phase-aqueous" {
        &[
            "source-individual-ion",
            "source-aqueous-unit-activity-infinite-dilution",
            "source-ionic-reference-hydrogen-zero",
        ]
    } else {
        &["source-pure-species"]
    };
    for axis in state_axes {
        point = point.with_quantity(*axis, QuantitySpec::dimensional(Dims::NONE), 1.0)?;
    }
    if pack_id == "calcite-solid-usgs-b2131-thermochemistry" {
        point = point.with_quantity(
            "source-phase-calcite",
            QuantitySpec::dimensional(Dims::NONE),
            1.0,
        )?;
    }
    if pack_id == "dolomite-solid-usgs-b2131-thermochemistry" {
        point = point.with_quantity(
            "source-phase-dolomite",
            QuantitySpec::dimensional(Dims::NONE),
            1.0,
        )?;
    }
    if matches!(
        pack_id,
        "pyrite-solid-janaf-thermochemistry" | "pyrite-solid-usgs-b2131-thermochemistry"
    ) {
        point = point.with_quantity(
            "source-phase-pyrite",
            QuantitySpec::dimensional(Dims::NONE),
            1.0,
        )?;
    }
    Ok(point)
}

fn malachite_inventory(args: &[String]) -> Result<(), Box<dyn Error>> {
    if args.len() != 5 {
        return Err("usage: mineral_processing --malachite-inventory MALACHITE_MOL PROTON_MOL WATER_MOL EXTENT_MOL".into());
    }
    let reaction = MineralReaction::MalachiteAcidLeach.build()?;
    let initial = [
        (
            SpeciesId::new("Cu2CO3(OH)2(s)")?,
            Amount::new(args[1].parse()?),
        ),
        (SpeciesId::new("H+(aq)")?, Amount::new(args[2].parse()?)),
        (SpeciesId::new("H2O(l)")?, Amount::new(args[3].parse()?)),
        (SpeciesId::new("Cu[2+](aq)")?, Amount::new(0.0)),
        (SpeciesId::new("CO2(aq)")?, Amount::new(0.0)),
    ];
    let extent = Amount::new(args[4].parse()?);
    let bound = reaction.limiting_extent(&initial[..2])?;
    let final_inventory = reaction.advance_inventory(&initial, extent)?;
    println!("Cu2CO3(OH)2(s) + 4 H+(aq) -> 2 Cu[2+](aq) + CO2(aq) + 3 H2O(l)");
    println!(
        "material-availability extent bound [mol]: {}",
        bound.value()
    );
    println!("caller-declared extent [mol]: {}", extent.value());
    for (species, remaining) in final_inventory {
        let before = initial
            .iter()
            .find(|(name, _)| *name == species)
            .ok_or("final species is absent from initial inventory")?
            .1;
        println!(
            "{} [mol]: {} -> {}",
            species.as_str(),
            before.value(),
            remaining.value()
        );
    }
    println!(
        "Extent is an input, not predicted recovery. Proton inventory is not a commercial acid dose."
    );
    Ok(())
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.as_slice() == ["--help"] {
        println!("{USAGE}");
        println!("Default without --reaction: lime-hydration. Source pack order:");
        for name in REACTIONS {
            let selected = case(name)?;
            println!("{name}: {}", selected.packs.join(" "));
        }
        return Ok(());
    }
    let show_equilibrium = args.first().is_some_and(|arg| arg == "--equilibrium");
    let args = if show_equilibrium {
        &args[1..]
    } else {
        args.as_slice()
    };
    if args
        .first()
        .is_some_and(|arg| arg == "--malachite-inventory")
    {
        if show_equilibrium {
            return Err("--equilibrium requires a complete source-pack reaction mode".into());
        }
        return malachite_inventory(args);
    }
    let (case, args) = if args.len() >= 2 && args[0] == "--reaction" {
        (case(&args[1])?, &args[2..])
    } else if args.len() == 4 {
        (case("lime-hydration")?, args)
    } else {
        return Err(USAGE.into());
    };
    if args.len() != case.packs.len() + 1 {
        return Err(format!(
            "expected {} source packs in this order, then reacted kg: {}",
            case.packs.len(),
            case.packs.join(" ")
        )
        .into());
    }
    let reacted_kg: f64 = args[case.packs.len()].parse()?;
    let packs = args
        .iter()
        .zip(case.packs)
        .map(|(path, expected)| load(Path::new(path), expected))
        .collect::<Result<Vec<_>, _>>()?;
    let species = case
        .species
        .iter()
        .map(|id| SpeciesId::new(*id))
        .collect::<Result<Vec<_>, _>>()?;
    let points = case
        .phase_axes
        .iter()
        .zip(case.packs)
        .map(|(phase, pack)| reference_point(case.reference, phase, pack))
        .collect::<Result<Vec<_>, _>>()?;
    let sources = (0..packs.len())
        .map(|i| FormationEnthalpySource {
            species: &species[i],
            pack: &packs[i],
            at: &points[i],
            reference_axis: case.reference,
        })
        .collect::<Vec<_>>();
    let reaction = case.reaction.build()?;
    let equilibrium = if show_equilibrium {
        Some(reaction.standard_equilibrium(&sources)?)
    } else {
        None
    };
    let heat = match &equilibrium {
        Some(standard) => standard.enthalpy().clone(),
        None => reaction.standard_enthalpy_change(&sources)?,
    };
    let mut molar_masses = Vec::with_capacity(packs.len());
    for (pack, point) in packs.iter().zip(&points) {
        let answer = pack.claims().query_typed(
            &PropertyKey::with_quantity(
                "molar-mass",
                QuantitySpec::dimensional(Dims([0, 1, 0, 0, 0, -1])),
            ),
            point,
            SelectionPolicy::SingleClaimOnly,
        )?;
        pack.claims().verify_receipt(&answer.receipt)?;
        molar_masses.push(MolarMass::new(answer.evidence.value.value));
    }
    let extent =
        reaction.extent_from_reacted_mass(&species[0], Mass::new(reacted_kg), molar_masses[0])?;
    let amounts = reaction.amounts(extent)?;

    println!("{}, 298.15 K, 100000 Pa", case.equation);
    println!("reacted pure {} mass [kg]: {reacted_kg}", case.species[0]);
    println!(
        "source reactant molar mass [kg/mol]: {}",
        molar_masses[0].value()
    );
    println!(
        "reaction extent [mol of written equation]: {}",
        extent.value()
    );
    let mut consumed_mass = 0.0;
    let mut formed_mass = 0.0;
    for amount in &amounts {
        let index = species
            .iter()
            .position(|id| id == &amount.term.species)
            .ok_or("reaction participant lacks an associated source")?;
        let mass = MassAmountBasis::from_amount(amount.amount, molar_masses[index])?;
        let direction = if amount.term.coefficient < 0 {
            consumed_mass += mass.mass().value();
            "consumed"
        } else {
            formed_mass += mass.mass().value();
            "formed"
        };
        println!(
            "{} {direction}: {} mol, {} kg",
            case.species[index],
            amount.amount.value(),
            mass.mass().value()
        );
    }
    if !consumed_mass.is_finite() || !formed_mass.is_finite() {
        return Err("summed source-molar-mass conversion is not representable".into());
    }
    println!(
        "Source-molar-mass residual, formed minus consumed [kg]: {}",
        formed_mass - consumed_mass
    );
    println!(
        "Atoms and charge close exactly; independently rounded source molar masses are retained without adjustment."
    );
    println!(
        "standard reaction enthalpy [J/mol extent]: {}",
        heat.molar().value()
    );
    println!(
        "standard enthalpy change [J]: {}",
        heat.at_extent(extent)?.value()
    );
    if let Some(standard) = equilibrium {
        println!(
            "standard reaction entropy [J/(mol extent K)]: {}",
            standard.molar_entropy().value()
        );
        println!(
            "standard reaction Gibbs energy [J/mol extent]: {}",
            standard.molar_gibbs().value()
        );
        println!(
            "ln K (dimensionless standard activities): {}",
            standard.log_equilibrium_constant().value()
        );
        println!(
            "log10 K (dimensionless standard activities): {}",
            standard.log_equilibrium_constant().value() / core::f64::consts::LN_10
        );
        println!(
            "These constants use the selected standard states; activities, pH and equilibrium composition require separate inputs and models."
        );
    }
    for pack in &packs {
        println!("source pack: {} ({})", pack.pack_id(), pack.content_hash());
    }
    println!("Negative enthalpy is exothermic. Source confidence probabilities are unstated.");
    println!(
        "Reported source uncertainty magnitudes, where available, are not a propagated reaction interval."
    );
    println!(
        "This uses actually reacted pure material; it does not predict conversion, plant temperature, kiln fuel duty or operating reagent doses."
    );
    Ok(())
}
