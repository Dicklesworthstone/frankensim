//! Finite crystalline hydrated lime and cyanide at a shared, source-backed pH.
//!
//! Run with --materials PACK_DIR; --help lists the required source profiles.
//! The supplied lime mass is a finite calcium component, not an assertion of
//! complete dissolution. The original mineral/acid-base owner jointly solves
//! dissolved Ca, remaining solid, HCN/CN and water. NaCN mass is explicitly
//! declared already dissolved, with its sodium counterion retained.
#![deny(unsafe_code)]

use std::{error::Error, path::Path};

use fs_matdb::{NormalizedPack, PropertyKey, QueryPoint, SelectionPolicy};
use fs_qty::{
    Amount, Dims, Mass, MolarMass, QuantitySpec,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    SpeciesId,
    stoichiometry::{
        minerals::MineralReaction,
        thermodynamics::{
            FormationGibbsSource, StandardReactionGibbs,
            acid_base::{
                AcidDissociation, AcidFamily, IdealAcidBaseSystem, SpectatorIon,
                WaterAutoprotolysis,
                precipitation::{
                    HydrolysisPrecipitation, PrecipitationOptions, PrecipitationSpeciation,
                    PrecipitationStatus,
                },
            },
        },
    },
};

type Result<T> = std::result::Result<T, Box<dyn Error>>;
const REFERENCE: &str = "source-elemental-reference-nbs1989";
const FAMILY: &str = "uncomplexed-cyanide";
const USAGE: &str = "usage: lime_buffered_cyanide --materials PACK_DIR [--lime-mass-kg VALUE] [--cyanide-mass-kg VALUE] [--solvent-mass-kg VALUE]";
const PACKS: [(&str, &str); 9] = [
    ("hydrogen-cyanide-aqueous-nbs1989", "HCN(aq)"),
    ("hydrogen-ion-aqueous-nbs1989", "H+(aq)"),
    ("cyanide-ion-aqueous-nbs1989", "CN[-](aq)"),
    ("hydroxide-ion-aqueous-nbs1989", "OH[-](aq)"),
    ("water-liquid-nbs1989-thermochemistry", "H2O(l)"),
    ("calcium-ion-aqueous-nbs1989", "Ca[2+](aq)"),
    (
        "calcium-hydroxide-solid-nbs1989-thermochemistry",
        "Ca(OH)2(s)",
    ),
    ("sodium-ion-aqueous-nbs1989", "Na[+](aq)"),
    (
        "sodium-cyanide-cubic-solid-nbs1989-thermochemistry",
        "NaCN(s)",
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
        4 => flags.extend([("source-phase-liquid", 1.0), ("source-pure-species", 1.0)]),
        6 | 8 => flags.extend([("source-phase-solid", 1.0), ("source-pure-species", 1.0)]),
        _ => {
            flags.extend([
                ("source-phase-aqueous", 1.0),
                ("source-aqueous-unit-activity-infinite-dilution", 1.0),
            ]);
            if index != 0 {
                flags.extend([
                    ("source-individual-ion", 1.0),
                    ("source-ionic-reference-hydrogen-zero", 1.0),
                ]);
            }
        }
    }
    if index == 8 {
        flags.push(("source-phase-sodium-cyanide-cubic", 1.0));
    }
    for (axis, value) in flags {
        at = at.with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), value)?;
    }
    Ok(at)
}

struct Materials {
    packs: Vec<NormalizedPack>,
    points: Vec<QueryPoint>,
    species: Vec<SpeciesId>,
}
impl Materials {
    fn load(directory: &Path) -> Result<Self> {
        let packs = PACKS
            .iter()
            .map(|(name, _)| {
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
            species: PACKS
                .iter()
                .map(|(_, id)| SpeciesId::new(*id))
                .collect::<std::result::Result<Vec<_>, _>>()?,
        })
    }

    fn sources(&self, indices: &[usize]) -> Vec<FormationGibbsSource<'_>> {
        indices
            .iter()
            .map(|&i| FormationGibbsSource {
                species: &self.species[i],
                pack: &self.packs[i],
                at: &self.points[i],
                reference_axis: REFERENCE,
            })
            .collect()
    }

    fn molar_mass(&self, index: usize) -> Result<MolarMass> {
        let answer = self.packs[index].claims().query_typed(
            &PropertyKey::new("molar-mass", MolarMass::DIMS),
            &self.points[index],
            SelectionPolicy::SingleClaimOnly,
        )?;
        self.packs[index].claims().verify_receipt(&answer.receipt)?;
        let value = answer.evidence.value.value;
        if !value.is_finite() || value <= 0.0 {
            return Err("positive finite source molar mass required".into());
        }
        println!(
            "molar_mass_species={} kg_per_mol={:.12e} receipt={}",
            self.species[index].as_str(),
            value,
            answer.receipt.content_hash(),
        );
        Ok(MolarMass::new(value))
    }

    fn lime_model(&self) -> Result<HydrolysisPrecipitation> {
        let water = WaterAutoprotolysis::from_sources(
            MineralReaction::WaterDissociation.build()?,
            self.species[4].clone(),
            self.species[1].clone(),
            self.species[3].clone(),
            &self.sources(&[4, 1, 3]),
        )?;
        let cyanide = AcidDissociation::from_sources(
            MineralReaction::HydrogenCyanideDissociation.build()?,
            self.species[0].clone(),
            self.species[2].clone(),
            self.species[1].clone(),
            &self.sources(&[0, 1, 2]),
        )?;
        print_standard(water.standard());
        print_standard(cyanide.standard());
        let aqueous =
            IdealAcidBaseSystem::new(water, vec![AcidFamily::new(FAMILY, vec![cyanide])?])?;
        let model = HydrolysisPrecipitation::from_sources(
            aqueous,
            MineralReaction::CalciumHydroxidePrecipitation.build()?,
            self.species[5].clone(),
            self.species[6].clone(),
            &self.sources(&[5, 4, 6, 1]),
        )?;
        print_standard(model.standard());
        Ok(model)
    }
}

fn print_standard(standard: &StandardReactionGibbs) {
    println!(
        "reaction={} delta_G_J_per_mol={:.9} ln_K={:.12} T_K={} p_ref_Pa={}",
        standard.reaction().as_str(),
        standard.molar().value(),
        standard.log_equilibrium_constant().value(),
        standard.temperature().value(),
        standard.reference_pressure().value(),
    );
    for term in standard.terms() {
        println!(
            "  species={} source={} receipt={}",
            term.term.species.as_str(),
            term.pack_id,
            term.answer.receipt.content_hash(),
        );
    }
}

fn case(
    model: &HydrolysisPrecipitation,
    calcium_total: Amount,
    lime_molar_mass: MolarMass,
    cyanide_total: Amount,
    sodium: SpeciesId,
    solvent_mass: Mass,
) -> Result<PrecipitationSpeciation> {
    // The NBS aqueous reference is unit molality. One kg of solvent therefore
    // corresponds to one mol at unit ideal activity; this is not liquid volume.
    let basis = Amount::new(solvent_mass.value());
    let result = model.speciate(
        calcium_total,
        &[(FAMILY.to_owned(), cyanide_total)],
        &[SpectatorIon {
            species: sodium,
            charge: 1,
            amount: cyanide_total,
        }],
        basis,
        PrecipitationOptions::default(),
        || false,
    )?;
    let family = result
        .aqueous
        .families
        .iter()
        .find(|family| family.id == FAMILY)
        .ok_or("missing cyanide analytical family")?;
    let amount = |id| {
        family
            .forms
            .iter()
            .find(|form| form.species.as_str() == id)
            .map(|form| form.amount.value())
            .ok_or("missing cyanide family member")
    };
    let free = amount("CN[-](aq)")?;
    let hcn = amount("HCN(aq)")?;
    println!(
        "case: solvent_kg={:.9} total_lime_kg={:.12e} calcium_component_mol={:.12e} dissolved_NaCN_mol={:.12e}",
        solvent_mass.value(),
        calcium_total.value() * lime_molar_mass.value(),
        calcium_total.value(),
        cyanide_total.value(),
    );
    println!(
        "  phase={:?} dissolved_Ca_mol={:.12e} solid_CaOH2_mol={:.12e} solid_CaOH2_kg={:.12e}",
        result.status,
        result.dissolved_amount.value(),
        result.solid_amount.value(),
        result.solid_amount.value() * lime_molar_mass.value(),
    );
    println!(
        "  pH={:.9} free_CN_mol={:.12e} HCN_mol={:.12e} OH_mol={:.12e}",
        result.aqueous.ph.value(),
        free,
        hcn,
        result.aqueous.hydroxide_amount.value(),
    );
    println!(
        "  residual_Ca_mol={:.3e} residual_CN_mol={:.3e} residual_charge_mol_equiv={:.3e} ln_dissolution_IAP_over_K={:?}",
        result.component_balance_residual.value(),
        free + hcn - cyanide_total.value(),
        result.aqueous.charge_balance_residual.value(),
        result.log_saturation_ratio.map(|value| value.value()),
    );
    println!(
        "  dissolved_bracket_mol=[{:.12e},{:.12e}] pH_bracket=[{:.12},{:.12}] mineral_iterations={} aqueous_solves={} aqueous_iterations={}",
        result.dissolved_bracket[0].value(),
        result.dissolved_bracket[1].value(),
        result.aqueous.ph_bracket[0].value(),
        result.aqueous.ph_bracket[1].value(),
        result.iterations,
        result.aqueous_solves,
        result.aqueous_iterations,
    );
    Ok(result)
}

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.is_empty() || args.as_slice() == ["--help"] {
        println!("{USAGE}\nCompile these profiles as PACK_DIR/PACK_ID.fsmatpk:");
        for (name, _) in PACKS {
            println!("  {name}");
        }
        println!(
            "Without --lime-mass-kg, compare no lime, 0.001 mol and 0.1 mol of pure crystalline Ca(OH)2. Defaults: 1 kg solvent and 0.000490077 kg already-dissolved cubic NaCN reference equivalent."
        );
        return Ok(());
    }
    if args.len() < 2 || args.len() % 2 != 0 || args[0] != "--materials" {
        return Err(USAGE.into());
    }
    let (mut lime_mass, mut cyanide_mass, mut solvent_mass) = (None, None, None);
    for pair in args[2..].chunks_exact(2) {
        let slot = match pair[0].as_str() {
            "--lime-mass-kg" => &mut lime_mass,
            "--cyanide-mass-kg" => &mut cyanide_mass,
            "--solvent-mass-kg" => &mut solvent_mass,
            _ => return Err(USAGE.into()),
        };
        if slot.is_some() {
            return Err(format!("duplicate option {}", pair[0]).into());
        }
        let value: f64 = pair[1].parse()?;
        if !value.is_finite() || value < 0.0 {
            return Err("mass inputs must be finite and nonnegative".into());
        }
        *slot = Some(value);
    }
    let solvent = Mass::new(solvent_mass.unwrap_or(1.0));
    if solvent.value() <= 0.0 {
        return Err("solvent mass must be positive".into());
    }
    let materials = Materials::load(Path::new(&args[1]))?;
    let model = materials.lime_model()?;
    let lime_molar_mass = materials.molar_mass(6)?;
    let salt_molar_mass = materials.molar_mass(8)?;
    let dissolution = MineralReaction::SodiumCyanideDissolution.build()?;
    print_standard(&dissolution.standard_gibbs_change(&materials.sources(&[8, 7, 2]))?);
    let cyanide = dissolution.extent_from_reacted_mass(
        &materials.species[8],
        Mass::new(cyanide_mass.unwrap_or(0.000490077)),
        salt_molar_mass,
    )?;
    let totals = match lime_mass {
        Some(value) => vec![value / lime_molar_mass.value()],
        None => vec![0.0, 0.001, 0.1],
    };
    let mut results = Vec::new();
    for total in totals {
        if !total.is_finite() {
            return Err("lime formula-unit amount is not representable".into());
        }
        results.push(case(
            &model,
            Amount::new(total),
            lime_molar_mass,
            cyanide,
            materials.species[7].clone(),
            solvent,
        )?);
    }
    if lime_mass.is_none() && cyanide_mass.is_none() && solvent_mass.is_none() {
        let expected = [
            PrecipitationStatus::EmptyComponent,
            PrecipitationStatus::SolidAbsent,
            PrecipitationStatus::SolidPresent,
        ];
        if results
            .iter()
            .zip(expected)
            .any(|(r, status)| r.status != status)
            || results
                .windows(2)
                .any(|r| r[0].aqueous.ph.value() >= r[1].aqueous.ph.value())
        {
            return Err(
                "default finite-lime comparison did not exhibit depletion and buffered coexistence"
                    .into(),
            );
        }
    }
    println!(
        "Model boundary: fixed solvent, 298.15 K, ideal activities, one crystalline hydrated-lime phase. NaCN is already dissolved. No carbonate, sulfate, calcium complexes, other precipitates, dissolution kinetics, gas loss, commercial purity, heat feedback or plant dose calibration is inferred."
    );
    Ok(())
}
