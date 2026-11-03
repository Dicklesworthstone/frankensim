//! Run: cargo run -p fs-material --example grinding_to_leaching -- --materials PACK_DIR
//!
//! Compile these source directories with `xtask matdb-pack` into PACK_DIR:
//! copper-ii-oxide-solid-usgs-b2131-thermochemistry and
//! calcite-solid-usgs-b2131-thermochemistry (both under data/matdb/seed-v1).
//! The normalized files must retain those names with the `.fsmatpk` suffix.
//!
//! Only molar masses come from these actual USGS packs. All mill rates, feed
//! sizes/flows, 50/50 mineral fractions, liberated-sphere interpretations,
//! reaction resistances, solid molar densities, and bath conditions are authored
//! illustrations. They are not measured ore kinetics or plant predictions.
//! No elemental assays are supplied or inferred. The formal H2CO3 endpoint
//! does not solve carbonate speciation, and shared H+ does not solve activities.

use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, NormalizedPack, Provenance,
    QueryPoint, SelectionPolicy,
};
use fs_material::{
    comminution::{
        breakage::{
            BreakageColumn, BreakageTestContext, BreakageTransfer, SizeClass,
            SizeRedistributionCard,
        },
        kinetics::{
            CONTACT_TIME_AXIS, RATE_PREFIX, SELECTION_BREAKAGE_LAW, SelectionBreakageModel,
            UniformizationBudget,
        },
        leach_charge::{
            CollectedLeachCharge, MineralClassAllocation, MolarMassSource, ReactiveRadius,
        },
    },
    reaction_batch::{
        ELAPSED_TIME_AXIS, TEMPERATURE_AXIS, VOLUME_AXIS,
        shrinking_core::{
            CONCENTRATION_AXIS, DIFFUSION_RESISTANCE, DiffusionResistance, FILM_RESISTANCE, LAW,
            LOADING_AXIS, RADIUS_AXIS, SOLID_MOLAR_DENSITY, SURFACE_RESISTANCE,
            ShrinkingCoreReaction, TransferResistance,
            population::{LeachPopulation, PopulationNumerics, ShrinkingCoreBatch},
        },
    },
    separation::{ParticleClass, ParticleStream},
};
use fs_qty::{
    Amount, AmountConcentration, Dimensionless, Dims, Frequency, Length, MassFlowRate,
    QuantitySpec, Temperature, Time, Volume,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    SpeciesId,
    stoichiometry::{minerals::MineralReaction, network::ReactionNetwork},
};
use std::{collections::BTreeMap, error::Error, path::Path};

type Result<T> = std::result::Result<T, Box<dyn Error>>;
const COLLECTION_SECONDS: f64 = 0.1;
const ACID_MOL: f64 = 40.0;
const VOLUME_M3: f64 = 0.01;
const SIZES: [(&str, f64, f64); 3] = [
    ("coarse", 300e-6, 20.0),
    ("medium", 150e-6, 10.0),
    ("fine", 75e-6, 5.0),
];
const DECLARATION: &str = "Authored grinding-to-leaching scenario, not experimental data: \
    300/150/75 um classes, 20/10/5 kg/s solids, 7 kg/s water; coarse/medium/fine \
    selection rates 0.02/0.01/0 per second; next-finer daughters; each class \
    explicitly 50 mass% CuO and 50 mass% calcite, liberated spheres. CuO \
    film/surface/diffusion resistances 1000 s/m, 12000 s/m, 1e8 s/m2 and \
    solid molar density 60000 mol/m3; calcite resistances half those values \
    and molar density 27000 mol/m3. Bath: 40 mol H+, 0.01 m3, 298.15 K.";

fn provenance() -> Provenance {
    Provenance {
        source: DECLARATION.into(),
        license: "project-authored example".into(),
        artifact: Some(hash_domain(
            "grinding-to-leaching-example-v1",
            DECLARATION.as_bytes(),
        )),
    }
}
fn card(
    law: &str,
    parameters: impl IntoIterator<Item = (String, f64, Dims)>,
    validity: ValidityDomain,
) -> ConstitutiveModelCard {
    ConstitutiveModelCard {
        law: LawId(law.into()),
        law_version: 1,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        parameters: parameters
            .into_iter()
            .map(|(key, value, dims)| (key, LawParameter { value, dims }))
            .collect(),
        validity,
        sources: vec![provenance().artifact.expect("authored artifact")],
        provenance: provenance(),
    }
}
fn material(directory: &Path, name: &str, phase_axis: &str) -> Result<MolarMassSource> {
    let bytes = std::fs::read(directory.join(format!("{name}.fsmatpk")))?;
    let decoded = NormalizedPack::from_bytes(&bytes)?;
    let pack = NormalizedPack::from_bytes_verified(decoded.content_hash(), &bytes)?;
    if pack.pack_id() != name || pack.to_bytes() != bytes {
        return Err(format!("expected canonical source pack {name}").into());
    }
    let mut query = QueryPoint::new()
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
    for axis in [
        "source-standard-state",
        "source-elemental-reference-usgs-b2131",
        "source-phase-solid",
        phase_axis,
        "source-pure-species",
        "engineering-reference-model",
    ] {
        query = query.with_quantity(axis, QuantitySpec::dimensional(Dims::NONE), 1.0)?;
    }
    query = query.with_quantity(
        "source-specimen-qualified",
        QuantitySpec::dimensional(Dims::NONE),
        0.0,
    )?;
    let source = MolarMassSource::resolve(&pack, &query, SelectionPolicy::SingleClaimOnly)?;
    println!(
        "molar_mass_source={} hash={} kg_per_mol={:.9}; typed usage receipt verified",
        source.pack_id(),
        source.pack_hash(),
        source.molar_mass().value()
    );
    Ok(source)
}
fn grinding() -> Result<(SelectionBreakageModel, ParticleStream)> {
    let daughters = SizeRedistributionCard::new(
        SIZES
            .iter()
            .map(|(id, size, _)| SizeClass {
                id: (*id).into(),
                diameter: Length::new(*size),
            })
            .collect(),
        [("coarse", "medium"), ("medium", "fine"), ("fine", "fine")]
            .into_iter()
            .map(|(input, output)| BreakageColumn {
                input: input.into(),
                transfers: vec![BreakageTransfer {
                    output: output.into(),
                    mass_fraction: Dimensionless::new(1.0),
                }],
            })
            .collect(),
        BreakageTestContext {
            specimen: "authored CuO/calcite mixture".into(),
            method: "synthetic congruent next-finer daughter chain".into(),
            provenance: provenance(),
        },
        ValidityDomain::unconstrained(),
    )?;
    let source = card(
        SELECTION_BREAKAGE_LAW,
        [("coarse", 0.02), ("medium", 0.01), ("fine", 0.0)]
            .into_iter()
            .map(|(id, rate)| (format!("{RATE_PREFIX}{id}"), rate, Frequency::DIMS)),
        ValidityDomain::unconstrained().with_quantity(
            CONTACT_TIME_AXIS,
            QuantitySpec::dimensional(Time::DIMS),
            0.0,
            120.0,
        ),
    );
    let stream = ParticleStream::new(
        SIZES
            .iter()
            .map(|(id, size, flow)| {
                ParticleClass::new(
                    *id,
                    Length::new(*size),
                    MassFlowRate::new(*flow),
                    BTreeMap::new(),
                )
            })
            .collect::<std::result::Result<Vec<_>, _>>()?,
        MassFlowRate::new(7.0),
    )?;
    Ok((
        SelectionBreakageModel::from_card(daughters, source)?,
        stream,
    ))
}
fn local_law(
    kind: MineralReaction,
    solid: &str,
    scale: f64,
    density: f64,
) -> Result<ShrinkingCoreReaction> {
    let mut validity = ValidityDomain::unconstrained();
    for (axis, dims, low, high) in [
        (ELAPSED_TIME_AXIS, Time::DIMS, 0.0, 100.0),
        (TEMPERATURE_AXIS, Temperature::DIMS, 298.15, 298.15),
        (VOLUME_AXIS, Volume::DIMS, VOLUME_M3, VOLUME_M3),
        (RADIUS_AXIS, Length::DIMS, 1e-6, 1e-3),
        (
            CONCENTRATION_AXIS,
            AmountConcentration::DIMS,
            0.0,
            ACID_MOL / VOLUME_M3,
        ),
        (LOADING_AXIS, AmountConcentration::DIMS, 0.0, 5000.0),
    ] {
        validity = validity.with_quantity(axis, QuantitySpec::dimensional(dims), low, high);
    }
    let source = card(
        LAW,
        [
            (FILM_RESISTANCE, 1000.0 * scale, TransferResistance::DIMS),
            (
                SURFACE_RESISTANCE,
                12000.0 * scale,
                TransferResistance::DIMS,
            ),
            (DIFFUSION_RESISTANCE, 1e8 * scale, DiffusionResistance::DIMS),
            (SOLID_MOLAR_DENSITY, density, AmountConcentration::DIMS),
        ]
        .map(|(name, value, dims)| (name.into(), value, dims)),
        validity,
    );
    Ok(ShrinkingCoreReaction::from_card(
        kind.build()?,
        SpeciesId::new(solid)?,
        SpeciesId::new("H+(aq)")?,
        source,
    )?)
}
fn charge(stream: &ParticleStream, sources: &[MolarMassSource; 2]) -> Result<CollectedLeachCharge> {
    let mut allocations = Vec::new();
    for class in stream.classes() {
        for (index, solid) in ["CuO(s)", "CaCO3(s)"].into_iter().enumerate() {
            allocations.push(MineralClassAllocation {
                population_id: format!("{}/{solid}", class.id()),
                class_id: class.id().into(),
                solid: SpeciesId::new(solid)?,
                mass_fraction: Dimensionless::new(0.5),
                fraction_source: provenance(),
                molar_mass: sources[index].clone(),
                radius: ReactiveRadius::LiberatedSphere {
                    provenance: provenance(),
                },
            });
        }
    }
    Ok(CollectedLeachCharge::collect(
        stream,
        Time::new(COLLECTION_SECONDS),
        allocations,
    )?)
}
fn bath(
    charge: &CollectedLeachCharge,
    models: &[ShrinkingCoreReaction; 2],
) -> Result<ShrinkingCoreBatch> {
    // Same canonical population order and compensated sum as the batch owner.
    let mut totals: BTreeMap<SpeciesId, (f64, f64)> = BTreeMap::new();
    let mut populations = Vec::new();
    for p in charge.populations() {
        let (sum, correction) = totals.entry(p.solid().clone()).or_default();
        let amount = p.amount().value();
        let next = *sum + amount;
        *correction += if sum.abs() >= amount.abs() {
            (*sum - next) + amount
        } else {
            (amount - next) + *sum
        };
        *sum = next;
        let model = if p.solid().as_str() == "CuO(s)" {
            &models[0]
        } else {
            &models[1]
        };
        // The batch supplies original mineral-total loading and all owned axes.
        populations.push(LeachPopulation::new(
            p.id(),
            model.clone(),
            p.radius(),
            p.amount(),
            QueryPoint::new(),
        )?);
    }
    let network = ReactionNetwork::new(models.iter().map(|m| m.reaction().clone()).collect())?;
    let inventory = network
        .species()
        .iter()
        .map(|id| {
            let amount = if id.as_str() == "H+(aq)" {
                ACID_MOL
            } else {
                totals
                    .get(id)
                    .map_or(0.0, |(sum, correction)| sum + correction)
            };
            (id.clone(), Amount::new(amount))
        })
        .collect::<Vec<_>>();
    // Zero product stocks define reaction-product increments. Carrier water is
    // retained independently in kg; no solvent-mass/volume relation is inferred.
    Ok(ShrinkingCoreBatch::new(
        populations,
        SpeciesId::new("H+(aq)")?,
        &inventory,
        Volume::new(VOLUME_M3),
        Temperature::new(298.15),
    )?)
}
fn stock(inventory: &[(SpeciesId, Amount)], name: &str) -> f64 {
    inventory
        .iter()
        .find(|(id, _)| id.as_str() == name)
        .expect("complete admitted inventory")
        .1
        .value()
}
fn run(
    seconds: f64,
    mill: &SelectionBreakageModel,
    feed: &ParticleStream,
    sources: &[MolarMassSource; 2],
    models: &[ShrinkingCoreReaction; 2],
) -> Result<()> {
    let pass = mill.prepare(
        Time::new(seconds),
        &QueryPoint::new(),
        UniformizationBudget {
            maximum_l1_tail: Dimensionless::new(1e-12),
            max_poisson_terms: 256,
            max_coefficient_products: 1_000_000,
        },
        || false,
    )?;
    let charge = charge(&pass.transform(feed)?, sources)?;
    println!(
        "\ngrinding_s={seconds} collection_s={COLLECTION_SECONDS} nominal_L1_tail_bound={:.3e}",
        2.0 * pass.poisson_tail_bound.value()
    );
    for (name, before, after) in [
        (
            "dry",
            feed.dry_solids()?.value() * COLLECTION_SECONDS,
            charge.dry_mass().value(),
        ),
        (
            "carrier_water",
            feed.water().value() * COLLECTION_SECONDS,
            charge.water_mass().value(),
        ),
    ] {
        println!(
            "collection_balance={name} expected_kg={before:.12} collected_kg={after:.12} residual_kg={:.3e}",
            after - before
        );
        if (after - before).abs() > 1e-12 * before.max(1.0) {
            return Err("collection mass balance failed".into());
        }
    }
    println!(
        "unallocated_dry_kg={:.12}; no elemental assay reconciliation is asserted",
        charge.unallocated_mass().value()
    );
    let mut bath = bath(&charge, models)?;
    let original_cuo = stock(bath.inventory(), "CuO(s)");
    let original_calcite = stock(bath.inventory(), "CaCO3(s)");
    println!(
        "initial_elapsed_s={} CuO_mol={original_cuo:.9} calcite_mol={original_calcite:.9} Hplus_mol={ACID_MOL}",
        bath.elapsed().value()
    );
    for step in 1..=2 {
        let advance = bath.advance(Time::new(3.0), PopulationNumerics::default(), || false)?;
        let copper = stock(bath.inventory(), "Cu[2+](aq)");
        let calcium = stock(bath.inventory(), "Ca[2+](aq)");
        let acid = stock(bath.inventory(), "H+(aq)");
        let acid_residual = acid + 2.0 * (copper + calcium) - ACID_MOL;
        let copper_residual = stock(bath.inventory(), "CuO(s)") + copper - original_cuo;
        println!(
            "step={step} previous_elapsed_s={} cumulative_elapsed_s={} Cu_formed_mol={copper:.9} Ca_formed_mol={calcium:.9} Hplus_remaining_mol={acid:.9}",
            advance.previous_elapsed.value(),
            advance.elapsed.value()
        );
        println!(
            "acid_balance_residual_mol={acid_residual:.3e} copper_balance_residual_mol={copper_residual:.3e} estimated_time_error_s={:.3e}",
            advance.estimated_time_error.value()
        );
        for p in bath.populations() {
            println!(
                "  population={} solid={} conversion={:.7} remaining_mol={:.9} core_radius_um={:.5}",
                p.id,
                p.solid.as_str(),
                p.conversion.value(),
                p.remaining.value(),
                p.core_radius.value() * 1e6
            );
        }
    }
    Ok(())
}
fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    let directory = match args.as_slice() {
        [flag, directory] if flag == "--materials" => Path::new(directory),
        _ => return Err("usage: grinding_to_leaching --materials PACK_DIR (see example module docs for required packs)".into()),
    };
    println!("{DECLARATION}");
    println!(
        "Collected carrier water remains a separate mass account; network H2O records reaction-formed water from a zero product baseline."
    );
    let sources = [
        material(
            directory,
            "copper-ii-oxide-solid-usgs-b2131-thermochemistry",
            "source-phase-tenorite",
        )?,
        material(
            directory,
            "calcite-solid-usgs-b2131-thermochemistry",
            "source-phase-calcite",
        )?,
    ];
    let models = [
        local_law(
            MineralReaction::CupricOxideAcidLeach,
            "CuO(s)",
            1.0,
            60_000.0,
        )?,
        local_law(
            MineralReaction::CalciteCarbonicAcidNeutralization,
            "CaCO3(s)",
            0.5,
            27_000.0,
        )?,
    ];
    let (mill, feed) = grinding()?;
    for seconds in [30.0, 120.0] {
        run(seconds, &mill, &feed, &sources, &models)?;
    }
    println!(
        "\nEach second advance retains the original radii, loading, depleted acid and accepted core history. Tail/time estimates exclude source uncertainty and are not certified global bounds."
    );
    Ok(())
}
