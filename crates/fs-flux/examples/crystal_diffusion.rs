//! Source-bound crystal dissolution/growth coupled to molecular diffusion.
//!
//! Run: cargo run -p fs-flux --example crystal_diffusion
//! Add `-- --water` for coupled retaining-water membranes, and `--nonideal`
//! for explicitly bound virial activities. Both flags may be selected together.
//! The two finite baths exchange dissolved molecules through the original
//! stationary-solvent Maxwell--Stefan face. Crystals remain local; their actual
//! surface areas evolve with solid stock while reference geometry stays fixed.
//! Saturation and film kinetics come from separately selected immutable model
//! cards. All numbers are authored demonstration data, not measured properties
//! of a named chemical or a calibrated experiment.

use std::sync::Arc;

use fs_blake3::hash_bytes;
use fs_evidence::ValidityDomain;
use fs_flux::osmosis::{
    OsmoticCompartment, OsmoticMembrane, OsmoticNetwork, OsmoticPressureVolume,
    OsmoticStepControls,
    coupled::CoupledSolutionControls,
    diffusion::{SolutionDiffusionControls, SolutionDiffusionLink},
};
use fs_matdb::{
    ClaimSet, ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, MaterialCard,
    MaterialStateId, Provenance,
};
use fs_material::osmosis::{
    HydraulicPermeance, MembraneTransport, ReflectionCoefficient, SolutePermeance,
    diffusion::{SolutionDiffusionSolve, SolutionMaxwellStefan},
    nonideal::VirialSolutionModel,
    solubility::{
        kinetics::KineticSolubilitySolution,
        source::{
            DISSOLUTION_ENTHALPY_PARAMETER, NEUTRAL_SOLUBILITY_LAW,
            NeutralSolubilityKineticsSource, NeutralSolubilitySource, PARTICLE_FILM_KINETICS_LAW,
            REFERENCE_SOLID_PARAMETER, REFERENCE_SOLUBILITY_PARAMETER,
            REFERENCE_TEMPERATURE_PARAMETER, TEMPERATURE_AXIS, TRANSFER_CONDUCTANCE_PARAMETER,
        },
    },
};
use fs_qty::{
    Amount, AmountConcentration, Dims, MolarMass, QuantitySpec, Temperature, Volume,
    VolumetricFlowRate,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermochem::{
    MolarEnergyQuantityV1, SpeciesId, solution::solubility::kinetics::SolubilityKineticStateV1,
    transport::maxwell_stefan::MaxwellStefanMolarSolveV1,
};

const T: f64 = 300.0;
const VOLUME: f64 = 1e-6;
const SOLVENT_MOLAR_VOLUME: f64 = 18e-6;
const SOLVENT: &str = "demonstration-solvent";
const SOLUTE: &str = "demonstration-neutral-solute";

fn authored_source() -> Result<NeutralSolubilityKineticsSource, Box<dyn std::error::Error>> {
    let validity = ValidityDomain::unconstrained()
        .with_quantity(
            TEMPERATURE_AXIS,
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            280.0,
            330.0,
        )
        .with_quantity(
            format!("source-solvent:{SOLVENT}"),
            QuantitySpec::dimensional(Dims::NONE),
            1.0,
            1.0,
        )
        .with_quantity(
            format!("source-solute:{SOLUTE}"),
            QuantitySpec::dimensional(Dims::NONE),
            1.0,
            1.0,
        );
    let make_model =
        |law: &str, stateful: bool, parameters: &[(&str, f64, Dims)]| ConstitutiveModelCard {
            law: LawId(law.into()),
            law_version: 1,
            parameters: parameters
                .iter()
                .map(|&(name, value, dims)| (name.into(), LawParameter { value, dims }))
                .collect(),
            state_schema_version: u32::from(stateful),
            initial_state: if stateful {
                InitialStatePolicy::RequiresDeclaredState
            } else {
                InitialStatePolicy::ZeroInternalState
            },
            validity: validity.clone(),
            sources: vec![hash_bytes(
                b"authored crystal-diffusion demonstration, not measured data",
            )],
            provenance: Provenance {
                source: "authored neutral pair; prescribed film geometry and particle reference"
                    .into(),
                license: "demonstration-only".into(),
                artifact: None,
            },
        };
    let equilibrium = make_model(
        NEUTRAL_SOLUBILITY_LAW,
        false,
        &[
            (REFERENCE_TEMPERATURE_PARAMETER, T, Temperature::DIMS),
            (
                REFERENCE_SOLUBILITY_PARAMETER,
                100.0,
                AmountConcentration::DIMS,
            ),
            (
                DISSOLUTION_ENTHALPY_PARAMETER,
                6000.0,
                MolarEnergyQuantityV1::DIMS,
            ),
        ],
    );
    let film = make_model(
        PARTICLE_FILM_KINETICS_LAW,
        true,
        &[
            (
                TRANSFER_CONDUCTANCE_PARAMETER,
                1e-8,
                VolumetricFlowRate::DIMS,
            ),
            (REFERENCE_SOLID_PARAMETER, 2e-4, Amount::DIMS),
        ],
    );
    let (equilibrium_id, film_id) = (equilibrium.content_hash(), film.content_hash());
    let material = Arc::new(MaterialCard::assemble(
        MaterialStateId {
            chemistry: format!("{SOLVENT}/{SOLUTE}"),
            phase: "ideal neutral solution and pure crystal".into(),
            process: "authored fixed-count particles in prescribed film".into(),
            revision: 0,
        },
        ClaimSet::new(),
        vec![equilibrium, film],
    )?);
    let equilibrium = NeutralSolubilitySource::from_material(
        material.clone(),
        equilibrium_id,
        SpeciesId::new(SOLVENT)?,
        SpeciesId::new(SOLUTE)?,
    )?;
    Ok(equilibrium.with_kinetics(material, film_id)?)
}

fn bath(
    source: &NeutralSolubilityKineticsSource,
    dissolved: f64,
    solid: f64,
    thermodynamics: Option<&Arc<VirialSolutionModel>>,
) -> Result<OsmoticCompartment, Box<dyn std::error::Error>> {
    let phases = KineticSolubilitySolution::new(
        Temperature::new(T),
        SOLVENT_MOLAR_VOLUME,
        Volume::new(VOLUME),
        vec![(
            source.kinetics().clone(),
            SolubilityKineticStateV1::from_phases(Amount::new(dissolved), Amount::new(solid))?,
        )],
    )?;
    let phases = match thermodynamics {
        Some(model) => phases.with_nonideal_thermodynamics(Arc::clone(model))?,
        None => phases,
    };
    Ok(OsmoticCompartment::with_kinetic_solubility(
        phases,
        OsmoticPressureVolume::new(VOLUME, 101_325.0, 0.0)?,
        0.1 * VOLUME,
    )?)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let arguments: Vec<_> = std::env::args().skip(1).collect();
    let water = arguments.iter().any(|arg| arg == "--water");
    let nonideal = arguments.iter().any(|arg| arg == "--nonideal");
    // Retain the source binding throughout the simulation, so both original
    // cards and their exact model identities remain available for inspection.
    let source = authored_source()?;
    let thermodynamics = Arc::new(VirialSolutionModel::new(
        T,
        SOLVENT_MOLAR_VOLUME,
        1.0,
        vec![300.0],
        vec![if nonideal { 0.002 } else { 0.0 }],
    )?);
    let diffusion = SolutionMaxwellStefan::new(
        thermodynamics.clone(),
        (SpeciesId::new(SOLVENT)?, MolarMass::new(0.018)),
        vec![(SpeciesId::new(SOLUTE)?, MolarMass::new(0.06))],
        vec![0.0, 1e-9, 1e-9, 0.0],
    )?;
    let binding = nonideal.then_some(&thermodynamics);
    let cells = vec![
        bath(&source, 0.0, 2e-4, binding)?,
        bath(&source, 1.75e-4, 1e-5, binding)?,
    ];
    let links = vec![SolutionDiffusionLink::new(0, 1, 5e-3, 1e-3, diffusion)?];
    let physical = SolutionDiffusionControls {
        constitutive: SolutionDiffusionSolve {
            friction: MaxwellStefanMolarSolveV1::default(),
            entropy_tolerance_w_m3_k: 1e-9,
        },
        maximum_species_fraction: 0.2,
        amount_balance_tolerance_mol: 1e-15,
        free_energy_tolerance_j: 1e-11,
    };
    let mut network = if water {
        let membrane = MembraneTransport::new(
            HydraulicPermeance::new(1e-11)?,
            vec![(ReflectionCoefficient::new(1.0)?, SolutePermeance::new(0.0)?)],
        )?;
        OsmoticNetwork::with_kinetic_coupled_solution_transport(
            vec![SOLUTE.into()],
            cells,
            vec![OsmoticMembrane::new(0, 1, 1e-4, membrane)?],
            links,
            thermodynamics,
            CoupledSolutionControls {
                diffusion: physical,
                water_balance_tolerance_m3: 1e-18,
            },
            1 << 22,
        )?
    } else {
        OsmoticNetwork::with_kinetic_solution_diffusion(
            vec![SOLUTE.into()],
            cells,
            links,
            thermodynamics,
            physical,
            1 << 22,
        )?
    };
    let numerical = OsmoticStepControls {
        maximum_step_s: 10.0,
        relative_tolerance: 1e-5,
        absolute_amount_tolerance_mol: 1e-13,
        absolute_volume_tolerance_m3: 1e-18,
        maximum_steps: 20_000,
        maximum_rejections: 128,
        ..OsmoticStepControls::default()
    };
    let mut cumulative_potential = 0.0;
    println!(
        "time_s,bath,solvent_volume_m3,dissolved_mol,solid_mol,concentration_mol_m3,particle_kA_m3_s,cumulative_potential_change_J,total_component_mol"
    );
    for time in [0.0, 10.0, 30.0, 100.0, 300.0, 1000.0] {
        if time > network.time_s() {
            let report = network.advance(time - network.time_s(), numerical, || false)?;
            cumulative_potential += report
                .coupled_potential_change
                .map(|value| value.total_potential_change_j)
                .or(report.diffusion_free_energy_change_j)
                .ok_or("missing chemical/mechanical potential report")?;
        }
        let phases: Vec<_> = network
            .compartments()
            .iter()
            .map(|cell| {
                cell.kinetic_solubility_snapshot()
                    .ok_or("missing crystal inventory")
            })
            .collect::<Result<_, _>>()?;
        let total: f64 = phases
            .iter()
            .map(|cell| cell.last_steps()[0].state().total_solute().value())
            .sum();
        for (index, cell) in phases.iter().enumerate() {
            let state = cell.last_steps()[0].state();
            println!(
                "{time:.3},{index},{:.12e},{:.12e},{:.12e},{:.9},{:.12e},{cumulative_potential:.12e},{total:.12e}",
                cell.solution().solvent_volume_m3(),
                state.dissolved().value(),
                state.solid().value(),
                state.dissolved().value() / cell.solution().solvent_volume_m3(),
                cell.models()[0].conductance_at(state.solid())?.value(),
            );
        }
    }
    Ok(())
}
