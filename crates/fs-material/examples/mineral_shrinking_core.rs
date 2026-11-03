//! Finite-bath size sensitivity, with explicitly synthetic transport coefficients.
//! Run: cargo run -p fs-material --example mineral_shrinking_core
//! No rate below is a measurement or a prediction for a named mine or ore.

use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance, QueryPoint,
};
use fs_material::reaction_batch::{
    ELAPSED_TIME_AXIS, TEMPERATURE_AXIS, VOLUME_AXIS,
    shrinking_core::{
        CONCENTRATION_AXIS, DIFFUSION_RESISTANCE, DiffusionResistance, FILM_RESISTANCE, LAW,
        LOADING_AXIS, LeachNumerics, RADIUS_AXIS, SOLID_MOLAR_DENSITY, SURFACE_RESISTANCE,
        ShrinkingCoreReaction, TransferResistance,
    },
};
use fs_qty::{Amount, AmountConcentration, Length, QuantitySpec, Temperature, Time, Volume};
use fs_thermochem::{SpeciesId, stoichiometry::minerals::MineralReaction};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut validity = ValidityDomain::unconstrained();
    for (axis, dims, lo, hi) in [
        (ELAPSED_TIME_AXIS, Time::DIMS, 0.0, 86_400.0),
        (TEMPERATURE_AXIS, Temperature::DIMS, 298.15, 298.15),
        (VOLUME_AXIS, Volume::DIMS, 0.01, 0.01),
        (RADIUS_AXIS, Length::DIMS, 50e-6, 200e-6),
        (CONCENTRATION_AXIS, AmountConcentration::DIMS, 0.0, 800.0),
        (LOADING_AXIS, AmountConcentration::DIMS, 100.0, 100.0),
    ] {
        validity = validity.with_quantity(axis, QuantitySpec::dimensional(dims), lo, hi);
    }
    let source = ConstitutiveModelCard {
        law: LawId(LAW.into()),
        law_version: 1,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        parameters: [
            (FILM_RESISTANCE, 1e5, TransferResistance::DIMS),
            (SURFACE_RESISTANCE, 1e7, TransferResistance::DIMS),
            (DIFFUSION_RESISTANCE, 1e10, DiffusionResistance::DIMS),
            (SOLID_MOLAR_DENSITY, 20_000.0, AmountConcentration::DIMS),
        ]
        .into_iter()
        .map(|(key, value, dims)| (key.into(), LawParameter { value, dims }))
        .collect(),
        validity,
        sources: vec![hash_domain(
            "mineral-shrinking-core-example",
            b"synthetic coefficients; no fitted ore kinetics",
        )],
        provenance: Provenance {
            source: "authored size-sensitivity scenario; not calibrated plant data".into(),
            license: "project example".into(),
            artifact: None,
        },
    };
    let solid = SpeciesId::new("CuFeS2(s)")?;
    let reagent = SpeciesId::new("Fe[3+](aq)")?;
    let model = ShrinkingCoreReaction::from_card(
        MineralReaction::FerricChalcopyriteLeach.build()?,
        solid.clone(),
        reagent.clone(),
        source,
    )?;
    println!(
        "SYNTHETIC transport coefficients; separate monodisperse batches, no plant calibration"
    );
    println!(
        "Fixed 10 L bath, 1 mol chalcopyrite, 4 mol ferric reagent, 298.15 K; particle sizes are radii"
    );
    println!(
        "radius_um,time_s,conversion,core_radius_um,Cu_formed_mol,FeIII_remaining_mol,time_error_s"
    );
    for radius in [50e-6, 100e-6, 200e-6] {
        let inventory: Vec<_> = model
            .reaction()
            .terms()
            .iter()
            .map(|term| {
                let n = if term.species == solid {
                    1.0
                } else if term.species == reagent {
                    4.0
                } else {
                    0.0
                };
                (term.species.clone(), Amount::new(n))
            })
            .collect();
        for time in [3600.0, 18_000.0, 86_400.0] {
            let out = model.advance(
                &inventory,
                Volume::new(0.01),
                Length::new(radius),
                Time::new(time),
                Temperature::new(298.15),
                &QueryPoint::new(),
                LeachNumerics::default(),
                || false,
            )?;
            let copper = out
                .final_inventory
                .iter()
                .find(|(id, _)| id.as_str() == "Cu[2+](aq)")
                .ok_or("missing conserved copper product")?
                .1
                .value();
            println!(
                "{:.0},{time:.0},{:.9},{:.6},{copper:.9},{:.9},{:.3e}",
                radius * 1e6,
                out.conversion.value(),
                out.core_radius.value() * 1e6,
                out.reactants[1].remaining.value(),
                out.estimated_time_error.value()
            );
        }
    }
    Ok(())
}
