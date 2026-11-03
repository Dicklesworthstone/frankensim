//! Run: cargo run -p fs-material --example competing_mineral_leaching
//!
//! Synthetic CuO/calcite competition for a finite acid inventory. Both
//! empirical rates see the same instantaneous H+ stock. Constants are authored
//! examples, not measured mineral rates, surface mechanisms or ore forecasts.
//! Calcite uses the explicit formal H2CO3 endpoint, not retained bicarbonate;
//! this bookkeeping choice does not solve carbon speciation or gas transfer.

use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance, QueryPoint,
};
use fs_material::reaction_batch::{
    BILINEAR_BATCH_LAW, BilinearBatchReaction, BilinearRateConstant, ELAPSED_TIME_AXIS,
    RATE_PARAMETER, REACTANT_A_CONCENTRATION_AXIS, REACTANT_B_CONCENTRATION_AXIS, TEMPERATURE_AXIS,
    VOLUME_AXIS,
    network::{BatchIntegrationBudget, BilinearBatchNetwork},
};
use fs_qty::{Amount, AmountConcentration, Dimensionless, QuantitySpec, Temperature, Time, Volume};
use fs_thermochem::{SpeciesId, stoichiometry::minerals::MineralReaction};
use std::collections::BTreeMap;

fn model(
    kind: MineralReaction,
    mineral: &str,
    k: f64,
) -> Result<BilinearBatchReaction, Box<dyn std::error::Error>> {
    let declaration =
        format!("authored {mineral} bilinear k={k} m3/(mol s); no mineral calibration");
    let card = ConstitutiveModelCard {
        law: LawId(BILINEAR_BATCH_LAW.into()),
        law_version: 1,
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        parameters: BTreeMap::from([(
            RATE_PARAMETER.into(),
            LawParameter {
                value: k,
                dims: BilinearRateConstant::DIMS,
            },
        )]),
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                ELAPSED_TIME_AXIS,
                QuantitySpec::dimensional(Time::DIMS),
                0.0,
                100.0,
            )
            .with_quantity(
                VOLUME_AXIS,
                QuantitySpec::dimensional(Volume::DIMS),
                1.0,
                1.0,
            )
            .with_quantity(
                TEMPERATURE_AXIS,
                QuantitySpec::dimensional(Temperature::DIMS),
                298.15,
                298.15,
            )
            .with_quantity(
                REACTANT_A_CONCENTRATION_AXIS,
                QuantitySpec::dimensional(AmountConcentration::DIMS),
                0.0,
                10.0,
            )
            .with_quantity(
                REACTANT_B_CONCENTRATION_AXIS,
                QuantitySpec::dimensional(AmountConcentration::DIMS),
                0.0,
                12.0,
            ),
        sources: vec![hash_domain(
            "competing-mineral-leaching-example-v1",
            declaration.as_bytes(),
        )],
        provenance: Provenance {
            source: declaration,
            license: "project example".into(),
            artifact: None,
        },
    };
    Ok(BilinearBatchReaction::from_card(
        kind.build()?,
        SpeciesId::new(mineral)?,
        SpeciesId::new("H+(aq)")?,
        card,
    )?)
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let batch = BilinearBatchNetwork::new(
        vec![
            model(MineralReaction::CupricOxideAcidLeach, "CuO(s)", 0.01)?,
            model(
                MineralReaction::CalciteCarbonicAcidNeutralization,
                "CaCO3(s)",
                0.002,
            )?,
        ],
        SpeciesId::new("H+(aq)")?,
    )?;
    let initial = batch
        .reaction_network()
        .species()
        .iter()
        .map(|id| {
            let mol = match id.as_str() {
                "CuO(s)" | "CaCO3(s)" => 10.0,
                "H+(aq)" => 12.0,
                _ => 0.0,
            };
            (id.clone(), Amount::new(mol))
        })
        .collect::<Vec<_>>();
    let conditions = batch
        .models()
        .iter()
        .map(|model| (model.reaction().id().clone(), QueryPoint::new()))
        .collect::<Vec<_>>();
    let result = batch.advance(
        &initial,
        Volume::new(1.0),
        Time::new(30.0),
        Temperature::new(298.15),
        &conditions,
        BatchIntegrationBudget {
            maximum_step: Time::new(1.0),
            absolute_extent_error: Amount::new(1e-8),
            relative_extent_error: Dimensionless::new(1e-7),
            maximum_attempts: 4096,
        },
    )?;
    println!(
        "Synthetic simultaneous mineral/acid batch; source admission does not establish calibration."
    );
    println!(
        "Inventory/volume is empirical availability, not pure-solid activity or exposed surface area."
    );
    println!("Initial: 10 mol CuO, 10 mol calcite, 12 mol H+; 1 m3; 298.15 K; 30 s.");
    println!(
        "Calcite endpoint: CaCO3 + 2 H+ -> Ca2+ + formal H2CO3; no carbon speciation or gas transfer is inferred."
    );
    for (reaction, (_, error)) in result
        .stage
        .reactions
        .iter()
        .zip(&result.estimated_extent_errors)
    {
        println!(
            "{}: {:.9} mol extent; accumulated local error estimate {:.3e} mol",
            reaction.reaction.as_str(),
            reaction.extent.value(),
            error.value()
        );
    }
    println!("Complete final inventory:");
    for (species, amount) in &result.stage.final_inventory {
        println!("  {}: {:.9} mol", species.as_str(), amount.value());
    }
    println!(
        "Accepted steps: {}; rejected trials: {}",
        result.accepted_steps, result.rejected_steps
    );
    println!(
        "Error estimates are not certified global bounds; equations retain exact atom/charge conservation."
    );
    for model in result.model().models() {
        println!(
            "Original kinetic source: {}",
            model.source().provenance.source
        );
    }
    Ok(())
}
