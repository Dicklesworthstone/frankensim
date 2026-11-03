//! Run: cargo run -p fs-material --example mineral_leach_train
//!
//! One ideal full leach tank, continued through changes in feed and flow. CuO
//! and calcite compete for the same finite H+ inventory; dissolved Cu/Ca and
//! reaction water are retained and discharged through the actual tank flow.
//! The authored local kinetic constants illustrate coupling, not measured ore
//! performance. All particles share the liquid residence in this ideal model.
//! Calcite uses an explicit bicarbonate endpoint: this does not infer pH,
//! speciation, gas transfer, electroneutrality or temperature feedback.

use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance, QueryPoint,
};
use fs_material::reaction_flow::{
    BilinearRateConstant, FlowIntegrationBudget, LOCAL_BILINEAR_LAW, LocalBilinearReaction,
    RATE_PARAMETER, REACTANT_A_CONCENTRATION_AXIS, REACTANT_B_CONCENTRATION_AXIS, ReactiveCstr,
    TEMPERATURE_AXIS, VOLUME_AXIS,
};
use fs_qty::{
    Amount, AmountConcentration, Dimensionless, QuantitySpec, Temperature, Time, Volume,
    VolumetricFlowRate,
};
use fs_thermochem::{SpeciesId, stoichiometry::minerals::MineralReaction};
use std::collections::BTreeMap;

fn model(
    kind: MineralReaction,
    mineral: &str,
    k: f64,
) -> Result<LocalBilinearReaction, Box<dyn std::error::Error>> {
    let declaration =
        format!("synthetic local {mineral}/H+ law: k={k} m3/(mol s); no ore calibration");
    let source = ConstitutiveModelCard {
        law: LawId(LOCAL_BILINEAR_LAW.into()),
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
                20.0,
            )
            .with_quantity(
                REACTANT_B_CONCENTRATION_AXIS,
                QuantitySpec::dimensional(AmountConcentration::DIMS),
                0.0,
                20.0,
            ),
        sources: vec![hash_domain(
            "mineral-leach-cstr-example-v1",
            declaration.as_bytes(),
        )],
        provenance: Provenance {
            source: declaration,
            license: "project example".into(),
            artifact: None,
        },
    };
    Ok(LocalBilinearReaction::from_card(
        kind.build()?,
        SpeciesId::new(mineral)?,
        SpeciesId::new("H+(aq)")?,
        source,
    )?)
}
fn main() -> Result<(), Box<dyn std::error::Error>> {
    let tank = ReactiveCstr::new(
        Volume::new(1.0),
        vec![
            model(MineralReaction::CupricOxideAcidLeach, "CuO(s)", 0.01)?,
            model(
                MineralReaction::CalciteBicarbonateNeutralization,
                "CaCO3(s)",
                0.002,
            )?,
        ],
        SpeciesId::new("H+(aq)")?,
    )?;
    let mut inventory = tank
        .reaction_network()
        .species()
        .iter()
        .map(|s| {
            let n = match s.as_str() {
                "CuO(s)" => 2.0,
                "CaCO3(s)" => 8.0,
                "H2O(l)" => 55_500.0,
                _ => 0.0,
            };
            (s.clone(), Amount::new(n))
        })
        .collect::<Vec<_>>();
    let conditions = tank
        .models()
        .iter()
        .map(|m| (m.reaction().id().clone(), QueryPoint::new()))
        .collect::<Vec<_>>();
    let budget = FlowIntegrationBudget {
        maximum_step: Time::new(1.0),
        absolute_amount_error: Amount::new(1e-8),
        relative_amount_error: Dimensionless::new(1e-7),
        maximum_attempts: 16_384,
    };
    println!("Synthetic CuO/calcite leaching: one 1 m3 full ideal tank at 298.15 K.");
    println!(
        "Local kinetic rates are explicitly authored examples, not measured mineral calibration."
    );
    println!("CuO + 2 H+ -> Cu[2+] + H2O; calcite + H+ -> Ca[2+] + HCO3[-].");
    println!("All network-active species are tracked; no pH or counterion speciation is inferred.");
    for (name, q, seconds, acid) in [
        ("start with no acid", 0.1, 2.0, 0.0),
        ("feed 8 mol/m3 H+", 0.1, 20.0, 8.0),
        ("stop flow; retained tank continues reacting", 0.0, 5.0, 0.0),
        ("restart at lower flow with acid-free feed", 0.05, 10.0, 0.0),
    ] {
        let inlet = tank
            .reaction_network()
            .species()
            .iter()
            .map(|s| {
                let c = match s.as_str() {
                    "CuO(s)" => 2.0,
                    "CaCO3(s)" => 8.0,
                    "H+(aq)" => acid,
                    "H2O(l)" => 55_500.0,
                    _ => 0.0,
                };
                (s.clone(), AmountConcentration::new(c))
            })
            .collect::<Vec<_>>();
        let result = tank.advance(
            &inventory,
            &inlet,
            VolumetricFlowRate::new(q),
            Time::new(seconds),
            Temperature::new(298.15),
            &conditions,
            budget,
        )?;
        println!("\n{name}: Q={q} m3/s, duration={seconds} s");
        for reaction in &result.reactions.reactions {
            println!(
                "  {}: {:.6} mol extent",
                reaction.reaction.as_str(),
                reaction.extent.value()
            );
        }
        println!("  species: inlet mol / outlet mol / retained mol");
        for (i, (s, n)) in result.final_inventory.iter().enumerate() {
            println!(
                "    {}: {:.6} / {:.6} / {:.6}",
                s.as_str(),
                result.inlet_amounts[i].1.value(),
                result.outlet_amounts[i].1.value(),
                n.value()
            );
        }
        let residual = result
            .balance_residuals
            .iter()
            .map(|(_, n)| n.value().abs())
            .fold(0.0_f64, f64::max);
        let extent_error = result
            .estimated_extent_errors
            .iter()
            .map(|(_, n)| n.value())
            .fold(0.0_f64, f64::max);
        println!("  maximum independent species balance residual: {residual:.3e} mol");
        println!(
            "  accepted/rejected numerical steps: {}/{}; largest extent error estimate: {extent_error:.3e} mol",
            result.accepted_steps, result.rejected_steps
        );
        inventory = result.final_inventory;
    }
    println!(
        "\nSecond-order symmetric splitting; error estimates are not certified global bounds."
    );
    Ok(())
}
