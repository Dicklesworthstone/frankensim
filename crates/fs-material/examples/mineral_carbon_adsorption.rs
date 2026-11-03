//! Finite solution/carbon contacts with explicitly authored equilibrium parameters.
//! Run: cargo run -p fs-material --example mineral_carbon_adsorption -- [DRY_CARBON_KG]
//! These parameters are synthetic. No gold-ore, carbon, pH or elution calibration
//! is claimed, and no contact time or kinetic recovery is supplied.

use fs_blake3::hash_domain;
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ConstitutiveModelCard, InitialStatePolicy, LawId, LawParameter, Provenance, QueryPoint,
};
use fs_material::aqueous_adsorption::{
    AQUEOUS_LANGMUIR_LAW, AdsorbentLoading, AdsorptionInventory, AqueousLangmuir,
    CAPACITY_PARAMETER, CONCENTRATION_AXIS, HALF_SATURATION_PARAMETER, TEMPERATURE_AXIS,
};
use fs_qty::{
    Amount, AmountConcentration, Mass, QuantitySpec, Temperature, Volume,
    chemistry::ElementId,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use std::collections::BTreeMap;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.len() > 1 {
        return Err("usage: mineral_carbon_adsorption [DRY_CARBON_KG]".into());
    }
    let carbon = args
        .first()
        .map(|s| s.parse::<f64>())
        .transpose()?
        .unwrap_or(100.0);
    let source = ConstitutiveModelCard {
        law: LawId(AQUEOUS_LANGMUIR_LAW.into()),
        law_version: 1,
        parameters: BTreeMap::from([
            (CAPACITY_PARAMETER.into(), LawParameter { value: 0.2, dims: AdsorbentLoading::DIMS }),
            (HALF_SATURATION_PARAMETER.into(), LawParameter { value: 0.05, dims: AmountConcentration::DIMS }),
        ]),
        state_schema_version: 0,
        initial_state: InitialStatePolicy::ZeroInternalState,
        validity: ValidityDomain::unconstrained()
            .with_quantity(CONCENTRATION_AXIS, QuantitySpec::dimensional(AmountConcentration::DIMS), 0.0, 0.2)
            .with_quantity(TEMPERATURE_AXIS, QuantitySpec::semantic(SemanticType::new(QuantityKind::AbsoluteTemperature, ValueForm::Static)), 298.15, 298.15),
        sources: vec![hash_domain("authored-carbon-contact-example", b"synthetic qmax0.2 half0.05")],
        provenance: Provenance {
            source: "Authored equilibrium example; no experimental gold/carbon fit. Isotherm form: USGS PHREEQC example 19 equation 35.".into(),
            license: "project example".into(),
            artifact: None,
        },
    };
    let model = AqueousLangmuir::from_card(ElementId::new("Au")?, source)?;
    let first = model.equilibrate(
        AdsorptionInventory {
            dissolved: Amount::new(10.0),
            adsorbed: Amount::new(0.0),
        },
        Volume::new(100.0),
        Mass::new(carbon),
        Temperature::new(298.15),
        &QueryPoint::new(),
    )?;
    println!("Authored finite-bath gold/carbon equilibrium; not a calibrated plant model");
    println!("dry carbon [kg]: {carbon}");
    println!(
        "initial solution metal [mol]: {}",
        first.initial.dissolved.value()
    );
    println!(
        "final solution metal [mol]: {}",
        first.final_inventory.dissolved.value()
    );
    println!(
        "loaded carbon metal [mol]: {}",
        first.final_inventory.adsorbed.value()
    );
    println!(
        "final solution concentration [mol/m3]: {}",
        first.equilibrium_concentration.value()
    );
    println!(
        "vacant carbon capacity [mol]: {}",
        first.vacant_capacity.value()
    );
    println!(
        "metal balance residual [mol]: {}",
        first.conservation_residual.value()
    );
    // Separate the liquid, then contact the same loaded carbon with a clean bath.
    let wash = model.equilibrate(
        AdsorptionInventory {
            dissolved: Amount::new(0.0),
            adsorbed: first.final_inventory.adsorbed,
        },
        Volume::new(100.0),
        Mass::new(carbon),
        Temperature::new(298.15),
        &QueryPoint::new(),
    )?;
    println!(
        "fresh-bath metal after second contact [mol]: {}",
        wash.final_inventory.dissolved.value()
    );
    println!(
        "remaining carbon metal [mol]: {}",
        wash.final_inventory.adsorbed.value()
    );
    println!(
        "second-contact transfer to carbon [mol]: {}",
        wash.transferred_to_adsorbent.value()
    );
    println!(
        "second-contact balance residual [mol]: {}",
        wash.conservation_residual.value()
    );
    Ok(())
}
