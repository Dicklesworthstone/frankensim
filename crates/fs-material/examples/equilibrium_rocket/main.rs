//! Run: cargo run -p fs-material --example equilibrium_rocket -- 300 3000000 100000 2
//!
//! Optional arguments: feed temperature [K], chamber pressure [Pa], ambient
//! pressure [Pa], exit/throat area ratio. This is an explicitly manufactured
//! A <-> B isomer system, not a real propellant model or an engine prediction.
//! Original NASA-9 cards carry variable cp and a declared product formation
//! enthalpy; the HP chamber adds no separate heat of reaction. Its equilibrium
//! products then expand with frozen chemistry through the existing nozzle.

use fs_material::thermochemical_gas::ThermochemicalGasModel;
use fs_material::thermochemical_gas::energy::GasEnergySolve;
use fs_material::thermochemical_gas::energy::equilibrium::GasEquilibriumChamberSolve;
use fs_material::thermochemical_gas::energy::parcel::transfer::nozzle::rocket::{
    RocketNozzle, RocketNozzleSolve,
};
use fs_qty::{Amount, Temperature};
use fs_thermochem::equilibrium::composition::{
    EquilibriumEnthalpySolveV1, EquilibriumSolveV1, IdealGasEquilibriumV1,
};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{
    ChargeVector, ElementId, ElementalMatrix, Nasa9RegionV1, Nasa9StandardStateModelV1, SpeciesId,
    UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K,
};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
#[path = "../material_gas_spring/source.rs"]
mod source;

fn main() -> Result<()> {
    let arguments = std::env::args().skip(1).collect::<Vec<_>>();
    if arguments.len() > 4 {
        return Err("expected at most: feed-K chamber-Pa ambient-Pa exit/throat-area".into());
    }
    let mut values: [f64; 4] = [300.0, 3e6, 1e5, 2.0];
    for (value, argument) in values.iter_mut().zip(arguments) {
        *value = argument.parse()?;
        if !value.is_finite() {
            return Err("finite example arguments required".into());
        }
    }
    let [inlet, chamber_pressure, ambient_pressure, area_ratio] = values;
    // Existing fixture grammar, original transport/source owners. Catalog
    // fractions cover every possible species; feed amounts are separate below.
    let catalog = source::load(
        "A 0.5 0.01 4 0.001 1.7e-5 0.025\nB 0.5 0.01 4 0.001 2.0e-5 0.15",
        "equilibrium_rocket: authored variable-cp A/B fixture",
    )?;
    let models = catalog
        .thermodynamics()
        .components()
        .iter()
        .map(|original| -> Result<Nasa9StandardStateModelV1> {
            let mut card = original.regions()[0].card().clone();
            let formation_j_mol = if original.species().as_str() == "B" {
                -24000.0
            } else {
                0.0
            };
            card.parameters
                .get_mut("a7")
                .ok_or("fixture missing NASA-9 a7")?
                .value = formation_j_mol / UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;
            card.provenance.source = format!(
                "SYNTHETIC equilibrium_rocket: {} formation offset {formation_j_mol} J/mol",
                original.species()
            );
            Ok(Nasa9StandardStateModelV1::new(
                original.species().clone(),
                original.molar_mass(),
                original.convention().clone(),
                vec![Nasa9RegionV1::from_card(card)?],
            )?)
        })
        .collect::<Result<Vec<_>>>()?;
    let bank = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            models.clone(),
            catalog.thermodynamics().declared_composition().clone(),
        )?,
        catalog.transport_materials().to_vec(),
        None,
    )?;
    let ids: Vec<_> = models.iter().map(|model| model.species().clone()).collect();
    // Explicit manufactured conserved atom inventory; names imply no formula.
    let chemistry = IdealGasEquilibriumV1::new(
        ElementalMatrix::new(vec![ElementId::new("H")?], ids.clone(), vec![vec![1, 1]])?,
        ChargeVector::new(ids, vec![0, 0])?,
        models,
    )?;
    let chamber = bank.equilibrium_chamber_hp(
        &chemistry,
        &source::point(inlet, chamber_pressure)?,
        &[(SpeciesId::new("A")?, Amount::new(2.0))],
        Temperature::new(inlet),
        GasEquilibriumChamberSolve {
            equilibrium: EquilibriumEnthalpySolveV1 {
                equilibrium: EquilibriumSolveV1::default(),
                minimum_temperature_k: 200.0,
                maximum_temperature_k: 1800.0,
                enthalpy_tolerance_j: 1e-6,
                maximum_temperature_iterations: 128,
            },
            mass_tolerance_kg: 1e-12,
        },
        || false,
    )?;
    let nozzle = RocketNozzle {
        throat_area_m2: 1e-4,
        exit_area_m2: 1e-4 * area_ratio,
    };
    let flow = chamber.parcel().rocket_nozzle_flow(
        ambient_pressure,
        nozzle,
        RocketNozzleSolve {
            nozzle: GasEnergySolve {
                minimum_temperature_k: 200.0,
                maximum_temperature_k: 1800.0,
                energy_tolerance_j_kg: 1e-6,
                maximum_iterations: 128,
            },
            relative_mass_flow_tolerance: 1e-10,
        },
        || false,
    )?;
    println!("Source: manufactured neutral A/B isomers; no real propellant accuracy claim.");
    println!("Chemistry: adiabatic HP chamber; frozen products in the nozzle.");
    println!(
        "chamber_temperature_K,{:.9}",
        chamber.equilibrium().temperature().value()
    );
    for product in chamber.equilibrium().species() {
        println!(
            "chamber_mole_fraction_{},{:.12}",
            product.species, product.mole_fraction
        );
    }
    println!("chamber_mass_kg,{:.12}", chamber.parcel().mass_kg());
    println!("chamber_volume_m3,{:.12}", chamber.parcel().volume_m3());
    println!(
        "chamber_mass_residual_kg,{:.6e}",
        chamber.balance().mass_residual_kg
    );
    println!(
        "chamber_enthalpy_residual_J,{:.6e}",
        chamber.balance().enthalpy_residual_j
    );
    println!(
        "freeze_enthalpy_residual_J,{:.6e}",
        chamber.balance().freeze_enthalpy_residual_j
    );
    println!(
        "maximum_mole_fraction_correction,{:.6e}",
        chamber.balance().maximum_mole_fraction_correction
    );
    println!("nozzle_regime,{:?}", flow.regime);
    println!("nozzle_mass_flow_kg_s,{:.9}", flow.mass_flow_kg_s);
    println!(
        "nozzle_exit_temperature_K,{:.9}",
        flow.exit().state().temperature
    );
    println!(
        "nozzle_exit_pressure_Pa,{:.9}",
        flow.exit().state().pressure
    );
    println!("nozzle_exit_mach,{:.9}", flow.exit_mach_number);
    println!("nozzle_exit_velocity_m_s,{:.9}", flow.exit_velocity_m_s);
    println!("thrust_N,{:.9}", flow.thrust_n);
    println!("specific_impulse_s,{:.9}", flow.specific_impulse_s);
    println!(
        "nozzle_energy_residual_J_kg,{:.6e}",
        flow.energy_residual_j_kg
    );
    println!(
        "nozzle_relative_mass_flow_residual,{:.6e}",
        flow.relative_mass_flow_residual
    );
    Ok(())
}
