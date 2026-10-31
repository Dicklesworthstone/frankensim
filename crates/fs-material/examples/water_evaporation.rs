//! Finite-atmosphere evaporation/condensation with original source calorics.
//!
//! Create a new PACK_DIRECTORY first; the pack compiler refuses overwrites.
//! Compile these existing data/matdb/seed-v1/NAME/manifest.tsv sources with
//! `cargo run -p xtask -- matdb-pack --manifest ... --out PACK_DIRECTORY/NAME.pack`:
//! nitrogen, nitrogen-nasa9-mcbride2002, nitrogen-transport-svehla1995,
//! water-vapor, water-vapor-nasa9-mcbride2002,
//! water-vapor-transport-iapws2008-2011.
//!
//! Run `cargo run -p fs-material --example water_evaporation -- PACK_DIRECTORY
//! [evaporate|condense] [duration_s=30] [max_dt_s=0.25] [htc_w_m2_k=12]
//! [mass_transfer_m_s=0.005]`. Refine max_dt_s to resolve the first-order
//! heat/phase split; use htc=0 to inspect the contribution of sensible heat.
//!
//! A 20 g water layer at 300 K exchanges with 10 L of N2/H2O at 310 K and 90 kPa.
//! Original species, NASA-9 and transport packs are decoded by the shared gas
//! loader. IF97 supplies the liquid chart and saturation. Film coefficients
//! and reference geometry are explicit inputs, not inferred from source names.
//! Initial vapor mole fraction is 0.01 (evaporate) or 0.06 (condense).
//! These transport packs admit 75..101.325 kPa; source excursions refuse.
//!
//! Liquid H + gas U is conserved in fixed reference geometry. Liquid pressure
//! remains 90 kPa; this is not a sealed moving two-phase vessel. There is no
//! external heater, dryout topology, freezing or boiling. The explicit
//! IF97/NASA compatibility budget is 1000 J/kg; neither original law is altered
//! to erase its temperature-dependent difference from the other source.

use fs_material::evaporation::{
    WaterEvaporation,
    gas::{
        GasEvaporationSolve, GasFilmHeatSolve, GasLiquidBinding,
        timeline::{GasEvaporationTimelineSurface, advance_deposition_with_gas_and_heat},
    },
};
use fs_material::gas_species::input::{
    ThermochemicalPhase,
    packs::{MAX_SOURCE_PACK_BYTES, SourcePackBytes, read_source_parcel},
};
use fs_material::phase::water::{IsobaricWater, WaterChartControls};
use fs_material::thermochemical_gas::energy::{
    GasEnergySolve,
    parcel::{ClosedGasParcel, transfer::species::boundary::GasSpeciesBoundarySolve},
};
use fs_qty::Pressure;
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionSimulation, Event,
    MaterialId, StepConfig, TimedEvent,
};
use fs_thermochem::{SpeciesId, water::if97::bounds::If97HeatCapacityBoundSolve};
use std::{error::Error, io::Read, path::Path};

type Result<T> = std::result::Result<T, Box<dyn Error>>;
const MASS_KG: f64 = 0.02;
const AREA_M2: f64 = 0.01;
const GAS_VOLUME_M3: f64 = 0.01;
const PRESSURE_PA: f64 = 90_000.0;
const ENERGY_TOLERANCE_J: f64 = 1e-6;
const PACKS: [&str; 6] = [
    "nitrogen",
    "nitrogen-nasa9-mcbride2002",
    "nitrogen-transport-svehla1995",
    "water-vapor",
    "water-vapor-nasa9-mcbride2002",
    "water-vapor-transport-iapws2008-2011",
];

fn load_gas(directory: &Path, vapor_fraction: f64) -> Result<ClosedGasParcel> {
    let mut bytes = Vec::with_capacity(PACKS.len());
    for name in PACKS {
        let mut source = Vec::new();
        std::fs::File::open(directory.join(format!("{name}.pack")))?
            .take((MAX_SOURCE_PACK_BYTES + 1) as u64)
            .read_to_end(&mut source)?;
        if source.len() > MAX_SOURCE_PACK_BYTES {
            return Err(format!("{name} exceeds the original gas-loader pack budget").into());
        }
        bytes.push(source);
    }
    let sources = [
        SourcePackBytes {
            species: &bytes[0],
            nasa9: &bytes[1],
            transport: &bytes[2],
        },
        SourcePackBytes {
            species: &bytes[3],
            nasa9: &bytes[4],
            transport: &bytes[5],
        },
    ];
    let declaration = format!(
        "thermochemical 310 {PRESSURE_PA} mole 2 N2 {} H2O {vapor_fraction}",
        1.0 - vapor_fraction
    );
    Ok(ClosedGasParcel::from_isochoric(read_source_parcel(
        &declaration,
        &sources,
        Some("H2O"),
        ThermochemicalPhase::FrozenGas,
        GAS_VOLUME_M3,
    )?))
}

fn main() -> Result<()> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.is_empty() || args.len() > 6 {
        return Err("usage: water_evaporation PACK_DIRECTORY [evaporate|condense] [duration_s=30] [max_dt_s=0.25] [htc_w_m2_k=12] [mass_transfer_m_s=0.005]".into());
    }
    let mode = args.get(1).map_or("evaporate", String::as_str);
    let vapor_fraction = match mode {
        "evaporate" => 0.01,
        "condense" => 0.06,
        _ => return Err("mode must be evaporate or condense".into()),
    };
    let number = |index: usize, fallback: f64| -> Result<f64> {
        Ok(args.get(index).map_or(Ok(fallback), |s| s.parse())?)
    };
    let duration_s = number(2, 30.0)?;
    let max_dt_s = number(3, 0.25)?;
    let htc = number(4, 12.0)?;
    let mass_transfer = number(5, 0.005)?;
    if !duration_s.is_finite()
        || duration_s <= 0.0
        || duration_s > 600.0
        || !max_dt_s.is_finite()
        || max_dt_s <= 0.0
        || !htc.is_finite()
        || htc < 0.0
        || !mass_transfer.is_finite()
        || mass_transfer < 0.0
    {
        return Err("require 0 < duration <= 600 s, positive finite dt, and nonnegative finite film coefficients".into());
    }
    let mut gas = load_gas(Path::new(&args[0]), vapor_fraction)?;
    let catalog = gas.model().clone();
    let water = WaterEvaporation::try_new(
        IsobaricWater::new(
            Pressure::new(PRESSURE_PA),
            WaterChartControls {
                minimum_temperature_k: 273.16,
                maximum_temperature_k: 450.0,
                enthalpy_tolerance_j_kg: 1e-7,
                maximum_iterations: 128,
                maximum_entropy_seam_j_kg_k: 0.1,
                heat_capacity: If97HeatCapacityBoundSolve {
                    maximum_intervals: 65535,
                },
            },
            || false,
        )?,
        [280.0, 330.0],
    )?;
    let liquid = water.liquid_state(300.0)?;
    let vapor = SpeciesId::new("H2O")?;
    let binding = GasLiquidBinding {
        species: vapor.clone(),
        liquid_enthalpy_offset_j_kg: catalog.species_specific_enthalpy_j_kg(&vapor, 300.0)?
            - water.vapor_specific_enthalpy_j_kg(300.0)?,
    };
    let initial_total_mass = MASS_KG + gas.mass_kg();
    let thermal = DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(1),
                segment_id: 1,
                body_id: 1,
                material_id: MaterialId(1),
                mass_kg: MASS_KG,
                specific_h_j_kg: liquid.specific_enthalpy_j_kg(),
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [0.1, 0.0, 0.0],
                    volume_m3: MASS_KG * liquid.specific_volume_m3_kg(),
                    surface_area_m2: AREA_M2,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        }],
        &[(MaterialId(1), &water)],
        DepositionBudget {
            max_cells: 1,
            max_links: 0,
            max_events: 1,
            max_substeps: 200_000,
        },
    )?;
    let thermal_solve = StepConfig {
        max_dt_s,
        energy_tolerance_j: ENERGY_TOLERANCE_J,
    };
    let (mut thermal, _) = thermal.advance_to(0.0, thermal_solve, || Ok(()))?;
    let phase_solve = GasEvaporationSolve {
        gas: GasSpeciesBoundarySolve {
            gas: GasEnergySolve {
                minimum_temperature_k: 280.0,
                maximum_temperature_k: 330.0,
                energy_tolerance_j_kg: 1e-7,
                maximum_iterations: 128,
            },
            mass_tolerance_kg: 1e-12,
            energy_tolerance_j: ENERGY_TOLERANCE_J,
        },
        maximum_inventory_fraction: 0.05,
        maximum_temperature_increment_k: 2.0,
        latent_compatibility_tolerance_j_kg: 1000.0,
    };
    let heat_solve = GasFilmHeatSolve {
        pure_liquid_enthalpy_range_j_kg: water.enthalpy_range_j_kg(),
        temperature_tolerance_k: 1e-8,
        maximum_iterations: 128,
    };
    let surface = GasEvaporationTimelineSurface {
        cell_id: CellId(1),
        law: &water,
        binding: &binding,
        pure_liquid_enthalpy_range_j_kg: water.enthalpy_range_j_kg(),
        mass_transfer_coefficient_m_s: mass_transfer,
    };
    println!(
        "# mode={mode}; liquid_mass={MASS_KG} kg; area={AREA_M2} m2; gas_volume={GAS_VOLUME_M3} m3; initial_pressure={PRESSURE_PA} Pa"
    );
    println!(
        "# htc={htc} W/(m2 K); mass_transfer={mass_transfer} m/s; max_dt={max_dt_s} s; IF97/NASA compatibility=1000 J/kg"
    );
    println!(
        "time_s,liquid_mass_kg,liquid_temperature_k,gas_temperature_k,gas_pressure_pa,vapor_partial_pressure_pa,interface_saturation_ratio,gas_sound_speed_m_s,gas_impedance_pa_s_m,gas_viscosity_pa_s,gas_conductivity_w_m_k,transferred_water_kg,sensible_heat_into_liquid_j,energy_defect_j,mass_defect_kg,rejected_intervals"
    );
    let (mut transferred, mut sensible_heat, mut energy_defect, mut rejected) =
        (0.0, 0.0, 0.0, 0usize);
    for index in 0..=20 {
        let time_s = duration_s * index as f64 / 20.0;
        if index > 0 {
            let candidate = advance_deposition_with_gas_and_heat(
                &thermal,
                surface,
                &gas,
                &catalog,
                time_s,
                thermal_solve,
                phase_solve,
                htc,
                heat_solve,
                || Ok(()),
            )?;
            let (next_thermal, next_gas, _, audit) = candidate.into_parts();
            transferred += audit.transferred_mass_kg;
            sensible_heat += audit.sensible_heat_into_liquid_j;
            energy_defect += audit.energy_residual_j;
            rejected += audit.rejected_intervals;
            thermal = next_thermal;
            gas = next_gas;
        }
        let cell = &thermal.cells()[&CellId(1)];
        let temperature = water.temperature_at_specific_enthalpy_k(cell.specific_h_j_kg)?;
        let vapor_pressure = gas.species_partial_pressure_pa(&vapor)?;
        let saturation_ratio = vapor_pressure / water.saturation(temperature)?.pressure().value();
        let state = gas.state().state();
        let mass_defect = cell.mass_kg + gas.mass_kg() - initial_total_mass;
        println!(
            "{time_s:.8},{:.12e},{temperature:.9},{:.9},{:.9},{vapor_pressure:.9},{saturation_ratio:.9},{:.9},{:.9},{:.12e},{:.12e},{transferred:.12e},{sensible_heat:.12e},{energy_defect:.12e},{mass_defect:.12e},{rejected}",
            cell.mass_kg,
            state.temperature,
            state.pressure,
            state.sound_speed,
            state.characteristic_impedance,
            state.dynamic_viscosity,
            state.thermal_conductivity
        );
    }
    Ok(())
}
