//! Oxygen-limited methane burning coupled to the original pressure/vent network.
//!
//! cargo run -p fs-airflow --features compartment-fire --example ventilated_compartment_fire --
//!   [end-s=60] [maximum-dt-s=0.02]
//!
//! This is a G0/G1 balance and ventilation demonstration, not an experimental
//! fire prediction. It uses the well-mixed, one-active-layer limit of TwoZone,
//! adiabatic rigid walls, a continuous explicitly supplied methane feed, and
//! instantaneous inventory-limited complete oxidation after transport.
//! There is no prescribed combustion heat. Original formation energies release
//! sensible heat exactly once. No ignition, flame structure, CO/soot production,
//! dissociation, radiation, or building-failure prediction is made.

#[path = "ventilated_compartment_fire/nasa9.rs"]
mod nasa9;

use fs_airflow::compartment::{Ambient, Endpoint, NetworkControl, Room, WallVent, ZoneNetwork};
use fs_airflow::stratified::{DischargeCoefficient, WallOpening};
use fs_airflow::{SourceProvenance, ToleranceBasis};
use fs_gas::multicomponent::MixtureGas;
use fs_gas::multicomponent::compartment::burning::InventoryLimitedReaction;
use fs_gas::multicomponent::compartment::{Layer, TwoZone};
use fs_qty::{
    Acceleration, Amount, Area, Dimensionless, Length, MassFlowRate, Pressure, Temperature, Time,
};
use fs_thermochem::ReactionId;
use std::error::Error;

const FUEL_RATE_KG_S: f64 = 0.002;
const AMBIENT_K: f64 = 300.0;
const ROOM_AREA_M2: f64 = 1.0;
const ROOM_HEIGHT_M: f64 = 1.0;

#[derive(Default)]
struct AcceptedChemistry {
    extent_mol: f64,
    released_enthalpy_j: f64,
    released_power_w: f64,
    chemical_species_change_kg: [f64; 5],
    chemical_physical_source_j: f64,
    chemical_physical_energy_defect_j: f64,
    oxygen_limited: bool,
}

struct RunReport {
    name: &'static str,
    duration_s: f64,
    steps: usize,
    halvings: usize,
    temperature_k: f64,
    pressure_pa: f64,
    final_species_kg: [f64; 5],
    boundary_species_kg: [f64; 5],
    species_balance_kg: [f64; 5],
    element_balance_mol: [f64; 4],
    injected_fuel_kg: f64,
    reacted_fuel_kg: f64,
    initial_oxygen_limited_fuel_kg: f64,
    first_oxygen_limited_s: Option<f64>,
    last_released_power_w: f64,
    released_enthalpy_j: f64,
    physical_energy_balance_j: f64,
    maximum_chemical_energy_defect_j: f64,
    maximum_mass_residual_kg: f64,
    maximum_species_residual_kg: f64,
    maximum_pressure_residual_pa: f64,
    maximum_split_pressure_change_pa: f64,
    maximum_fraction_sum_error: f64,
}

fn inventory(room: &TwoZone) -> [f64; 5] {
    std::array::from_fn(|s| {
        room.layer(Layer::Lower).species_masses()[s] + room.layer(Layer::Upper).species_masses()[s]
    })
}

fn run(
    name: &'static str,
    width_m: f64,
    gas: &MixtureGas,
    burning: &InventoryLimitedReaction,
    air: &[f64],
    end_s: f64,
    maximum_dt_s: f64,
) -> Result<RunReport, Box<dyn Error>> {
    let room = TwoZone::new(
        gas.clone(),
        Area::new(ROOM_AREA_M2),
        Length::new(ROOM_HEIGHT_M),
        Length::new(0.0),
        Pressure::new(101325.0),
        [Temperature::new(AMBIENT_K); 2],
        [air.to_vec(), air.to_vec()],
        || false,
    )?;
    let initial_species = inventory(&room);
    let initial_energy = room.physical_internal_energy().value();
    let discharge = DischargeCoefficient::new(
        Dimensionless::new(0.7),
        SourceProvenance::new(
            "NIST TN 1889v1, section 4.1; nominal CFAST wall-opening Cd=0.7",
            "NIST-TN-1889v1:4.1",
        ),
        0.0,
        ToleranceBasis::EngineeringAllowance,
    )?;
    let mut network = ZoneNetwork::new(
        vec![Room {
            gas: room,
            floor_elevation: Length::new(0.0),
        }],
        Ambient {
            pressure: Pressure::new(101325.0),
            reference_height: Length::new(0.0),
            temperature: Temperature::new(AMBIENT_K),
            mass_fractions: air.to_vec(),
        },
        vec![WallVent {
            left: Endpoint::Room(0),
            right: Endpoint::Ambient,
            opening: WallOpening::new(
                Length::new(0.0),
                Length::new(ROOM_HEIGHT_M),
                Length::new(width_m),
                discharge,
            )?,
        }],
        vec![],
        Acceleration::new(9.81),
        || false,
    )?;
    let controls = NetworkControl::default();
    let mut out = RunReport {
        name,
        duration_s: 0.0,
        steps: 0,
        halvings: 0,
        temperature_k: AMBIENT_K,
        pressure_pa: 101325.0,
        final_species_kg: initial_species,
        boundary_species_kg: [0.0; 5],
        species_balance_kg: [0.0; 5],
        element_balance_mol: [0.0; 4],
        injected_fuel_kg: 0.0,
        reacted_fuel_kg: 0.0,
        initial_oxygen_limited_fuel_kg: initial_species[4] * nasa9::MASSES[0]
            / (2.0 * nasa9::MASSES[4]),
        first_oxygen_limited_s: None,
        last_released_power_w: 0.0,
        released_enthalpy_j: 0.0,
        physical_energy_balance_j: 0.0,
        maximum_chemical_energy_defect_j: 0.0,
        maximum_mass_residual_kg: 0.0,
        maximum_species_residual_kg: 0.0,
        maximum_pressure_residual_pa: 0.0,
        maximum_split_pressure_change_pa: 0.0,
        maximum_fraction_sum_error: 0.0,
    };
    let mut boundary_energy_j = 0.0;
    let mut total_extent_mol = 0.0;
    while network.rooms()[0].gas.time().value() < end_s {
        if out.steps >= 200_000 {
            return Err("step budget exhausted; increase the maximum timestep".into());
        }
        let remaining = end_s - network.rooms()[0].gas.time().value();
        // Reserve an admissible last interval instead of leaving a roundoff-sized
        // remainder below the network's declared minimum duration.
        let duration =
            if remaining > maximum_dt_s && remaining - maximum_dt_s < controls.minimum_step_s {
                0.5 * remaining
            } else {
                maximum_dt_s.min(remaining)
            };
        let (report, chemistry) = network.step_with_stages(
            Time::new(duration),
            controls,
            |rooms, _dt, rates, _cancelled| {
                // The inlet supplies original partial ENTHALPY, including the
                // feed's retained formation-energy gauge. It enters once.
                let feed = rooms[0].gas.inflow(
                    MassFlowRate::new(FUEL_RATE_KG_S),
                    &[1.0, 0.0, 0.0, 0.0, 0.0],
                    Temperature::new(AMBIENT_K),
                )?;
                rates[0][Layer::Upper.index()].add(&feed)?;
                Ok(AcceptedChemistry::default())
            },
            |rooms, dt, payload, additional, cancelled| {
                let room = &mut rooms[0].gas;
                let before = room.physical_internal_energy().value();
                // Availability is evaluated AFTER inflow and outflow. Inventory
                // that has left through a vent cannot also feed the reaction.
                let requested = room.layer(Layer::Upper).species_masses()[0] / nasa9::MASSES[0];
                let accepted = burning.apply_split(
                    room,
                    Layer::Upper,
                    Amount::new(requested),
                    dt,
                    Temperature::new(AMBIENT_K),
                    controls.zone,
                    &mut *cancelled,
                )?;
                let source = accepted.burning.chemical_rates;
                additional.add(&source)?;
                payload.extent_mol = accepted.burning.actual_extent.value();
                payload.released_enthalpy_j = accepted.burning.released_enthalpy.value();
                payload.released_power_w = accepted.burning.released_power.value();
                payload.oxygen_limited = accepted.burning.limiting_species == Some(4);
                for s in 0..5 {
                    payload.chemical_species_change_kg[s] = dt.value() * source.species_kg_s[s];
                }
                payload.chemical_physical_source_j = dt.value()
                    * (source.referenced_energy_w
                        + source.species_kg_s[..5]
                            .iter()
                            .zip(room.gas().energy_references())
                            .map(|(m, e)| m * e)
                            .sum::<f64>());
                payload.chemical_physical_energy_defect_j = room.physical_internal_energy().value()
                    - before
                    - payload.chemical_physical_source_j;
                Ok(())
            },
            || false,
        )?;
        out.steps += 1;
        out.halvings += report.halvings;
        out.injected_fuel_kg += FUEL_RATE_KG_S * report.elapsed.value();
        total_extent_mol += chemistry.extent_mol;
        out.released_enthalpy_j += chemistry.released_enthalpy_j;
        out.last_released_power_w = chemistry.released_power_w;
        boundary_energy_j +=
            report.external_physical_energy_j - chemistry.chemical_physical_source_j;
        for s in 0..5 {
            // Network's aggregate source includes chemistry. Remove that exact
            // source here to retain actual external species for the atom audit.
            out.boundary_species_kg[s] +=
                report.external_species_change_kg[s] - chemistry.chemical_species_change_kg[s];
        }
        let room = &network.rooms()[0].gas;
        if chemistry.oxygen_limited && out.first_oxygen_limited_s.is_none() {
            out.first_oxygen_limited_s = Some(room.time().value());
        }
        out.maximum_chemical_energy_defect_j = out
            .maximum_chemical_energy_defect_j
            .max(chemistry.chemical_physical_energy_defect_j.abs());
        out.maximum_mass_residual_kg = out
            .maximum_mass_residual_kg
            .max(report.mass_residual_kg.abs());
        out.maximum_species_residual_kg = out.maximum_species_residual_kg.max(
            report.species_residual_kg[..5]
                .iter()
                .map(|x| x.abs())
                .fold(0.0, f64::max),
        );
        out.maximum_pressure_residual_pa = out
            .maximum_pressure_residual_pa
            .max(report.maximum_pressure_residual_pa);
        out.maximum_split_pressure_change_pa = out
            .maximum_split_pressure_change_pa
            .max(report.maximum_split_pressure_change_pa);
        let active = room.layer(Layer::Upper);
        if active
            .species_masses()
            .iter()
            .any(|m| !m.is_finite() || *m < 0.0)
            || room.layer(Layer::Lower).mass_kg() != 0.0
        {
            return Err("invalid well-mixed species inventory".into());
        }
        let sum_y = active
            .species_masses()
            .iter()
            .map(|m| m / active.mass_kg())
            .sum::<f64>();
        out.maximum_fraction_sum_error = out.maximum_fraction_sum_error.max((sum_y - 1.0).abs());
    }
    let room = &network.rooms()[0].gas;
    out.duration_s = room.time().value();
    out.temperature_k = room.layer(Layer::Upper).temperature().value();
    out.pressure_pa = room.pressure().value();
    out.final_species_kg = inventory(room);
    out.reacted_fuel_kg = total_extent_mol * nasa9::MASSES[0];
    out.physical_energy_balance_j =
        room.physical_internal_energy().value() - initial_energy - boundary_energy_j;
    for s in 0..5 {
        out.species_balance_kg[s] = out.final_species_kg[s]
            - initial_species[s]
            - out.boundary_species_kg[s]
            - nasa9::STOICHIOMETRY[s] * nasa9::MASSES[s] * total_extent_mol;
    }
    let elemental = gas.source().elemental_matrix();
    for e in 0..4 {
        out.element_balance_mol[e] = (0..5)
            .map(|s| {
                (out.final_species_kg[s] - initial_species[s] - out.boundary_species_kg[s])
                    / nasa9::MASSES[s]
                    * elemental.get(e, s).unwrap() as f64
            })
            .sum();
    }
    if width_m == 0.0 {
        let expected = out.injected_fuel_kg.min(out.initial_oxygen_limited_fuel_kg);
        if (out.reacted_fuel_kg - expected).abs() > 1e-10 {
            return Err("closed-room stoichiometric oxygen/fuel limit failed".into());
        }
    }
    if out.physical_energy_balance_j.abs() > 0.01
        || out.species_balance_kg.iter().any(|v| v.abs() > 1e-9)
        || out.element_balance_mol.iter().any(|v| v.abs() > 1e-7)
        || out.maximum_fraction_sum_error > 1e-12
    {
        return Err("conserved reacting-room balance tolerance exceeded".into());
    }
    Ok(out)
}

impl RunReport {
    fn print_json(&self) {
        println!("{{");
        println!(
            "  \"case\": \"{}\", \"duration_s\": {:.12}, \"steps\": {}, \"halvings\": {},",
            self.name, self.duration_s, self.steps, self.halvings
        );
        println!(
            "  \"temperature_k\": {:.12}, \"pressure_pa\": {:.12},",
            self.temperature_k, self.pressure_pa
        );
        println!("  \"final_species_kg\": {:?},", self.final_species_kg);
        println!(
            "  \"net_boundary_species_kg\": {:?},",
            self.boundary_species_kg
        );
        println!("  \"species_balance_kg\": {:?},", self.species_balance_kg);
        println!("  \"element_balance_mol\": {:?},", self.element_balance_mol);
        println!(
            "  \"injected_fuel_kg\": {:.12}, \"reacted_fuel_kg\": {:.12},",
            self.injected_fuel_kg, self.reacted_fuel_kg
        );
        println!(
            "  \"initial_oxygen_limited_fuel_kg\": {:.12},",
            self.initial_oxygen_limited_fuel_kg
        );
        match self.first_oxygen_limited_s {
            Some(t) => println!("  \"first_oxygen_limited_s\": {t:.12},"),
            None => println!("  \"first_oxygen_limited_s\": null,"),
        }
        println!(
            "  \"last_released_power_w\": {:.12}, \"released_enthalpy_j\": {:.12},",
            self.last_released_power_w, self.released_enthalpy_j
        );
        println!(
            "  \"physical_energy_balance_j\": {:.12e},",
            self.physical_energy_balance_j
        );
        println!(
            "  \"maximum_chemical_energy_defect_j\": {:.12e},",
            self.maximum_chemical_energy_defect_j
        );
        println!(
            "  \"maximum_mass_residual_kg\": {:.12e},",
            self.maximum_mass_residual_kg
        );
        println!(
            "  \"maximum_species_residual_kg\": {:.12e},",
            self.maximum_species_residual_kg
        );
        println!(
            "  \"maximum_pressure_residual_pa\": {:.12e},",
            self.maximum_pressure_residual_pa
        );
        println!(
            "  \"maximum_split_pressure_change_pa\": {:.12e},",
            self.maximum_split_pressure_change_pa
        );
        println!(
            "  \"maximum_fraction_sum_error\": {:.12e}",
            self.maximum_fraction_sum_error
        );
        println!("}}");
    }
}

fn positive(value: &str) -> Result<f64, Box<dyn Error>> {
    let parsed: f64 = value.parse()?;
    if !parsed.is_finite() || parsed <= 0.0 {
        return Err("arguments must be finite and positive".into());
    }
    Ok(parsed)
}

fn main() -> Result<(), Box<dyn Error>> {
    let args = std::env::args().skip(1).collect::<Vec<_>>();
    if args.len() > 2 {
        return Err("usage: ventilated_compartment_fire [end-s=60] [maximum-dt-s=0.02]".into());
    }
    let end_s = args
        .first()
        .map(|s| positive(s))
        .transpose()?
        .unwrap_or(60.0);
    let maximum_dt_s = args
        .get(1)
        .map(|s| positive(s))
        .transpose()?
        .unwrap_or(0.02);
    if end_s > 120.0 || !(0.001..=0.25).contains(&maximum_dt_s) {
        return Err("supported demonstration: end <=120 s, maximum dt 0.001..=0.25 s".into());
    }
    let (gas, mechanism, air) = nasa9::source()?;
    if !gas
        .source()
        .models()
        .iter()
        .zip(nasa9::NAMES)
        .all(|(m, n)| m.species().as_str() == n)
    {
        return Err("source species axis differs from the authored reaction".into());
    }
    let burning =
        InventoryLimitedReaction::new(&gas, mechanism, &ReactionId::new(nasa9::REACTION)?, || {
            false
        })?;
    let closed = run("closed", 0.0, &gas, &burning, &air, end_s, maximum_dt_s)?;
    let ventilated = run(
        "doorway-ventilated",
        0.2,
        &gas,
        &burning,
        &air,
        end_s,
        maximum_dt_s,
    )?;
    println!("{{");
    println!(
        "\"model\": \"estimated well-mixed ideal-gas inventory-limited complete methane oxidation\","
    );
    println!(
        "\"source_commit\": \"{}\", \"source_blob\": \"{}\",",
        nasa9::SOURCE_COMMIT,
        nasa9::SOURCE_BLOB
    );
    println!(
        "\"species\": {:?}, \"elements\": [\"C\", \"H\", \"N\", \"O\"],",
        nasa9::NAMES
    );
    println!("\"maximum_dt_s\": {maximum_dt_s}, \"fuel_rate_kg_s\": {FUEL_RATE_KG_S},");
    println!("\"room_volume_m3\": 1.0, \"ambient_temperature_k\": {AMBIENT_K},");
    println!("\"caloric_anchor_k\": 300.0, \"high_region_corrections\": [");
    for (i, model) in gas.source().models().iter().enumerate() {
        let correction = &model
            .caloric_reference()
            .ok_or("missing explicit caloric reference")?
            .regions()[1];
        println!(
            "{{\"species\": \"{}\", \"enthalpy_j_mol\": {:.12e}, \"entropy_j_mol_k\": {:.12e}}}{}",
            model.species().as_str(),
            correction.enthalpy_offset().value(),
            correction.entropy_offset().value(),
            if i + 1 < 5 { "," } else { "" }
        );
    }
    println!("], \"cases\": [");
    closed.print_json();
    println!(",");
    ventilated.print_json();
    println!("]}}");
    Ok(())
}
