//! Data-only material substitution in native buoyant Stokes/enthalpy evolution.
//!
//! Run with --features enthalpy-coupling. The same source card supplies local
//! density, viscosity, conductivity and caloric response. This coarse closed
//! box demonstrates composition, not a calibrated liquid or converged CFD case.

use std::{error::Error, io::Read};

use fs_conduction::transient::deposition::*;
use fs_exec::CancelGate;
use fs_lbm::conjugate::{FvBoundary, SimpleConfig, ViscousHeatControls, VoxelDomain};
use fs_lbm::conjugate::enthalpy::advection::liquid::LiquidFlowEnthalpyInput;
use fs_lbm::material::transport::boussinesq::BoussinesqConfig;
use fs_matdb::{PropertyKey, SelectionPolicy};
use fs_material::liquid_temperature::{LiquidTemperatureLaw, LiquidViscosityBasis};
use fs_material::phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection,
    SPECIFIC_HEAT_CAPACITY_DIMS, SPECIFIC_HEAT_CAPACITY_PROPERTY,
};
use fs_qty::Density;

#[path = "material_buoyant_flow/material.rs"]
mod material;
use material::{MaterialInput, point};

#[derive(Clone, Copy, Debug)]
struct Observation {
    time_s: f64,
    maximum_speed_m_s: f64,
    minimum_temperature_k: f64,
    maximum_temperature_k: f64,
    thermal_change_j: f64,
    physical_heat_j: f64,
    kinetic_energy_j: f64,
    numerical_loss_j: f64,
    buoyancy_work_j: f64,
    combined_defect_j: f64,
}

fn initial<'a>(
    domain: &VoxelDomain, input: &MaterialInput, caloric: &'a ResolvedHeatCapacityCurve,
    hot_k: f64,
) -> Result<DepositionSimulation<'a>, Box<dyn Error>> {
    let mut events = Vec::new();
    for c in 0..4 {
        let [x, y, z] = domain.coords(c);
        let centre = domain.center(x, y, z);
        let mut start = centre; let mut end = centre;
        start[0] -= 0.5; end[0] += 0.5;
        events.push(TimedEvent { time_s: 0.0, event: Event::Birth(DepositedCell {
            id: CellId(c as u64), segment_id: c as u64, body_id: 1,
            material_id: MaterialId(1), mass_kg: input.reference_density * domain.volume(c),
            specific_h_j_kg: if x == 0 { 0.0 } else { input.cp * (hot_k - input.temperatures[0]) },
            geometry: CellGeometry { start_m: start, end_m: end,
                volume_m3: domain.volume(c), surface_area_m2: 6.0 },
            boundary: Boundary { ambient_k: input.temperatures[0], htc_w_m2_k: 0.0,
                reservoir: None, source_w: 0.0 },
        }) });
    }
    for (i, [a, b]) in [[0, 1], [2, 3], [0, 2], [1, 3]].into_iter().enumerate() {
        events.push(TimedEvent { time_s: 0.0, event: Event::Contact(ThermalLink {
            id: 100 + i as u64, a: CellId(a), b: CellId(b), area_m2: 1.0,
            // The selected source owner replaces this descriptor at every epoch.
            conductance_w_k: 0.0,
        }) });
    }
    let state = DepositionSimulation::new(events,
        &[(MaterialId(1), caloric.curve() as &dyn ThermalLaw)],
        DepositionBudget { max_cells: 4, max_links: 4, max_events: 8, max_substeps: 100_000 })?;
    Ok(state.advance_to(0.0, thermal_step(), || Ok(()))?.0)
}

fn thermal_step() -> StepConfig {
    StepConfig { max_dt_s: 0.05, energy_tolerance_j: 1e-8 }
}

fn simulate(
    input: &MaterialInput, steps: usize, hot_k: f64, mut emit: impl FnMut(Observation),
) -> Result<Observation, Box<dyn Error>> {
    if !(1..=2000).contains(&steps) || !hot_k.is_finite()
        || hot_k <= input.temperatures[0] || hot_k >= input.temperatures[1] {
        return Err("need 1..2000 steps and a hot temperature strictly inside the material span".into());
    }
    let reference = point(input.temperatures[0]);
    let upper = point(input.temperatures[1]);
    let liquid = LiquidTemperatureLaw::try_new(&input.card, &reference, "T",
        LiquidViscosityBasis::Dynamic, MaterialPropertySelection::SingleClaimOnly)?;
    let cp = PropertyKey::new(SPECIFIC_HEAT_CAPACITY_PROPERTY, SPECIFIC_HEAT_CAPACITY_DIMS);
    let rho = PropertyKey::new(DENSITY_PROPERTY, Density::DIMS);
    let caloric = ResolvedHeatCapacityCurve::try_new(&input.card, HeatCapacitySourceQuery {
        heat_capacity: &cp, density: &rho, temperature_axis: "T", lower: &reference,
        upper: &upper, selection: SelectionPolicy::SingleClaimOnly,
    }, SolidLiquidPhase::Liquid, 0.0, 1e-8)?;
    let domain = VoxelDomain::new(2, 2, 1, 1.0)?;
    let mut thermal = initial(&domain, input, &caloric, hot_k)?;
    let starting_energy = thermal.energy_j();
    let starting_mass = thermal.mass_kg();
    let ids = [Some(CellId(0)), Some(CellId(1)), Some(CellId(2)), Some(CellId(3))];
    let mut faces = [FvBoundary::wall(); 6];
    faces[4] = FvBoundary::Symmetry; faces[5] = FvBoundary::Symmetry;
    let mut flow = SimpleConfig::new(faces);
    flow.tolerance = 1e-13; flow.pressure_tolerance = 1e-12;
    let gate = CancelGate::new();
    let assignment = LiquidFlowEnthalpyInput {
        cell_ids: &ids, material_id: MaterialId(1), liquid: &liquid, caloric: &caloric,
        inlet_specific_h_j_kg: [None; 6], max_cells: 4, max_connections: 8,
    };
    let mut observer = assignment.feedback_boussinesq(&thermal, &domain, &flow, (), &gate, 100_000,
        BoussinesqConfig { gravity_m_s2: [0.0, -9.81, 0.0],
            max_relative_density_departure: 0.05, max_relative_specific_heat_departure: 0.0 })?
        .with_fluid_conduction(&thermal, &[100, 101, 102, 103])?
        .with_unsteady_viscous_heating(ViscousHeatControls {
            maximum_reynolds: 0.1, power_tolerance_w: 1e-7,
        })?;
    let (mut physical_heat, mut numerical_loss, mut buoyancy_work, mut supplied_work) = (0.0, 0.0, 0.0, 0.0);
    let mut last = None;
    for step in 1..=steps {
        let time = step as f64 * thermal_step().max_dt_s;
        let (next, accepted, heat) = thermal.advance_observed(&observer, time, thermal_step(), || Ok(()))?;
        let work = accepted.unsteady_heating_audit().ok_or("missing accepted mechanical account")?;
        physical_heat += work.generated_heat_j;
        numerical_loss += work.backward_euler_loss_j;
        buoyancy_work += work.supplied_buoyancy_work_j;
        // Buoyancy is a SUBSET of body work; never add the two together.
        supplied_work += work.supplied_wall_work_j + work.supplied_body_work_j + work.pressure_constraint_work_j;
        if (heat.external_heat_j - work.generated_heat_j).abs() > 1e-8
            || next.mass_kg().to_bits() != starting_mass.to_bits() {
            return Err("the original thermal or mass account did not close".into());
        }
        let temperatures = (0..4).map(|c| next.temperature_k(CellId(c)))
            .collect::<Result<Vec<_>, _>>()?;
        let maximum_speed = accepted.unsteady_flow().ok_or("missing accepted momentum")?
            .velocity_m_s.iter().map(|v| v.iter().map(|x| x*x).sum::<f64>().sqrt())
            .fold(0.0_f64, f64::max);
        let thermal_change = next.energy_j() - starting_energy;
        let observation = Observation {
            time_s: next.time_s(), maximum_speed_m_s: maximum_speed,
            minimum_temperature_k: temperatures.iter().copied().fold(f64::INFINITY, f64::min),
            maximum_temperature_k: temperatures.iter().copied().fold(f64::NEG_INFINITY, f64::max),
            thermal_change_j: thermal_change, physical_heat_j: physical_heat,
            kinetic_energy_j: work.final_kinetic_energy_j, numerical_loss_j: numerical_loss,
            buoyancy_work_j: buoyancy_work,
            combined_defect_j: thermal_change + work.final_kinetic_energy_j + numerical_loss - supplied_work,
        };
        if observation.combined_defect_j.abs() > 1e-6 {
            return Err("thermal plus kinetic plus numerical-loss account exceeded 1e-6 J".into());
        }
        emit(observation); last = Some(observation);
        thermal = next; observer = accepted;
    }
    Ok(last.expect("positive admitted step count"))
}

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if !(1..=3).contains(&args.len()) {
        return Err("usage: material_buoyant_flow MATERIAL.txt [STEPS=40] [HOT_K=400]".into());
    }
    let steps = args.get(1).map_or(Ok(40), |s| s.parse())?;
    let hot_k = args.get(2).map_or(Ok(400.0), |s| s.parse())?;
    let mut text = String::new();
    std::fs::File::open(&args[0])?.take(65_537).read_to_string(&mut text)?;
    let input = MaterialInput::parse(&text)?;
    println!("# synthetic material card; 2x2x1 unit-metre cells; closed low-Re Boussinesq/Stokes model");
    println!("# buoyancy work is an explicit mechanical exchange; no gravitational-potential store is modeled");
    println!("time_s,max_speed_m_s,min_K,max_K,thermal_change_J,physical_heat_J,kinetic_J,numerical_loss_J,buoyancy_work_J,combined_defect_J");
    simulate(&input, steps, hot_k, |r| println!(
        "{:.6},{:.9e},{:.9},{:.9},{:.9e},{:.9e},{:.9e},{:.9e},{:.9e},{:.3e}",
        r.time_s, r.maximum_speed_m_s, r.minimum_temperature_k, r.maximum_temperature_k,
        r.thermal_change_j, r.physical_heat_j, r.kinetic_energy_j, r.numerical_loss_j,
        r.buoyancy_work_j, r.combined_defect_j,
    ))?;
    Ok(())
}

#[cfg(test)]
#[path = "material_buoyant_flow/tests.rs"]
mod tests;
