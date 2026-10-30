//! Complete subcritical and supercritical Rankine cycles on original IF97.
//!
//! Run: `cargo run -p fs-material --example water_regional_rankine`.
//! A single prepared isobar spans cold liquid, Region 3 and hot Region 2 steam.
//! The two illustrated boiler pressures are 20 MPa and 25 MPa; both cycles
//! include the feed pump, boiler, turbine and complete condensation at 10 kPa.
//! All energy quantities are per kilogram of steady throughput. Heat exchangers
//! have no pressure drop; the pump and turbine are adiabatic, with explicitly
//! supplied efficiencies, and kinetic/potential changes are neglected.
//!
//! Original source h/u/s are unchanged. Finite IF97 regional matching defects
//! are reported separately from sensible and physical latent contributions.
//! They are not extra heat sources, an artificial phase plateau or machine
//! entropy production. Reservoir temperatures and equipment calibration are not
//! supplied; the efficiencies below describe these illustrative cycles only.

use fs_material::phase::water::WaterEntropySolve;
use fs_material::phase::water::flow::across::{WaterFlowChart, WaterFlowState};
use fs_material::phase::water::flow::{WaterFlowSolve, WaterShaftDevice, WaterShaftSolve};
use fs_material::phase::water::regional::{
    RegionalJunctionControls, RegionalWater, RegionalWaterControls, RegionalWaterState,
};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::water::if97::region3::pressure::DensePressureSolve;
use fs_thermochem::water::if97::{SOURCE_RELEASE, SOURCE_URL};

const CONDENSER_PA: f64 = 10_000.0;
const BOILER_OUTLET_K: f64 = 823.15;
const PUMP_EFFICIENCY: f64 = 0.80;
const TURBINE_EFFICIENCY: f64 = 0.85;
const ENERGY_TOLERANCE_J_KG: f64 = 2e-5;

fn chart(pressure_pa: f64) -> Result<RegionalWater, Box<dyn std::error::Error>> {
    Ok(RegionalWater::new(
        Pressure::new(pressure_pa),
        RegionalWaterControls {
            minimum_temperature_k: 273.16,
            maximum_temperature_k: 873.15,
            enthalpy_tolerance_j_kg: 1e-6,
            maximum_iterations: 96,
            maximum_saturation_entropy_defect_j_kg_k: 0.1,
            pressure: DensePressureSolve {
                pressure_tolerance_pa: 1e-4,
                density_tolerance_kg_m3: 1e-8,
                maximum_iterations: 128,
            },
            junctions: RegionalJunctionControls {
                maximum_enthalpy_jump_j_kg: 200.0,
                maximum_entropy_jump_j_kg_k: 0.25,
                maximum_volume_jump_m3_kg: 1e-5,
                maximum_temperature_gap_k: 1e-8,
                maximum_iterations: 96,
            },
        },
        || false,
    )?)
}

fn shaft(efficiency: f64) -> WaterShaftSolve {
    WaterShaftSolve {
        isentropic_efficiency: efficiency,
        entropy: WaterEntropySolve {
            entropy_tolerance_j_kg_k: 1e-8,
            maximum_iterations: 96,
        },
        balance: WaterFlowSolve {
            energy_tolerance_j_kg: ENERGY_TOLERANCE_J_KG,
            entropy_tolerance_j_kg_k: 2e-8,
        },
    }
}

fn regional(state: WaterFlowState) -> Result<RegionalWaterState, Box<dyn std::error::Error>> {
    match state {
        WaterFlowState::Regional(state) => Ok(state),
        _ => Err("cycle endpoint lost its regional pressure chart".into()),
    }
}

fn row(label: &str, state: RegionalWaterState) {
    let source = state.single_phase_source().map_or_else(
        || "coexistence".to_owned(),
        |source| format!("{:?}", source.region()),
    );
    let quality = state
        .vapor_quality()
        .map_or_else(|| "n/a".to_owned(), |x| format!("{x:.9}"));
    println!(
        "{label},{source},{:.6},{:.9},{:.9},{:.9},{:.9},{:.12},{quality}",
        state.pressure_pa(),
        state.temperature_k(),
        state.specific_enthalpy_j_kg(),
        state.specific_internal_energy_j_kg(),
        state.specific_entropy_j_kg_k(),
        state.specific_volume_m3_kg(),
    );
}

fn run_cycle(boiler_pressure_pa: f64) -> Result<(), Box<dyn std::error::Error>> {
    let low = chart(CONDENSER_PA)?;
    let high = chart(boiler_pressure_pa)?;
    let liquid = low.state_at_quality(0.0, || false)?;
    let pump = WaterFlowState::from(liquid).adiabatic_shaft_to(
        &WaterFlowChart::from(high),
        WaterShaftDevice::Compressor,
        shaft(PUMP_EFFICIENCY),
        || false,
    )?;
    let feed = regional(pump.outlet)?;
    let target = high.state_at_temperature(Temperature::new(BOILER_OUTLET_K), || false)?;
    let heat_in = target.specific_enthalpy_j_kg() - feed.specific_enthalpy_j_kg();
    // One accepted heat operation traverses the source regions and any actual
    // saturation plateau. The caller does not split or reconcile the sources.
    let boiler = high.advance_heat(feed, 1.0, heat_in, ENERGY_TOLERANCE_J_KG, || false)?;
    let steam = boiler.state;
    let turbine = WaterFlowState::from(steam).adiabatic_shaft_to(
        &WaterFlowChart::from(low),
        WaterShaftDevice::Expander,
        shaft(TURBINE_EFFICIENCY),
        || false,
    )?;
    let exhaust = regional(turbine.outlet)?;
    let heat_out = liquid.specific_enthalpy_j_kg() - exhaust.specific_enthalpy_j_kg();
    let condenser = low.advance_heat(exhaust, 1.0, heat_out, ENERGY_TOLERANCE_J_KG, || false)?;
    let returned = condenser.state;
    let net_work = pump.balance.shaft_work_j_kg + turbine.balance.shaft_work_j_kg;
    let efficiency = net_work / heat_in;
    if !(heat_in > 0.0
        && heat_out < 0.0
        && net_work > 0.0
        && efficiency > 0.0
        && efficiency < 1.0
        && returned.vapor_quality() == Some(0.0))
    {
        return Err("Rankine heat, work or complete-condensation admission failed".into());
    }

    println!("# Complete Rankine cycle: boiler_pressure={boiler_pressure_pa:.0} Pa");
    println!(
        "state,source,pressure_Pa,temperature_K,enthalpy_J_kg,internal_energy_J_kg,source_entropy_J_kg_K,specific_volume_m3_kg,vapor_quality"
    );
    row("condenser_outlet", liquid);
    row("pump_outlet", feed);
    for (label, t) in [
        ("boiler_600K", 600.0),
        ("boiler_630K", 630.0),
        ("boiler_650K", 650.0),
    ] {
        row(
            label,
            high.state_at_temperature(Temperature::new(t), || false)?,
        );
    }
    row("boiler_outlet", steam);
    row("turbine_outlet", exhaust);
    row("returned_condensate", returned);

    println!(
        "junction,lower_temperature_K,upper_temperature_K,source_delta_h_J_kg,source_delta_s_J_kg_K,source_delta_v_m3_kg"
    );
    for junction in high.junctions() {
        let [lower, upper] = junction.temperature_bounds_k();
        println!(
            "{:?},{lower:.12},{upper:.12},{:.9},{:.9},{:.12e}",
            junction.kind(),
            junction.enthalpy_change().value(),
            junction.entropy_change().value(),
            junction.volume_change().value(),
        );
    }

    let stages = [
        ("pump", liquid, feed, 0.0, pump.balance.shaft_work_j_kg),
        ("boiler", feed, steam, heat_in, 0.0),
        (
            "turbine",
            steam,
            exhaust,
            0.0,
            turbine.balance.shaft_work_j_kg,
        ),
        ("condenser", exhaust, returned, heat_out, 0.0),
    ];
    let mut maximum_residual = 0.0_f64;
    let mut cycle = [0.0; 4];
    println!(
        "stage,heat_J_kg,shaft_work_J_kg,delta_h_J_kg,delta_u_J_kg,delta_pv_J_kg,delta_source_s_J_kg_K,h_balance_J_kg,u_pv_balance_J_kg"
    );
    for (label, before, after, heat, work) in stages {
        let dh = after.specific_enthalpy_j_kg() - before.specific_enthalpy_j_kg();
        let du = after.specific_internal_energy_j_kg() - before.specific_internal_energy_j_kg();
        let dpv = after.pressure_pa() * after.specific_volume_m3_kg()
            - before.pressure_pa() * before.specific_volume_m3_kg();
        let ds = after.specific_entropy_j_kg_k() - before.specific_entropy_j_kg_k();
        let h_residual = dh + work - heat;
        let u_residual = du + dpv + work - heat;
        for residual in [h_residual, u_residual] {
            maximum_residual = maximum_residual.max(residual.abs());
            if !residual.is_finite() || residual.abs() > ENERGY_TOLERANCE_J_KG {
                return Err(format!("{label} energy residual {residual} J/kg").into());
            }
        }
        for (sum, change) in cycle.iter_mut().zip([dh, du, dpv, ds]) {
            *sum += change;
        }
        println!(
            "{label},{heat:.9},{work:.9},{dh:.9},{du:.9},{dpv:.9},{ds:.9},{h_residual:.9e},{u_residual:.9e}"
        );
    }
    let heat_work_residual = heat_in + heat_out - net_work;
    let first_law_residual = cycle[1] + cycle[2] + net_work - heat_in - heat_out;
    for residual in [
        cycle[0],
        cycle[1],
        cycle[2],
        heat_work_residual,
        first_law_residual,
    ] {
        if !residual.is_finite() || residual.abs() > 4.0 * ENERGY_TOLERANCE_J_KG {
            return Err(format!("closed-cycle energy residual {residual} J/kg").into());
        }
    }
    if !cycle[3].is_finite() || cycle[3].abs() > 1e-7 {
        return Err("closed-cycle fluid entropy failed to return".into());
    }
    println!(
        "# boiler_sensible={:.9} J/kg; boiler_latent={:.9} J/kg; regional_matching={:.9} J/kg",
        boiler.sensible_enthalpy_change_j,
        boiler.latent_enthalpy_change_j,
        boiler.regional_enthalpy_matching_change_j,
    );
    println!(
        "# pump_input={:.9} J/kg; turbine_output={:.9} J/kg; net_work={net_work:.9} J/kg",
        -pump.balance.shaft_work_j_kg, turbine.balance.shaft_work_j_kg,
    );
    println!(
        "# heat_in={heat_in:.9} J/kg; heat_out={heat_out:.9} J/kg; cycle_efficiency={efficiency:.9}"
    );
    println!(
        "# heat_work_residual={heat_work_residual:.9e} J/kg; first_law_residual={first_law_residual:.9e} J/kg; maximum_stage_residual={maximum_residual:.9e} J/kg"
    );
    println!(
        "# fluid_entropy_cycle_change={:.9e} J/(kg*K); pump_entropy_generation={:.9}; turbine_entropy_generation={:.9}",
        cycle[3],
        pump.balance.source_entropy_change_j_kg_k,
        turbine.balance.source_entropy_change_j_kg_k,
    );
    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("# {SOURCE_RELEASE}: {SOURCE_URL}");
    println!("# Heat positive into water; shaft work positive outward; all amounts per kg.");
    println!(
        "# eta_pump={PUMP_EFFICIENCY}; eta_turbine={TURBINE_EFFICIENCY}; boiler_outlet={BOILER_OUTLET_K} K."
    );
    println!(
        "# Regional matching defects remain source-fit diagnostics; no reservoir entropy claim."
    );
    for pressure in [20e6, 25e6] {
        run_cycle(pressure)?;
    }
    Ok(())
}
