//! A complete water/steam Rankine cycle using the original IF97 caloric source.
//!
//! Run: `cargo run -p fs-material --example water_rankine`.
//! A pump raises saturated liquid from 10 kPa to 4 MPa; a constant-pressure
//! boiler heats it to 773.15 K; a turbine expands it back to 10 kPa; a condenser
//! returns the wet exhaust to saturated liquid. All amounts are per kilogram
//! of steady throughput. The pump and turbine are adiabatic, kinetic/potential
//! energy changes are neglected, and the heat exchangers have no pressure drop.
//! The declared efficiencies are illustrative inputs, not calibrated machines.
//!
//! Pressure-changing devices use raw IF97 entropy. The source's small regional
//! saturation matching defects remain visible and are not machine entropy
//! generation. Heat-exchanger reservoir temperatures are not supplied, so this
//! example reports fluid entropy changes without claiming a complete reservoir
//! entropy budget. Flow work is already in h=u+p*v and is never added to shaft
//! work again.

use fs_material::phase::water::flow::{WaterFlowSolve, WaterShaftDevice, WaterShaftSolve};
use fs_material::phase::water::{
    IsobaricWater, WaterChartControls, WaterEntropySolve, WaterError, WaterState,
};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::water::if97::bounds::If97HeatCapacityBoundSolve;
use fs_thermochem::water::if97::{SOURCE_RELEASE, SOURCE_URL};

const CONDENSER_PA: f64 = 10_000.0;
const BOILER_PA: f64 = 4_000_000.0;
const BOILER_OUTLET_K: f64 = 773.15;
const PUMP_EFFICIENCY: f64 = 0.8;
const TURBINE_EFFICIENCY: f64 = 0.85;
const DEVICE_ENERGY_TOLERANCE_J_KG: f64 = 2e-7;
const CYCLE_ENERGY_TOLERANCE_J_KG: f64 = 8e-7;

fn chart(pressure_pa: f64) -> Result<IsobaricWater, WaterError> {
    IsobaricWater::new(
        Pressure::new(pressure_pa),
        WaterChartControls {
            minimum_temperature_k: 273.16,
            maximum_temperature_k: 1073.15,
            enthalpy_tolerance_j_kg: 1e-7,
            maximum_iterations: 96,
            maximum_entropy_seam_j_kg_k: 0.1,
            heat_capacity: If97HeatCapacityBoundSolve {
                maximum_intervals: 8191,
            },
        },
        || false,
    )
}

fn shaft(efficiency: f64) -> WaterShaftSolve {
    WaterShaftSolve {
        isentropic_efficiency: efficiency,
        entropy: WaterEntropySolve {
            entropy_tolerance_j_kg_k: 1e-9,
            maximum_iterations: 96,
        },
        balance: WaterFlowSolve {
            energy_tolerance_j_kg: DEVICE_ENERGY_TOLERANCE_J_KG,
            entropy_tolerance_j_kg_k: 2e-9,
        },
    }
}

fn row(name: &str, state: WaterState) {
    let quality = state
        .vapor_quality()
        .map_or_else(|| "n/a".to_owned(), |value| format!("{value:.9}"));
    println!(
        "{name},{:?},{:.6},{:.6},{:.9},{:.9},{:.9},{:.12},{quality}",
        state.phase(),
        state.pressure_pa(),
        state.temperature_k(),
        state.specific_enthalpy_j_kg(),
        state.specific_internal_energy_j_kg(),
        state.specific_entropy_j_kg_k(),
        state.specific_volume_m3_kg(),
    );
}

fn pv(state: WaterState) -> f64 {
    state.pressure_pa() * state.specific_volume_m3_kg()
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let low = chart(CONDENSER_PA)?;
    let high = chart(BOILER_PA)?;
    let liquid = low
        .state_at_specific_enthalpy(
            low.saturation().liquid().specific_enthalpy().value(),
            || false,
        )?
        .into_state();
    let pump = liquid.adiabatic_shaft_to(
        &high,
        WaterShaftDevice::Compressor,
        shaft(PUMP_EFFICIENCY),
        || false,
    )?;
    let steam = high.state_at_temperature(Temperature::new(BOILER_OUTLET_K), || false)?;
    let turbine = steam.adiabatic_shaft_to(
        &low,
        WaterShaftDevice::Expander,
        shaft(TURBINE_EFFICIENCY),
        || false,
    )?;

    let heat_in = steam.specific_enthalpy_j_kg() - pump.outlet.specific_enthalpy_j_kg();
    let heat_out = liquid.specific_enthalpy_j_kg() - turbine.outlet.specific_enthalpy_j_kg();
    let net_work = pump.balance.shaft_work_j_kg + turbine.balance.shaft_work_j_kg;
    let efficiency = net_work / heat_in;
    if !(heat_in > 0.0 && heat_out < 0.0 && net_work > 0.0 && efficiency < 1.0) {
        return Err("Rankine cycle did not produce admitted positive net shaft work".into());
    }

    // Each row independently resolves Δu + Δ(p*v). Summing these changes must
    // close around the actual four-state cycle, while Q_in+Q_out equals W_net.
    let stages = [
        (
            "pump",
            liquid,
            pump.outlet,
            0.0,
            pump.balance.shaft_work_j_kg,
        ),
        ("boiler", pump.outlet, steam, heat_in, 0.0),
        (
            "turbine",
            steam,
            turbine.outlet,
            0.0,
            turbine.balance.shaft_work_j_kg,
        ),
        ("condenser", turbine.outlet, liquid, heat_out, 0.0),
    ];
    let mut cycle_dh = 0.0;
    let mut cycle_du = 0.0;
    let mut cycle_dpv = 0.0;
    let mut cycle_ds = 0.0;
    let mut maximum_residual = 0.0_f64;
    println!("# {SOURCE_RELEASE}: {SOURCE_URL}");
    println!(
        "# Steady flow; heat positive into water; shaft work positive outward; all amounts per kg."
    );
    println!(
        "# eta_pump={PUMP_EFFICIENCY}; eta_turbine={TURBINE_EFFICIENCY}; no calibrated machine model."
    );
    println!(
        "state,phase,pressure_Pa,temperature_K,enthalpy_J_kg,internal_energy_J_kg,source_entropy_J_kg_K,specific_volume_m3_kg,vapor_quality"
    );
    row("condenser_outlet", liquid);
    row("pump_outlet", pump.outlet);
    row("boiler_outlet", steam);
    row("turbine_outlet", turbine.outlet);
    println!(
        "stage,heat_J_kg,shaft_work_J_kg,delta_h_J_kg,delta_u_J_kg,delta_pv_J_kg,delta_source_s_J_kg_K,h_balance_J_kg,u_pv_balance_J_kg"
    );
    for (name, before, after, heat, work) in stages {
        let dh = after.specific_enthalpy_j_kg() - before.specific_enthalpy_j_kg();
        let du = after.specific_internal_energy_j_kg() - before.specific_internal_energy_j_kg();
        let dpv = pv(after) - pv(before);
        let ds = after.specific_entropy_j_kg_k() - before.specific_entropy_j_kg_k();
        let h_residual = dh + work - heat;
        let u_residual = du + dpv + work - heat;
        for residual in [h_residual, u_residual] {
            maximum_residual = maximum_residual.max(residual.abs());
            if !residual.is_finite() || residual.abs() > DEVICE_ENERGY_TOLERANCE_J_KG {
                return Err(format!("{name} energy residual {residual} J/kg").into());
            }
        }
        cycle_dh += dh;
        cycle_du += du;
        cycle_dpv += dpv;
        cycle_ds += ds;
        println!(
            "{name},{heat:.9},{work:.9},{dh:.9},{du:.9},{dpv:.9},{ds:.9},{h_residual:.9e},{u_residual:.9e}"
        );
    }
    let heat_work_residual = heat_in + heat_out - net_work;
    let first_law_residual = cycle_du + cycle_dpv + net_work - heat_in - heat_out;
    for residual in [
        cycle_dh,
        cycle_du,
        cycle_dpv,
        heat_work_residual,
        first_law_residual,
    ] {
        if !residual.is_finite() || residual.abs() > CYCLE_ENERGY_TOLERANCE_J_KG {
            return Err(format!("closed cycle energy residual {residual} J/kg").into());
        }
    }
    println!(
        "# heat_in={heat_in:.9} J/kg; heat_out={heat_out:.9} J/kg; net_shaft_work={net_work:.9} J/kg; efficiency={efficiency:.9}"
    );
    println!(
        "# cycle_delta_h={cycle_dh:.9e} J/kg; cycle_delta_u={cycle_du:.9e} J/kg; cycle_delta_pv={cycle_dpv:.9e} J/kg"
    );
    println!(
        "# heat_work_residual={heat_work_residual:.9e} J/kg; first_law_residual={first_law_residual:.9e} J/kg; maximum_stage_residual={maximum_residual:.9e} J/kg"
    );
    println!(
        "# raw_fluid_entropy_cycle_change={cycle_ds:.9e} J/(kg*K); pump_generation={:.9}; turbine_generation={:.9}",
        pump.balance.source_entropy_change_j_kg_k, turbine.balance.source_entropy_change_j_kg_k
    );
    println!(
        "# source_saturation_entropy_defect: condenser={:.9e} J/(kg*K); boiler={:.9e} J/(kg*K)",
        low.saturation_entropy_defect_j_kg_k(),
        high.saturation_entropy_defect_j_kg_k()
    );
    println!(
        "# Source defects are separate from entropy generation; heat-exchanger reservoir entropy is unspecified."
    );
    Ok(())
}
