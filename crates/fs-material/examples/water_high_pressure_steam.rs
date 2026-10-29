//! High-pressure boiling, steam expansion and condensation on original IF97.
//!
//! Run: `cargo run -p fs-material --example water_high_pressure_steam`.
//! One kilogram of steady throughput enters at 20 MPa and 630 K, boils and
//! superheats to 645 K, expands through an 85% isentropic-efficiency turbine
//! to 10 kPa, then condenses to saturated liquid. The dense boiler uses Region
//! 3; the low-pressure condenser uses Regions 1/2/4. This is a partial steam
//! train: feed-pump/preheater and return to the inlet are not modeled.
//!
//! Heat exchangers have no pressure drop. The turbine is adiabatic, and all
//! stages neglect kinetic/potential energy changes. Efficiency is an explicit
//! illustrative input. Original h/u/s and finite pressure/entropy closure
//! defects remain visible. Flow work is already in h=u+p*v and is not added
//! again to shaft work. Reservoir temperatures and phase-change kinetics are
//! not specified, so no reservoir entropy or transient boiling claim is made.

use fs_material::phase::water::flow::across::{WaterFlowChart, WaterFlowState};
use fs_material::phase::water::flow::{WaterFlowSolve, WaterShaftDevice, WaterShaftSolve};
use fs_material::phase::water::pressure::{DenseIsobaricControls, DenseIsobaricWater};
use fs_material::phase::water::{IsobaricWater, WaterChartControls, WaterEntropySolve};
use fs_qty::{Pressure, Temperature};
use fs_thermochem::water::if97::bounds::If97HeatCapacityBoundSolve;
use fs_thermochem::water::if97::region3::pressure::DensePressureSolve;
use fs_thermochem::water::if97::{SOURCE_RELEASE, SOURCE_URL};

const BOILER_PRESSURE_PA: f64 = 20e6;
const CONDENSER_PRESSURE_PA: f64 = 10_000.0;
const TURBINE_EFFICIENCY: f64 = 0.85;
const ENERGY_TOLERANCE_J_KG: f64 = 2e-5;

fn row(label: &str, state: WaterFlowState) {
    let source = match state {
        WaterFlowState::Gibbs(_) => "Gibbs",
        WaterFlowState::Dense(_) => "Helmholtz",
    };
    let quality = state
        .vapor_quality()
        .map_or_else(|| "n/a".to_owned(), |value| format!("{value:.9}"));
    println!(
        "{label},{source},{:.6},{:.9},{:.9},{:.9},{:.9},{:.12},{quality},{:.9e}",
        state.pressure_pa(),
        state.temperature_k(),
        state.specific_enthalpy_j_kg(),
        state.specific_internal_energy_j_kg(),
        state.specific_entropy_j_kg_k(),
        state.specific_volume_m3_kg(),
        state.pressure_closure_defect_j_kg(),
    );
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let boiler = DenseIsobaricWater::new(
        Pressure::new(BOILER_PRESSURE_PA),
        DenseIsobaricControls {
            minimum_temperature_k: 625.0,
            maximum_temperature_k: 646.0,
            enthalpy_tolerance_j_kg: 1e-5,
            maximum_iterations: 96,
            maximum_entropy_seam_j_kg_k: 0.1,
            pressure: DensePressureSolve {
                pressure_tolerance_pa: 1e-4,
                density_tolerance_kg_m3: 1e-8,
                maximum_iterations: 128,
            },
        },
        || false,
    )?;
    let condenser = IsobaricWater::new(
        Pressure::new(CONDENSER_PRESSURE_PA),
        WaterChartControls {
            minimum_temperature_k: 273.16,
            maximum_temperature_k: 1073.15,
            enthalpy_tolerance_j_kg: 1e-6,
            maximum_iterations: 96,
            maximum_entropy_seam_j_kg_k: 0.1,
            heat_capacity: If97HeatCapacityBoundSolve {
                maximum_intervals: 8191,
            },
        },
        || false,
    )?;
    let inlet = boiler.state_at_temperature(Temperature::new(630.0), || false)?;
    let liquid = boiler.state_at_quality(0.0, || false)?;
    let vapor = boiler.state_at_quality(1.0, || false)?;
    let steam = boiler.state_at_temperature(Temperature::new(645.0), || false)?;
    let targets = [liquid, vapor, steam];
    let mut boiler_states = [inlet; 4];
    let mut heats = [0.0; 3];
    for (i, target) in targets.into_iter().enumerate() {
        let before = boiler_states[i];
        let heat = target.specific_enthalpy_j_kg() - before.specific_enthalpy_j_kg();
        let heated = boiler.advance_heat(before, 1.0, heat, ENERGY_TOLERANCE_J_KG, || false)?;
        boiler_states[i + 1] = heated.state;
        heats[i] = heat;
    }
    let turbine_inlet = WaterFlowState::from(boiler_states[3]);
    let turbine = turbine_inlet.adiabatic_shaft_to(
        &WaterFlowChart::from(condenser),
        WaterShaftDevice::Expander,
        WaterShaftSolve {
            isentropic_efficiency: TURBINE_EFFICIENCY,
            entropy: WaterEntropySolve {
                entropy_tolerance_j_kg_k: 1e-8,
                maximum_iterations: 96,
            },
            balance: WaterFlowSolve {
                energy_tolerance_j_kg: ENERGY_TOLERANCE_J_KG,
                entropy_tolerance_j_kg_k: 2e-8,
            },
        },
        || false,
    )?;
    let WaterFlowState::Gibbs(wet) = turbine.outlet else {
        return Err("turbine outlet did not retain its low-pressure source".into());
    };
    let final_h = condenser.saturation().liquid().specific_enthalpy().value();
    let condenser_heat = final_h - wet.specific_enthalpy_j_kg();
    let cooled =
        condenser.advance_heat(wet, 1.0, condenser_heat, ENERGY_TOLERANCE_J_KG, || false)?;
    let final_state = WaterFlowState::from(cooled.state);
    if !(heats.into_iter().all(|q| q > 0.0)
        && condenser_heat < 0.0
        && turbine.balance.shaft_work_j_kg > 0.0
        && final_state.vapor_quality() == Some(0.0))
    {
        return Err("steam train failed heat, work or complete-condensation admission".into());
    }

    println!("# {SOURCE_RELEASE}: {SOURCE_URL}");
    println!(
        "# Partial steady steam train; all energies per kg of throughput; no feed-pump return."
    );
    println!(
        "# Heat positive into water; shaft work positive outward; eta_turbine={TURBINE_EFFICIENCY}."
    );
    println!(
        "state,source,pressure_Pa,temperature_K,enthalpy_J_kg,internal_energy_J_kg,source_entropy_J_kg_K,specific_volume_m3_kg,vapor_quality,pressure_closure_J_kg"
    );
    for (name, state) in [
        ("boiler_inlet", boiler_states[0].into()),
        ("saturated_liquid", boiler_states[1].into()),
        ("saturated_vapor", boiler_states[2].into()),
        ("boiler_outlet", turbine_inlet),
        ("turbine_isentropic_reference", turbine.isentropic_outlet),
        ("turbine_outlet", turbine.outlet),
        ("condenser_outlet", final_state),
    ] {
        row(name, state);
    }
    let stages = [
        (
            "liquid_heating",
            boiler_states[0].into(),
            boiler_states[1].into(),
            heats[0],
            0.0,
        ),
        (
            "boiling",
            boiler_states[1].into(),
            boiler_states[2].into(),
            heats[1],
            0.0,
        ),
        (
            "steam_superheating",
            boiler_states[2].into(),
            turbine_inlet,
            heats[2],
            0.0,
        ),
        (
            "turbine",
            turbine_inlet,
            turbine.outlet,
            0.0,
            turbine.balance.shaft_work_j_kg,
        ),
        (
            "condenser",
            turbine.outlet,
            final_state,
            condenser_heat,
            0.0,
        ),
    ];
    println!(
        "stage,heat_J_kg,shaft_work_J_kg,delta_h_J_kg,delta_u_J_kg,delta_pv_J_kg,delta_source_s_J_kg_K,h_balance_J_kg,u_pv_balance_J_kg"
    );
    let mut maximum_residual = 0.0_f64;
    for (name, before, after, heat, work) in stages {
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
                return Err(format!("{name} energy residual {residual} J/kg").into());
            }
        }
        println!(
            "{name},{heat:.9},{work:.9},{dh:.9},{du:.9},{dpv:.9},{ds:.9},{h_residual:.9e},{u_residual:.9e}"
        );
    }
    let heat_in: f64 = heats.into_iter().sum();
    let work = turbine.balance.shaft_work_j_kg;
    let train_residual = final_state.specific_enthalpy_j_kg() - inlet.specific_enthalpy_j_kg()
        + work
        - heat_in
        - condenser_heat;
    if !train_residual.is_finite() || train_residual.abs() > 5.0 * ENERGY_TOLERANCE_J_KG {
        return Err(format!("steam-train energy residual {train_residual} J/kg").into());
    }
    println!(
        "# boiler_heat={heat_in:.9} J/kg; condenser_heat={condenser_heat:.9} J/kg; outward_shaft_work={work:.9} J/kg"
    );
    println!(
        "# train_energy_residual={train_residual:.9e} J/kg; maximum_stage_residual={maximum_residual:.9e} J/kg"
    );
    println!(
        "# turbine_source_entropy_generation={:.9} J/(kg*K); isentropic_entropy_residual={:.9e} J/(kg*K)",
        turbine.balance.source_entropy_change_j_kg_k, turbine.isentropic_entropy_residual_j_kg_k
    );
    println!(
        "# original_saturation_entropy_defects: boiler={:.9e}; condenser={:.9e} J/(kg*K)",
        boiler.saturation_entropy_defect_j_kg_k(),
        condenser.saturation_entropy_defect_j_kg_k()
    );
    Ok(())
}
