//! Heat and cool a sealed 5-litre water vessel with evolving pressure and phase.
//!
//! Run: `cargo run -p fs-material --example water_rigid_vessel`.
//! The fixed 50 g water inventory starts in liquid/vapor equilibrium at 300 K.
//! Prescribed heat raises it to approximately 600 K and is then removed. The
//! original IF97 U,V equilibrium inverse determines pressure, temperature and
//! vapor fraction after every accepted heat interval. No mass leaves the vessel.
//!
//! This is quasistatic spatially uniform equilibrium, with a rigid, massless
//! wall and an explicitly prescribed heater/cooler. Wall strength/compliance,
//! finite phase-change rates, wall heat storage and reservoir entropy are not
//! modeled. In particular, Delta H contains V Delta p, but boundary work is zero.

use fs_material::phase::water::isochoric::{
    IsochoricWater, IsochoricWaterControls, IsochoricWaterState,
};
use fs_qty::Temperature;
use fs_thermochem::water::MassSpecificVolume;
use fs_thermochem::water::if97::{SOURCE_RELEASE, SOURCE_URL};

const MASS_KG: f64 = 0.05;
const VOLUME_M3: f64 = 0.005;
const DT_S: f64 = 5.0;
const HEATING_STEPS: u32 = 20;
const STEP_ENERGY_TOLERANCE_J: f64 = 1e-7;
const CYCLE_ENERGY_TOLERANCE_J: f64 = 1e-5;

fn row(time: f64, state: IsochoricWaterState, heat: f64, initial: IsochoricWaterState) {
    let du =
        MASS_KG * (state.specific_internal_energy_j_kg() - initial.specific_internal_energy_j_kg());
    let dh = MASS_KG * (state.specific_enthalpy_j_kg() - initial.specific_enthalpy_j_kg());
    let vdp = VOLUME_M3 * (state.pressure_pa() - initial.pressure_pa());
    println!(
        "{time:.1},{:?},{:.9},{:.6},{:.9},{:.9},{:.12},{heat:.9},{du:.9},{dh:.9},{vdp:.9}",
        state.phase(),
        state.temperature_k(),
        state.pressure_pa(),
        state.vapor_mass_fraction(),
        MASS_KG * state.vapor_mass_fraction(),
        state.volume_residual_m3_kg() * MASS_KG,
    );
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let law = IsochoricWater::new(
        MassSpecificVolume::new(VOLUME_M3 / MASS_KG),
        IsochoricWaterControls {
            minimum_temperature_k: 295.0,
            maximum_temperature_k: 623.15,
            minimum_pressure_pa: 100.0,
            maximum_pressure_pa: 100e6,
            volume_tolerance_m3_kg: 1e-14,
            internal_energy_tolerance_j_kg: 1e-6,
            maximum_pressure_iterations: 96,
            maximum_temperature_iterations: 96,
        },
        || false,
    )?;
    let initial = law.state_at_temperature(Temperature::new(300.0), || false)?;
    let hot = law.state_at_temperature(Temperature::new(600.0), || false)?;
    let heating_j =
        MASS_KG * (hot.specific_internal_energy_j_kg() - initial.specific_internal_energy_j_kg());
    let per_step_j = heating_j / f64::from(HEATING_STEPS);
    let power_w = per_step_j / DT_S;
    println!("# {SOURCE_RELEASE}: {SOURCE_URL}");
    println!(
        "# fixed_mass={MASS_KG} kg; fixed_volume={VOLUME_M3} m3; heater={power_w:.9} W; heating_time={} s",
        f64::from(HEATING_STEPS) * DT_S
    );
    println!("# Q=Delta U; Delta H=Delta U+V Delta p; boundary work=0 J.");
    println!(
        "time_s,phase,temperature_K,pressure_Pa,vapor_mass_fraction,vapor_kg,source_volume_residual_m3,cumulative_heat_J,delta_U_J,delta_H_J,V_delta_p_J"
    );
    row(0.0, initial, 0.0, initial);
    let mut current = initial;
    let mut net_heat = 0.0;
    let mut maximum_step_residual = 0.0_f64;
    let mut maximum_cycle_residual = 0.0_f64;
    let mut source_entropy_change = 0.0;
    for index in 1..=2 * HEATING_STEPS {
        let heat = if index <= HEATING_STEPS {
            per_step_j
        } else {
            -per_step_j
        };
        let change = law.advance_heat(current, MASS_KG, heat, STEP_ENERGY_TOLERANCE_J, || false)?;
        if change.boundary_work_j != 0.0 || change.fixed_volume_m3 != VOLUME_M3 {
            return Err("rigid-vessel geometry or boundary work changed".into());
        }
        for residual in [
            change.internal_energy_residual_j,
            change.enthalpy_residual_j,
        ] {
            maximum_step_residual = maximum_step_residual.max(residual.abs());
        }
        net_heat += heat;
        source_entropy_change += change.source_entropy_change_j_k;
        current = change.state;
        // Independently reconstruct from the complete actual endpoint states.
        let du = MASS_KG
            * (current.specific_internal_energy_j_kg() - initial.specific_internal_energy_j_kg());
        let dh = MASS_KG * (current.specific_enthalpy_j_kg() - initial.specific_enthalpy_j_kg());
        let vdp = VOLUME_M3 * (current.pressure_pa() - initial.pressure_pa());
        for residual in [du - net_heat, dh - vdp - net_heat] {
            maximum_cycle_residual = maximum_cycle_residual.max(residual.abs());
            if !residual.is_finite() || residual.abs() > CYCLE_ENERGY_TOLERANCE_J {
                return Err(
                    format!("rigid-vessel complete-history energy residual {residual} J").into(),
                );
            }
        }
        if index % 2 == 0 {
            row(f64::from(index) * DT_S, current, net_heat, initial);
        }
    }
    if net_heat.abs() > CYCLE_ENERGY_TOLERANCE_J
        || (current.temperature_k() - initial.temperature_k()).abs() > 1e-6
        || (current.pressure_pa() - initial.pressure_pa()).abs() > 1e-3
        || (current.vapor_mass_fraction() - initial.vapor_mass_fraction()).abs() > 1e-10
    {
        return Err("the heat/cool cycle did not recover its initial equilibrium state".into());
    }
    println!(
        "# heating={heating_j:.9} J; hot_pressure={:.9} Pa; hot_vapor_fraction={:.9}",
        hot.pressure_pa(),
        hot.vapor_mass_fraction()
    );
    println!(
        "# maximum_step_energy_residual={maximum_step_residual:.9e} J; maximum_history_energy_residual={maximum_cycle_residual:.9e} J; cycle_net_heat={net_heat:.9e} J"
    );
    println!(
        "# raw_source_entropy_cycle_change={source_entropy_change:.9e} J/K; initial_saturation_gibbs_gap={:.9e} J/kg",
        initial.saturation_gibbs_gap_j_kg().unwrap_or(0.0)
    );
    println!(
        "# Entropy is the raw IF97 value; source matching errors are not heat or entropy production."
    );
    Ok(())
}
