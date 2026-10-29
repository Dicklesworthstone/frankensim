//! Heat and cool dense water at fixed mass and volume with original IF97 Region 3.
//!
//! Run: `cargo run -p fs-material --example water_dense_vessel`.
//! A 1-litre rigid volume contains 0.5 kg of water at 650 K. Prescribed heat
//! raises its temperature to approximately 750 K over 100 seconds, then the
//! same power is reversed for 100 seconds. Density remains 500 kg/m3, while
//! pressure follows the original density/temperature Helmholtz equation.
//!
//! This is an illustrative, spatially uniform quasistatic equilibrium path.
//! It supplies neither a container-strength prediction nor wall compliance,
//! wall heat storage, finite-rate transport or a heater-reservoir entropy model.
//! The positive Cv bound and source-domain admission cover the entire prepared
//! interval. Pressure enthalpy V Delta p is reported separately from the exactly
//! zero work at the rigid boundary.

use fs_material::phase::water::dense::{DenseIsochore, DenseWaterControls, DenseWaterState};
use fs_qty::{Density, Temperature};
use fs_thermochem::water::if97::region3::increments::DenseIsochoricIntervalSolve;
use fs_thermochem::water::if97::{SOURCE_RELEASE, SOURCE_URL};

const MASS_KG: f64 = 0.5;
const VOLUME_M3: f64 = 0.001;
const DT_S: f64 = 5.0;
const HEATING_STEPS: u32 = 20;
const STEP_ENERGY_TOLERANCE_J: f64 = 1e-5;
const HISTORY_ENERGY_TOLERANCE_J: f64 = 5e-4;

fn row(time: f64, state: DenseWaterState, heat: f64, initial: DenseWaterState) {
    let du =
        MASS_KG * (state.specific_internal_energy_j_kg() - initial.specific_internal_energy_j_kg());
    let dh = MASS_KG * (state.specific_enthalpy_j_kg() - initial.specific_enthalpy_j_kg());
    let vdp = VOLUME_M3 * (state.pressure_pa() - initial.pressure_pa());
    let ds = MASS_KG * (state.specific_entropy_j_kg_k() - initial.specific_entropy_j_kg_k());
    println!(
        "{time:.1},{:.9},{:.9},{:.6},{heat:.9},{du:.9},{dh:.9},{vdp:.9},{ds:.12}",
        state.temperature_k(),
        state.pressure_pa() * 1e-6,
        state.density_kg_m3(),
    );
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let law = DenseIsochore::new(
        Density::new(MASS_KG / VOLUME_M3),
        DenseWaterControls {
            // The authored heat/cool trajectory is 650..750 K. Preparing a
            // slightly wider valid interval leaves its source endpoints free
            // from accumulated floating-point inverse residuals.
            minimum_temperature_k: 649.0,
            maximum_temperature_k: 751.0,
            internal_energy_tolerance_j_kg: 1e-6,
            maximum_iterations: 96,
            interval: DenseIsochoricIntervalSolve {
                maximum_intervals: 65_535,
            },
        },
        || false,
    )?;
    let initial = law.state_at_temperature(Temperature::new(650.0), || false)?;
    let hot = law.state_at_temperature(Temperature::new(750.0), || false)?;
    let heating_j =
        MASS_KG * (hot.specific_internal_energy_j_kg() - initial.specific_internal_energy_j_kg());
    let per_step_j = heating_j / f64::from(HEATING_STEPS);
    let power_w = per_step_j / DT_S;
    let interval = law.source_interval();
    println!("# {SOURCE_RELEASE}: {SOURCE_URL}");
    println!(
        "# fixed_mass={MASS_KG} kg; fixed_volume={VOLUME_M3} m3; density={} kg/m3; heater={power_w:.9} W",
        law.density().value(),
    );
    println!(
        "# prepared_temperature_interval={:?} K; proved_Cv=[{:.9},{:.9}] J/(kg K); interval_evaluations={}",
        interval.temperature_bounds_k(),
        interval.minimum_heat_capacity_j_kg_k(),
        interval.maximum_heat_capacity_j_kg_k(),
        interval.interval_evaluations(),
    );
    println!("# Q=Delta U; Delta H=Delta U+V Delta p; rigid-wall work=0 J.");
    println!(
        "time_s,temperature_K,pressure_MPa,density_kg_m3,cumulative_heat_J,delta_U_J,delta_H_J,V_delta_p_J,source_entropy_change_J_K"
    );
    row(0.0, initial, 0.0, initial);
    let mut current = initial;
    let mut net_heat = 0.0;
    let mut stable_delta_u = 0.0;
    let mut stable_delta_h = 0.0;
    let mut stable_vdp = 0.0;
    let mut source_entropy = 0.0;
    let mut maximum_step_residual = 0.0_f64;
    let mut maximum_history_residual = 0.0_f64;
    let mut inverse_iterations = 0;
    for index in 1..=2 * HEATING_STEPS {
        let heat = if index <= HEATING_STEPS {
            per_step_j
        } else {
            -per_step_j
        };
        let change = law.advance_heat(current, MASS_KG, heat, STEP_ENERGY_TOLERANCE_J, || false)?;
        if change.boundary_work_j != 0.0 || change.volume_m3 != VOLUME_M3 {
            return Err("dense-water rigid volume or boundary work changed".into());
        }
        for residual in [
            change.energy_residual_j,
            change.enthalpy_residual_j,
            change.represented_energy_residual_j,
            change.represented_enthalpy_residual_j,
        ] {
            maximum_step_residual = maximum_step_residual.max(residual.abs());
        }
        net_heat += heat;
        stable_delta_u += change.internal_energy_change_j;
        stable_delta_h += change.enthalpy_change_j;
        stable_vdp += change.pressure_enthalpy_change_j;
        source_entropy += change.entropy_change_j_k;
        inverse_iterations += change.inversion_iterations;
        current = change.state;
        let represented_du = MASS_KG
            * (current.specific_internal_energy_j_kg() - initial.specific_internal_energy_j_kg());
        let represented_dh =
            MASS_KG * (current.specific_enthalpy_j_kg() - initial.specific_enthalpy_j_kg());
        let represented_vdp = VOLUME_M3 * (current.pressure_pa() - initial.pressure_pa());
        for residual in [
            stable_delta_u - net_heat,
            stable_delta_h - stable_vdp - net_heat,
            represented_du - net_heat,
            represented_dh - represented_vdp - net_heat,
        ] {
            maximum_history_residual = maximum_history_residual.max(residual.abs());
            if !residual.is_finite() || residual.abs() > HISTORY_ENERGY_TOLERANCE_J {
                return Err(format!("dense-water history energy residual {residual} J").into());
            }
        }
        if index % 2 == 0 {
            row(f64::from(index) * DT_S, current, net_heat, initial);
        }
    }
    if net_heat.abs() > HISTORY_ENERGY_TOLERANCE_J
        || (current.temperature_k() - initial.temperature_k()).abs() > 1e-6
        || (current.pressure_pa() - initial.pressure_pa()).abs() > 1.0
    {
        return Err(
            "the dense-water heat/cool cycle did not recover its initial source state".into(),
        );
    }
    println!(
        "# heating={heating_j:.9} J; hot_pressure={:.9} MPa; hot_V_delta_p={:.9} J",
        hot.pressure_pa() * 1e-6,
        VOLUME_M3 * (hot.pressure_pa() - initial.pressure_pa()),
    );
    println!(
        "# maximum_step_energy_residual={maximum_step_residual:.9e} J; maximum_history_energy_residual={maximum_history_residual:.9e} J; cycle_net_heat={net_heat:.9e} J"
    );
    println!(
        "# source_entropy_cycle_change={source_entropy:.9e} J/K; total_caloric_inverse_iterations={inverse_iterations}"
    );
    println!(
        "# The source interval is reused across all heat steps; reservoir entropy is unspecified."
    );
    Ok(())
}
