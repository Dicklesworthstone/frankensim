//! The authored finite inventory and power history, recorded before playback.
use fs_couple::bernoulli_aperture::performance::AperturePerformance;
use fs_couple::dg_port::finite_gas::FiniteGasDgLoss;

pub(super) fn provenance(renderer: &AperturePerformance, clock_hz: u32) -> String {
    let Some(program) = renderer.finite_gas_heat_program() else { return String::new(); };
    let gas = renderer.finite_gas().expect("heat program retains its original inventory");
    let state = gas.state().state();
    let controls = program.controls.iter().map(|control| format!(
        "{{\"sample\":{},\"power_w\":{:e}}}", control.sample, control.power_w))
        .collect::<Vec<_>>().join(",");
    let (loss, error) = match program.coupling.loss {
        FiniteGasDgLoss::Lossless => ("lossless", String::new()),
        FiniteGasDgLoss::UniformViscothermal { maximum_relative_error } =>
            ("uniform-viscothermal", format!(",\"maximum_relative_error\":{maximum_relative_error:e}")),
        FiniteGasDgLoss::SpatialViscothermal { maximum_relative_error } =>
            ("spatial-viscothermal", format!(",\"maximum_relative_error\":{maximum_relative_error:e}")),
        FiniteGasDgLoss::UniformMolecularDiffusion { maximum_relative_error } =>
            ("molecular-diffusion", format!(",\"maximum_relative_error\":{maximum_relative_error:e}")),
    };
    let caloric = program.caloric;
    let coupling = program.coupling;
    format!(",\"finite_gas_heat\":{{\"clock_hz\":{clock_hz},\"phase\":\"frozen\",\"initial_mass_kg\":{:e},\"volume_m3\":{:e},\"initial_temperature_k\":{:e},\"initial_pressure_pa\":{:e},\"controls\":[{controls}],\"caloric\":{{\"temperature_bracket_k\":[{:e},{:e}],\"energy_tolerance_j_kg\":{:e},\"maximum_iterations\":{}}},\"coupling\":{{\"work_bracket_j\":[{:e},{:e}],\"work_tolerance_j\":{:e},\"energy_tolerance_j\":{:e},\"maximum_iterations\":{},\"maximum_temperature_increment_k\":{:e},\"loss\":\"{loss}\"{error}}},\"scope\":\"signed external power held on the original physical clock; fixed-volume closed material inventory pays acoustic parameter work; all gas snapshots describe the initial state; no mean mass transport or computed exterior pressure jump\"}}",
        gas.mass_kg(), gas.volume_m3(), state.temperature, state.pressure,
        caloric.minimum_temperature_k, caloric.maximum_temperature_k,
        caloric.energy_tolerance_j_kg, caloric.maximum_iterations,
        coupling.work_bracket_j[0], coupling.work_bracket_j[1], coupling.work_tolerance_j,
        coupling.energy_tolerance_j, coupling.maximum_iterations,
        coupling.maximum_temperature_increment_k)
}
