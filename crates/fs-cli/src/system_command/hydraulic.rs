//! Hydraulic scene elements owned outside fs-system: the fs-material laminar
//! pipe of a Carreau-Yasuda liquid.

use fs_material::fluid_transport::rheology::pipe::component::CarreauPipeComponent;
use fs_material::fluid_transport::rheology::pipe::{CarreauPipeControls, LaminarCarreauPipe};
use fs_material::fluid_transport::rheology::{
    CarreauYasuda, ThermalCarreauYasuda, ViscosityShiftMode,
};

use super::{J, Refusal, fields, number, positive, scene_shift};

/// A `carreau_pipe`: the authored viscosity curve, its temperature shift
/// and how the shift applies, and the pipe's radius, length and the highest
/// wall shear rate it may run at. The pipe is prepared at the shift's
/// reference temperature; each step reads the temperature at `heat`.
pub(super) fn carreau_pipe(definition: &J, what: &str) -> Result<CarreauPipeComponent, Refusal> {
    fields(
        definition,
        &[
            "kind",
            "viscosity_zero_pa_s",
            "viscosity_infinite_pa_s",
            "time_constant_s",
            "power_index",
            "transition",
            "rate_support_per_s",
            "shift",
            "shift_mode",
            "radius_m",
            "length_m",
            "maximum_wall_rate_per_s",
        ],
        what,
    )?;
    let support = match definition.get("rate_support_per_s").and_then(J::as_array) {
        Some([lowest, highest]) => lowest.as_f64().zip(highest.as_f64()),
        _ => None,
    }
    .ok_or_else(|| format!("{what}: rate_support_per_s is [lowest, highest] in 1/s"))?;
    let curve = CarreauYasuda::new(
        number(definition, "viscosity_zero_pa_s", what)?,
        number(definition, "viscosity_infinite_pa_s", what)?,
        number(definition, "time_constant_s", what)?,
        number(definition, "power_index", what)?,
        number(definition, "transition", what)?,
        [support.0, support.1],
    )
    .map_err(|e| format!("{what}: {e}"))?;
    let shift = scene_shift(definition, what)?.ok_or_else(|| {
        format!(
            "{what}: missing shift, the {{\"arrhenius\": ...}} or {{\"wlf\": ...}} \
             calibration of the viscosity in temperature"
        )
    })?;
    let mode = match definition.get("shift_mode").and_then(J::as_str) {
        Some("vertical") => ViscosityShiftMode::Vertical,
        Some("vertical-and-rate") => ViscosityShiftMode::VerticalAndRate,
        _ => {
            return Err(
                format!("{what}: shift_mode is \"vertical\" or \"vertical-and-rate\"").into(),
            );
        }
    };
    let pipe = LaminarCarreauPipe::new(
        ThermalCarreauYasuda::new(curve, shift, mode),
        shift.reference_temperature_k(),
        positive(definition, "radius_m", what)?,
        positive(definition, "length_m", what)?,
        positive(definition, "maximum_wall_rate_per_s", what)?,
    )
    .map_err(|e| format!("{what}: {e}"))?;
    CarreauPipeComponent::new(pipe, CarreauPipeControls::default())
        .map_err(|e| Refusal::Input(format!("{what}: {e}")))
}
