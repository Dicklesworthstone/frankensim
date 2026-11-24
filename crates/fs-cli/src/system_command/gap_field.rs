//! Scene admission and reporting for the native moving-gap field owner.
//!
//! The capacitor and inductor share `fs-system`'s `GapField`: one field
//! coordinate and one gap, with its original energy, discrete gradient,
//! tangent, gap limits and checkpoint identity. This boundary only names
//! parameters and coordinates; it evaluates no field or force law.

use fs_ports::component::Component;
use fs_system::elements::GapField;
use fs_system::{NodeId, System, SystemBuilder};

use super::{J, PortSpec, Refusal, add, fields, number, number_or};

pub(super) type State = fs_system::elements::GapFieldState;

/// The native field coordinate's scene spelling.
pub(super) const fn coordinate_name(inductor: bool) -> &'static str {
    if inductor {
        "flux_linkage_wb"
    } else {
        "charge_c"
    }
}

/// Add the original owner with its exact port tangent and admitted state.
pub(super) fn element(
    builder: &mut SystemBuilder,
    name: &str,
    definition: &J,
    inductor: bool,
    what: &str,
) -> Result<(NodeId, Vec<PortSpec>), Refusal> {
    let scale_name = if inductor {
        "mu_turns_area_h_m"
    } else {
        "epsilon_area_f_m"
    };
    let coordinate = coordinate_name(inductor);
    fields(
        definition,
        &["kind", scale_name, "gap_limits_m", "gap_m", coordinate],
        what,
    )?;
    let (lower, upper) = match definition.get("gap_limits_m").and_then(J::as_array) {
        Some([lower, upper]) => lower.as_f64().zip(upper.as_f64()),
        _ => None,
    }
    .ok_or_else(|| format!("{what}.gap_limits_m: expected [lower_m, upper_m]"))?;
    let scale = number(definition, scale_name, what)?;
    let owner = if inductor {
        GapField::inductor(scale, [lower, upper])
    } else {
        GapField::capacitor(scale, [lower, upper])
    }
    .map_err(|error| format!("{what}: {error}"))?;
    let state = State {
        coordinate: number_or(definition, coordinate, 0.0, what)?,
        gap_m: number(definition, "gap_m", what)?,
    };
    let energy = owner
        .stored_energy_j(&state)
        .map_err(|error| format!("{what}: {error}"))?;
    if !energy.is_finite() {
        return Err(format!("{what}: initial stored field energy must be finite").into());
    }
    // The native owner also refuses trial gaps outside this same interval.
    add(builder, name, owner, state)
}

/// Accepted state and the actual work-conjugate port values of the last step.
pub(super) fn values(system: &System, node: NodeId, inductor: bool) -> Vec<(&'static str, f64)> {
    let state = system.state::<State>(node);
    let port_value = |name: &str, output: bool| match system.port(node, name) {
        Some(port) if system.steps() > 0 => {
            let values = if output {
                system.last_output(port)
            } else {
                system.last_input(port)
            };
            values[0]
        }
        _ => f64::NAN,
    };
    vec![
        (
            coordinate_name(inductor),
            state.map_or(f64::NAN, |state| state.coordinate),
        ),
        ("gap_m", state.map_or(f64::NAN, |state| state.gap_m)),
        ("voltage_v", port_value("field", !inductor)),
        ("current_a", port_value("field", inductor)),
        ("attraction_n", port_value("gap", true)),
        ("gap_velocity_m_s", port_value("gap", false)),
    ]
}
