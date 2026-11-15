//! `modal_structure`: a linear structure in mass-normalized damped modes,
//! `fs-nlmodal`'s [`ModalStructure`], acted on through its attachments.
//!
//! The scene authors the modes (frequency and viscous damping ratio, and a
//! start in modal coordinates) and, for each attachment, one row of mode-shape
//! weights per port row. Every mode advances exactly under the effort it is
//! held at over a step; the work, the stored modal energy and the viscous heat
//! are the owner's. The structure publishes the end of each attachment's flow,
//! so a load line may sit beside it. Nothing here integrates anything.

use core::num::NonZeroUsize;

use fs_nlmodal::component::{
    ModalAcousticState, ModalAcousticTimeBudget, ModalAttachment, ModalStructure,
};
use fs_ports::StableId;
use fs_ports::component::{Carrier, PortRole, PortSpec};
use fs_system::{NodeId, System, SystemBuilder};

use super::{Domain, J, Refusal, add, fields, number, number_or};

/// What the kernel holds for a modal structure: one state per mode.
pub(super) type State = Vec<ModalAcousticState>;

/// A `modal_structure` element: the owner's model admitted at the scene's
/// clock, its attachments, and the start of every mode.
pub(super) fn structure(
    builder: &mut SystemBuilder,
    name: &str,
    definition: &J,
    what: &str,
) -> Result<(NodeId, Vec<PortSpec>), Refusal> {
    fields(
        definition,
        &["kind", "sample_rate_hz", "modes", "attachments"],
        what,
    )?;
    let rate = number(definition, "sample_rate_hz", what)?;
    if !(1.0..=f64::from(u32::MAX)).contains(&rate) || rate.fract() != 0.0 {
        return Err(format!("{what}: sample_rate_hz is a whole number of hertz").into());
    }
    let rows = definition
        .get("modes")
        .and_then(J::as_array)
        .filter(|rows| !rows.is_empty())
        .ok_or_else(|| format!("{what}: modes is a nonempty array of modes"))?;
    let mut modes = Vec::with_capacity(rows.len());
    let mut start = Vec::with_capacity(rows.len());
    for (index, row) in rows.iter().enumerate() {
        let inner = format!("{what}.modes[{index}]");
        fields(
            row,
            &["frequency_hz", "damping_ratio", "displacement", "velocity"],
            &inner,
        )?;
        let frequency_hz = number(row, "frequency_hz", &inner)?;
        modes.push((
            core::f64::consts::TAU * frequency_hz,
            number(row, "damping_ratio", &inner)?,
        ));
        start.push(ModalAcousticState {
            displacement_m_sqrt_kg: number_or(row, "displacement", 0.0, &inner)?,
            velocity_m_sqrt_kg_per_s: number_or(row, "velocity", 0.0, &inner)?,
        });
    }
    let table = definition
        .get("attachments")
        .and_then(J::as_object)
        .filter(|table| !table.is_empty())
        .ok_or_else(|| format!("{what}: attachments is an object of at least one attachment"))?;
    let mut attachments = Vec::with_capacity(table.len());
    for (id, attachment) in table {
        let inner = format!("{what}.attachments.{id}");
        fields(attachment, &["domain", "weights"], &inner)?;
        let domain = attachment
            .str_field("domain")
            .and_then(Domain::parse)
            .ok_or_else(|| {
                format!("{inner}: domain is one of translation, rotation, electrical")
            })?;
        let weight_rows = attachment
            .get("weights")
            .and_then(J::as_array)
            .unwrap_or_default();
        let width = NonZeroUsize::new(weight_rows.len()).ok_or_else(|| {
            format!("{inner}: weights is an array of rows, one weight per mode in each")
        })?;
        let mut weights = Vec::with_capacity(width.get() * modes.len());
        for (row, entries) in weight_rows.iter().enumerate() {
            let entries = entries
                .as_array()
                .filter(|entries| entries.len() == modes.len())
                .ok_or_else(|| {
                    format!(
                        "{inner}.weights[{row}]: one weight per mode, {} here",
                        modes.len()
                    )
                })?;
            for entry in entries {
                weights.push(
                    entry
                        .as_f64()
                        .filter(|weight| weight.is_finite())
                        .ok_or_else(|| format!("{inner}.weights[{row}]: finite numbers"))?,
                );
            }
        }
        attachments.push(ModalAttachment {
            id: StableId::new(id.as_str())
                .map_err(|_| format!("{inner}: {id:?} is not a port identifier"))?,
            kind: domain.port_kind(),
            width,
            weights,
        });
    }
    // The whole number checked above fits a u32.
    #[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
    let rate = rate as u32;
    let structure = ModalStructure::from_modes(
        rate,
        &modes,
        ModalAcousticTimeBudget::audible_reference(),
        attachments,
    )
    .and_then(ModalStructure::publishing_ends)
    .map_err(|error| format!("{what}: {error}"))?;
    add(builder, name, structure, start)
}

/// The report of a modal structure: every mode's coordinates, the flow at the
/// end of the last step at every row of every attachment, and the last step's
/// heat rate.
pub(super) fn values(system: &System, node: NodeId, ports: &[PortSpec]) -> Vec<(String, f64)> {
    let mut values = Vec::new();
    if let Some(state) = system.state::<State>(node) {
        for (mode, state) in state.iter().enumerate() {
            values.push((
                format!("mode{mode}.displacement"),
                state.displacement_m_sqrt_kg,
            ));
            values.push((
                format!("mode{mode}.velocity"),
                state.velocity_m_sqrt_kg_per_s,
            ));
        }
    }
    let last = |name: &str| {
        system
            .port(node, name)
            .filter(|_| system.steps() > 0)
            .map(|port| system.last_output(port))
    };
    for port in ports {
        let Carrier::Power(kind) = port.carrier() else {
            continue;
        };
        let (Some(domain), PortRole::EffortIn) = (Domain::of(kind), port.role()) else {
            continue;
        };
        let id = port.id().as_str();
        let ends = last(&format!("{id}.end"));
        for row in 0..port.width() {
            let key = if port.width() == 1 {
                format!("{id}.{}", domain.names().flow)
            } else {
                format!("{id}{row}.{}", domain.names().flow)
            };
            let end = ends
                .and_then(|ends| ends.get(row))
                .copied()
                .unwrap_or(f64::NAN);
            values.push((key, end));
        }
    }
    let heat_w = last("heat")
        .and_then(|values| values.first())
        .map_or(f64::NAN, |value| -value);
    values.push(("heat_w".to_owned(), heat_w));
    values
}
