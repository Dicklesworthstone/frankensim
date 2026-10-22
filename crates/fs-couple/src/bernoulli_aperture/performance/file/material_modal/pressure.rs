//! Explicit additional spatial pressure traces and mechanical-sample amplitudes.
use super::{PlateValveInputError, Reader, MAX_CONTROLS, bad, checkpoint};
use crate::bernoulli_aperture::modal::{ModalAperture, pressure::MaterialPressureShape};
use crate::bernoulli_aperture::two_sided::observer::schedule::spatial_pressures::MaterialPressureEvent;
use fs_exec::CancelGate;

pub(super) struct PressureInput {
    blade: usize,
    values: Vec<f64>,
    provenance: String,
    max_projection_terms: usize,
}
impl PressureInput {
    pub(super) fn bind(self, valve: ModalAperture, gate: &CancelGate)
        -> Result<ModalAperture, PlateValveInputError>
    {
        let plate = valve.material_plates().get(self.blade)
            .ok_or_else(|| bad(0, "spatial pressure names an absent original material blade"))?;
        let shape = MaterialPressureShape::from_plate(plate, self.values, self.provenance);
        valve.with_material_pressure_shape(self.blade, shape, self.max_projection_terms, gate)
            .map_err(PlateValveInputError::Physics)
    }
}

pub(super) fn shapes(rows: &mut Reader<'_>, blades: usize, gate: &CancelGate)
    -> Result<Vec<PressureInput>, PlateValveInputError>
{
    if rows.lines.clone().next().and_then(|line| line.split_ascii_whitespace().next())
        != Some("pressure_shapes") { return Ok(Vec::new()); }
    let count = rows.count("pressure_shapes", 2)?;
    let mut inputs = Vec::with_capacity(count);
    for _ in 0..count {
        checkpoint(gate)?;
        let mut row = rows.row("pressure_shape")?;
        let blade = row.parse()?;
        let max_projection_terms = row.count(16_777_216)?;
        let provenance = row.word()?.to_owned();
        if blade >= blades || inputs.iter().any(|p: &PressureInput| p.blade == blade) {
            return Err(bad(row.line, "pressure shape needs an existing uniquely addressed blade"));
        }
        row.finish()?;
        let nodes = rows.count("pressure_nodes", 512)?;
        let mut values = Vec::with_capacity(nodes);
        for _ in 0..nodes {
            checkpoint(gate)?;
            let mut row = rows.row("pressure_node")?;
            values.push(row.scalar()?);
            row.finish()?;
        }
        inputs.push(PressureInput { blade, values, provenance, max_projection_terms });
    }
    Ok(inputs)
}

pub(super) fn events(rows: &mut Reader<'_>, blades: usize, samples: u64, gate: &CancelGate)
    -> Result<[Vec<MaterialPressureEvent>; 2], PlateValveInputError>
{
    let mut events = [Vec::new(), Vec::new()];
    if rows.lines.clone().next().and_then(|line| line.split_ascii_whitespace().next())
        != Some("source_pressures") { return Ok(events); }
    let count = rows.count("source_pressures", MAX_CONTROLS)?;
    for _ in 0..count {
        checkpoint(gate)?;
        let mut row = rows.row("source_pressure")?;
        let sample = row.parse()?;
        let blade: usize = row.parse()?;
        let pressure_pa = row.scalar()?;
        if blade >= blades || blade >= events.len() || sample >= samples {
            return Err(bad(row.line, "spatial pressure needs an existing blade and mechanical sample within the window"));
        }
        row.finish()?;
        events[blade].push(MaterialPressureEvent { sample, pressure_pa });
    }
    Ok(events)
}
