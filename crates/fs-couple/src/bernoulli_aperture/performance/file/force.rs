//! Optional force records over the same specimen and source sample clock.
use super::*;
use crate::bernoulli_aperture::dynamic::force::PlateForceFootprint;
use crate::bernoulli_aperture::performance::force::{ApertureForceEvent, ApertureForceProgram};

pub(super) fn read(r: &mut Reader<'_>, count: usize, nodes: usize, triangles: usize,
    samples: u64, gate: &CancelGate) -> Result<ApertureForceProgram, PlateValveInputError>
{
    if count == 0 || count > 32 { return Err(bad(r.line, "mechanical program needs 1..=32 force ports")); }
    let mut footprints = Vec::with_capacity(count);
    for _ in 0..count {
        checkpoint(gate)?;
        let mut row = r.row("force_port")?;
        let footprint = match row.word()? {
            "node" => {
                let node: usize = row.parse()?;
                if node >= nodes { return Err(bad(r.line, "force port node is outside the original mesh")); }
                PlateForceFootprint::Node(node)
            }
            "patch" => {
                let n = row.count(triangles)?;
                if n == 0 { return Err(bad(r.line, "force port patch must not be empty")); }
                let mut indices = Vec::with_capacity(n);
                let mut seen = std::collections::BTreeSet::new();
                for _ in 0..n {
                    let index: usize = row.parse()?;
                    if index >= triangles || !seen.insert(index) { return Err(bad(r.line, "force patch repeats or names an absent triangle")); }
                    indices.push(index);
                }
                PlateForceFootprint::Patch(indices)
            }
            _ => return Err(bad(r.line, "force footprint must be a node or explicit triangle patch")),
        };
        row.finish()?; footprints.push(footprint);
    }
    let n = r.count("force_events", MAX_CONTROLS)?;
    let mut events = Vec::with_capacity(n);
    for _ in 0..n {
        checkpoint(gate)?;
        let mut row = r.row("force_event")?;
        let sample: u64 = row.parse()?; let port: usize = row.parse()?; let force_n = row.scalar()?;
        row.finish()?;
        if sample >= samples || port >= count { return Err(bad(r.line, "force event is outside its port list or source window")); }
        events.push(ApertureForceEvent { sample, port, force_n });
    }
    Ok(ApertureForceProgram { footprints, events })
}
