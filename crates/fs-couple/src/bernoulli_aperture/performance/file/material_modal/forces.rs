//! Mechanical-time assignments on the original source-defined plate ports.
use super::{PlateValveInputError, Reader, MAX_CONTROLS, bad, checkpoint};
use crate::bernoulli_aperture::performance::force::ApertureForceEvent;
use fs_exec::CancelGate;

pub(super) fn read(rows: &mut Reader<'_>, blades: usize, samples: u64, gate: &CancelGate)
    -> Result<[Vec<ApertureForceEvent>; 2], PlateValveInputError>
{
    let mut events = [Vec::new(), Vec::new()];
    if rows.lines.clone().next().and_then(|line| line.split_ascii_whitespace().next())
        != Some("source_forces") {
        return Ok(events);
    }
    let mut row = rows.row("source_forces")?;
    let count = match row.word()? {
        "none" => { row.finish()?; return Ok(events); }
        count => count.parse::<usize>()
            .map_err(|_| bad(row.line, "source_forces needs a bounded count or none"))?,
    };
    if count == 0 || count > MAX_CONTROLS {
        return Err(bad(row.line, "source_forces needs 1..=262144 assignments shared by both blades"));
    }
    row.finish()?;
    for _ in 0..count {
        checkpoint(gate)?;
        let mut row = rows.row("source_force")?;
        let sample = row.parse()?;
        let blade: usize = row.parse()?;
        let force_n = row.scalar()?;
        if blade >= blades || blade >= events.len() || sample >= samples {
            return Err(bad(row.line, "source force needs an existing blade and mechanical sample within the source window"));
        }
        row.finish()?;
        // Each original material plate already owns one complete force footprint.
        // Do not collapse its multimode projection into an opening coordinate.
        events[blade].push(ApertureForceEvent { sample, port: 0, force_n });
    }
    Ok(events)
}
