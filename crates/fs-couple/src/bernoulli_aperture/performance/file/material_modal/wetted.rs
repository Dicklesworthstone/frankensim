//! Original physical faces in the existing pressure/structure/flow junction.
use super::{PlateValveInputError, Reader, bad, checkpoint};
use crate::bernoulli_aperture::modal::wetted::MaterialWettedInput;
use fs_exec::CancelGate;
use fs_plate::shell::pressure::{OrientedShellFacet, ShellFace};

pub(super) fn read(rows: &mut Reader<'_>, blades: usize, gate: &CancelGate)
    -> Result<Vec<MaterialWettedInput>, PlateValveInputError>
{
    if rows.lines.clone().next().and_then(|line| line.split_ascii_whitespace().next())
        != Some("wetted_surfaces") { return Ok(Vec::new()); }
    let count = rows.count("wetted_surfaces", 2)?;
    let mut inputs = Vec::with_capacity(count);
    for _ in 0..count {
        checkpoint(gate)?;
        let mut row = rows.row("wetted_surface")?;
        let blade = row.parse()?;
        let max_projection_terms = row.count(16_777_216)?;
        let provenance = row.word()?.to_owned();
        if blade >= blades || inputs.iter().any(|input: &MaterialWettedInput| input.blade == blade) {
            return Err(bad(row.line, "wetted surface needs an existing uniquely addressed blade"));
        }
        row.finish()?;
        let count = rows.count("wetted_facets", 2048)?;
        let mut facets = Vec::with_capacity(count);
        for _ in 0..count {
            checkpoint(gate)?;
            let mut row = rows.row("wetted_facet")?;
            let triangle = row.parse()?;
            let side = match row.word()? {
                "positive" => ShellFace::Positive,
                "negative" => ShellFace::Negative,
                _ => return Err(bad(row.line, "upstream physical face must be positive or negative relative to original winding")),
            };
            row.finish()?;
            facets.push(OrientedShellFacet { triangle, side });
        }
        inputs.push(MaterialWettedInput { blade, facets, provenance, max_projection_terms });
    }
    Ok(inputs)
}
