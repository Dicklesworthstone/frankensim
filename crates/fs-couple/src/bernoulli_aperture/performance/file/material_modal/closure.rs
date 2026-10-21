//! Optional physical interface map after the original blade records.
use super::{PlateValveInputError, Reader, bad, checkpoint};
use crate::bernoulli_aperture::plate::closure::PlateClosureSpec;
use fs_exec::CancelGate;

pub(super) fn read(
    rows: &mut Reader<'_>,
    gate: &CancelGate,
) -> Result<Option<PlateClosureSpec>, PlateValveInputError> {
    if rows
        .lines
        .clone()
        .next()
        .and_then(|line| line.split_ascii_whitespace().next())
        != Some("closure")
    {
        return Ok(None);
    }
    let mut row = rows.row("closure")?;
    match row.word()? {
        "none" => {
            row.finish()?;
            return Ok(None);
        }
        "spatial" => {}
        _ => {
            return Err(bad(
                row.line,
                "closure must be none or spatial PRESSURE_STIFFNESS ALPHA LOSS MAX_PENETRATION SOURCE",
            ));
        }
    }
    let stiffness_pa_per_m_alpha = row.scalar()?;
    let alpha = row.scalar()?;
    let internal_loss_s_per_m = row.scalar()?;
    let max_penetration_m = row.scalar()?;
    let provenance = row.word()?.to_owned();
    row.finish()?;
    let count = rows.count("clearances", 512)?;
    let mut nodal_rest_gap_m = Vec::with_capacity(count);
    for _ in 0..count {
        checkpoint(gate)?;
        let mut row = rows.row("gap")?;
        nodal_rest_gap_m.push(row.scalar()?);
        row.finish()?;
    }
    let count = rows.count("lay_triangles", 2048)?;
    let mut lay_triangles = Vec::with_capacity(count);
    for _ in 0..count {
        checkpoint(gate)?;
        let mut row = rows.row("lay_triangle")?;
        lay_triangles.push(row.parse()?);
        row.finish()?;
    }
    Ok(Some(PlateClosureSpec {
        nodal_rest_gap_m,
        lay_triangles,
        stiffness_pa_per_m_alpha,
        alpha,
        internal_loss_s_per_m,
        provenance,
        max_penetration_m,
    }))
}
