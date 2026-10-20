//! Explicit acoustic roles on the original mechanical triangles. These records
//! select pressure/volume operators; they never remove material mass or stiffness.
use super::*;
use fs_couple::shear_shell::pressure::ShellPressureSurfaces;
use fs_plate::shell::pressure::{OrientedShellFacet, ShellFace};

pub(super) fn parse(
    rows: &mut Records<'_>,
    triangles: usize,
) -> Result<Option<ShellPressureSurfaces>> {
    if rows
        .lines
        .clone()
        .next()
        .and_then(|line| line.split_ascii_whitespace().next())
        != Some("pressure-surfaces")
    {
        return Ok(None);
    }
    let v = rows.row("pressure-surfaces", 2)?;
    let exterior = count(v[0], 0, triangles)?;
    let cavity = count(v[1], 0, triangles)?;
    // Empty is an explicit zero-pressure role. Absence of the entire block is
    // the only request for the original whole-midsurface compatibility path.
    Ok(Some(ShellPressureSurfaces {
        exterior: read_facets(rows, "exterior-facet", exterior, triangles)?,
        cavity: read_facets(rows, "cavity-facet", cavity, triangles)?,
    }))
}

fn read_facets(
    rows: &mut Records<'_>,
    record: &str,
    number: usize,
    triangles: usize,
) -> Result<Vec<OrientedShellFacet>> {
    let mut result = Vec::<OrientedShellFacet>::with_capacity(number);
    for _ in 0..number {
        let v = rows.row(record, 2)?;
        let triangle = index(v[0], triangles)?;
        if result.iter().any(|facet| facet.triangle == triangle) {
            return Err(
                "shell pressure surface must not repeat a triangle within one fluid role".into(),
            );
        }
        let side = match v[1] {
            "positive" => ShellFace::Positive,
            "negative" => ShellFace::Negative,
            _ => return Err("shell pressure orientation must be positive or negative".into()),
        };
        result.push(OrientedShellFacet { triangle, side });
    }
    Ok(result)
}
