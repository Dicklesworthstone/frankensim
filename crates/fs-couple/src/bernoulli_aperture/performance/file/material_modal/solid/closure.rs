//! Original solid exterior faces/edges, with no reconstructed or fitted shape.
use super::{PlateValveInputError, Reader, SolidImage, bad, checkpoint};
use crate::bernoulli_aperture::modal::solid::SolidSurfaceClosureSpec;
use crate::bernoulli_aperture::plate::closure::PlateClosureSpec;
use fs_exec::CancelGate;

pub(super) fn read(
    rows: &mut Reader<'_>,
    image: &SolidImage,
    gate: &CancelGate,
) -> Result<Option<SolidSurfaceClosureSpec>, PlateValveInputError> {
    if rows
        .lines
        .clone()
        .next()
        .and_then(|line| line.split_ascii_whitespace().next())
        != Some("solid_closure")
    {
        return Ok(None);
    }
    let mut row = rows.row("solid_closure")?;
    match row.word()? {
        "none" => {
            row.finish()?;
            return Ok(None);
        }
        "spatial" => {}
        _ => return Err(bad(row.line, "solid_closure requires none or spatial")),
    }
    let opening_sign = row.scalar()?;
    let stiffness_pa_per_m_alpha = row.scalar()?;
    let alpha = row.scalar()?;
    let internal_loss_s_per_m = row.scalar()?;
    let max_penetration_m = row.scalar()?;
    let max_projection_terms = row.count(16_777_216)?;
    let provenance = row.word()?.to_owned();
    row.finish()?;
    let surface = match image {
        SolidImage::Frozen(source) => source.surface_geometry(),
        SolidImage::Hereditary { source, .. } => source.surface_geometry(),
    }
    .ok_or_else(|| bad(rows.line, "solid closure requires its original P1 exterior"))?;
    let count = rows.count("solid_clearances", image.node_count())?;
    if count != image.node_count() {
        return Err(bad(
            rows.line,
            "solid clearances must cover every original volume node",
        ));
    }
    let mut nodal_rest_gap_m = Vec::with_capacity(count);
    for _ in 0..count {
        checkpoint(gate)?;
        let mut row = rows.row("solid_gap")?;
        nodal_rest_gap_m.push(row.scalar()?);
        row.finish()?;
    }
    let count = rows.count("solid_lay_faces", surface.faces().len())?;
    let mut lay_triangles = Vec::with_capacity(count);
    for _ in 0..count {
        checkpoint(gate)?;
        let mut row = rows.row("solid_lay_face")?;
        let mut nodes: [usize; 3] = [row.parse()?, row.parse()?, row.parse()?];
        nodes.sort_unstable();
        let id = surface
            .faces()
            .binary_search_by_key(&nodes, |face| face.nodes)
            .map_err(|_| bad(row.line, "solid lay must name an original exterior face"))?;
        lay_triangles.push(id);
        row.finish()?;
    }
    let count = rows.count("solid_slit_edges", 4096)?;
    let mut slit_edges = Vec::with_capacity(count);
    for _ in 0..count {
        checkpoint(gate)?;
        let mut row = rows.row("solid_slit_edge")?;
        slit_edges.push([row.parse()?, row.parse()?]);
        row.finish()?;
    }
    Ok(Some(SolidSurfaceClosureSpec {
        contact: PlateClosureSpec {
            nodal_rest_gap_m,
            lay_triangles,
            stiffness_pa_per_m_alpha,
            alpha,
            internal_loss_s_per_m,
            provenance,
            max_penetration_m,
        },
        slit_edges,
        opening_sign,
        max_projection_terms,
    }))
}
