//! Native section-skin selection and OBJ export for the existing exterior CLI.
use super::{Boundary, Specification, board_geometry, crowned_board};
use board_geometry::motion::{MotionSurface, skin::Skin};

pub const BODY: &str = "board-skin";
pub const CONTINUOUS_BODY: &str = "board-skin-continuous";

/// Selection is typed after path parsing: an invalid OBJ file whose CONTENT is
/// "board-skin" remains invalid OBJ and cannot trigger generated fallback.
pub fn read_body(path: &str) -> Result<(Option<String>, bool), String> {
    if path == BODY { Ok((None, false)) }
    else if path == CONTINUOUS_BODY { Ok((None, true)) }
    else { super::read_bounded(path, super::exterior_geometry::MAX_OBJ_BYTES).map(|s| (Some(s), false)) }
}
pub fn boundary(obj: Option<&str>, source: &str, spec: &Specification, motion: &MotionSurface, continuous: bool)
    -> Result<(Boundary, String), String> {
    if obj.is_some() && continuous { return Err("a supplied OBJ cannot also select thickness reconstruction".into()); }
    match obj {
        None => Boundary::from_board_skin(source, spec, motion, continuous),
        Some(text) => Ok((Boundary::from_obj(text, spec, motion)?, String::from("supplied OBJ acoustic boundary"))),
    }
}

/// Export the actual equilibrium skin, not the unloaded input dressed as a
/// loaded body. Structural admission/equilibrium/eigenanalysis remain with the
/// existing owner. Export invokes NO BEM or acoustic fit and claims neither.
pub fn export(board_path: &str, scale_path: &str, spec_path: &str, output: &str, continuous: bool) -> Result<(), String> {
    if std::path::Path::new(output).exists() { return Err("skin export output must be a fresh path".into()); }
    let spec = Specification::read(&super::read_bounded(spec_path, super::exterior_geometry::MAX_SPEC_BYTES)?)?;
    spec.require_board_skin()?;
    let courses = super::scale(scale_path)?;
    let keys: Vec<_> = courses.iter().map(|c| c.midi).collect();
    let source = super::read_bounded(board_path, 8 * 1024 * 1024)?;
    let prepared = if crowned_board::is_crowned(&source) {
        crowned_board::CrownedBoard::read(&source)?.prepare_with_motion(&keys, spec.board_band_hz)?
    } else { board_geometry::BoardGeometry::read(&source)?.prepare_with_motion(&keys, spec.board_band_hz)? };
    let motion = prepared.motion.as_ref().ok_or("skin export lacks the complete prepared motion geometry")?;
    // Export permits a larger inspectable asset than the online BEM wrapper's
    // 2048 panels. No modes/triangles are silently removed to meet either bound.
    let skin = if continuous { Skin::continuous_from_source(motion, &source, spec.offset_m, 250_000)? }
        else { Skin::from_source(motion, &source, spec.offset_m, 250_000)? };
    let image = if continuous { "explicit volume-preserving continuous thickness" } else { "exact piecewise-section columns" };
    let text = format!("# Skin image: {image}\n# Structural source: {}\n{}", prepared.provenance.replace(['\r','\n'], " "), skin.obj());
    super::publish(output, text.as_bytes())?;
    println!("Written {output} ({image}): {} vertices, {} panels, section/skin volume {:.12e}/{:.12e} m3, maximum facet mean-thickness change {:.3e} m. Actual prepared equilibrium; no acoustic solve/fit or scanned-cabinet claim. BEM playback retains its separate 2048-panel limit.",
        skin.vertices.len(), skin.triangles.len(), skin.section_volume_m3, skin.volume_m3, skin.maximum_thickness_change_m);
    Ok(())
}

#[cfg(test)]
#[path="section_skin_render_tests.rs"]
mod tests;
