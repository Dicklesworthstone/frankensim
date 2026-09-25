//! Actual paired soundboard skins -> existing grand_piano crowned-shell input.
//! See grand_piano/SOLID_SKINS.md for the explicit correspondence/material law.
#![allow(dead_code)]
#[path = "grand_piano/geometry.rs"] mod geometry;
#[path = "grand_piano/linear.rs"] mod linear;
#[path = "grand_piano/board.rs"] mod board;
#[path = "grand_piano/board_geometry.rs"] mod board_geometry;
#[path = "grand_piano/crowned_board.rs"] mod crowned_board;
#[path = "grand_piano/steinway_d.rs"] mod steinway_d;
#[path = "grand_piano/steinway_scale.rs"] mod steinway_scale;
#[path = "grand_piano/performance.rs"] mod performance;
#[path = "grand_piano/felt.rs"] mod felt;
#[path = "grand_piano/engine.rs"] mod engine;
#[path = "grand_piano/microphone.rs"] mod microphone;
#[path = "grand_piano/audio.rs"] mod audio;
#[path = "grand_piano/hammer_materials.rs"] mod hammer_materials;
#[path = "grand_piano/mesh_import.rs"] mod mesh_import;
#[path = "grand_piano/mesh_render.rs"] mod mesh_render;
use std::io::{Read, Write};

const USAGE: &str = "piano_solid_import INPUT.obj MATERIALS.fspi SKINS.fsps OUTPUT.fss

Select the UPPER skin with the existing FSPI part row. SKINS.fsps declares the
lower part and thickness,geometry. Supply complete one-based vertex pairs or
pairing,projected,TOLERANCE_IN_METRES for uniquely registered XY skin vertices.
Actual midpoint heights and mean normal separations replace a hand-authored
midsurface and nominal panel thickness. Material tensors, grain, supports,
bridges and beams remain explicit. No MTL parameters become physical constants.
Both skins must have matching triangle connectivity and pass shallow-shell
admission. No remeshing, welding, side-wall guessing or measured-fidelity claim.
Inspect labels first with piano_board_import inspect INPUT.obj.
Output must be a new path. See grand_piano/SOLID_SKINS.md.";

fn read(path: &str, cap: usize) -> Result<String, String> {
    let file = std::fs::File::open(path).map_err(|e| format!("{path}: {e}"))?;
    let mut bytes = Vec::new();
    file.take(cap as u64 + 1).read_to_end(&mut bytes).map_err(|e| format!("{path}: {e}"))?;
    if bytes.len() > cap { return Err(format!("{path}: exceeds {cap}-byte input cap")); }
    String::from_utf8(bytes).map_err(|_| format!("{path}: expected UTF-8"))
}
fn run(args: &[String]) -> Result<(), String> {
    if args.is_empty() || (args.len() == 1 && matches!(args[0].as_str(), "--help" | "-h")) {
        println!("{USAGE}"); return Ok(());
    }
    let [obj, spec, skins, output] = args else { return Err(USAGE.into()); };
    if std::path::Path::new(output).exists() { return Err(format!("{output}: output must be new")); }
    let imported = mesh_import::crowned::solid::import(
        &read(obj, mesh_import::MAX_OBJ_BYTES)?,
        &read(spec, mesh_import::MAX_SPEC_BYTES)?,
        &read(skins, mesh_import::MAX_SPEC_BYTES)?)?;
    let mut file = std::fs::OpenOptions::new().write(true).create_new(true).open(output)
        .map_err(|e| format!("{output}: cannot create new output: {e}"))?;
    file.write_all(imported.fsb.as_bytes()).and_then(|()| file.sync_all())
        .map_err(|e| format!("{output}: {e}"))?;
    println!("Written {output}: {} midsurface nodes, {} geometry-derived sections, maximum |midpoint height| {:.9e} m.\nUse grand_piano --preset steinway-d --board-geometry {output} --render piano.wav.\nThis imports supplied geometry; it does not certify a measured Steinway or converged acoustic fidelity.",
        imported.source_vertices.len(), imported.triangles, imported.max_projection_m);
    Ok(())
}
fn main() {
    if let Err(e) = run(&std::env::args().skip(1).collect::<Vec<_>>()) {
        eprintln!("piano_solid_import: {e}"); std::process::exit(1);
    }
}
#[cfg(test)]
mod cli_tests {
    use super::*;
    #[test]
    fn help_and_invalid_arity_do_not_open_assets() {
        assert!(run(&["--help".into()]).is_ok());
        assert!(run(&["missing.obj".into()]).is_err());
        assert!(run(&["a".into(),"b".into(),"c".into()]).is_err());
    }
}
