//! Editable OBJ geometry + physical sections -> the real grand_piano plate path.
//!
//! cargo run -p fs-couple --example piano_board_import -- inspect piano.obj
//! cargo run -p fs-couple --example piano_board_import -- import piano.obj board.fspi board.fsb
//! cargo run -p fs-couple --example grand_piano -- --board-geometry board.fsb --scale strings.csv --render piano.wav
//!
//! Offline utility: no new physical solver and no triangle mesh in audio.
#![allow(dead_code)] // Shared piano modules also expose runtime-only operations.
#[path = "grand_piano/geometry.rs"] mod geometry;
#[path = "grand_piano/linear.rs"] mod linear;
#[path = "grand_piano/board.rs"] mod board;
#[path = "grand_piano/board_geometry.rs"] mod board_geometry;
#[path = "grand_piano/steinway_d.rs"] mod steinway_d;
#[path = "grand_piano/steinway_scale.rs"] mod steinway_scale;
#[path = "grand_piano/performance.rs"] mod performance;
#[path = "grand_piano/felt.rs"] mod felt;
#[path = "grand_piano/engine.rs"] mod engine;
#[path = "grand_piano/microphone.rs"] mod microphone;
#[path = "grand_piano/audio.rs"] mod audio;
#[path = "grand_piano/hammer_materials.rs"] mod hammer_materials;
#[path = "grand_piano/mesh_import.rs"] mod mesh_import;
use std::io::{Read, Write};

const USAGE: &str = "piano_board_import inspect INPUT.obj
piano_board_import import INPUT.obj MATERIALS.fspi OUTPUT.fsb
piano_board_import export INPUT.fsb OUTPUT.obj OUTPUT.fspi

inspect reports the actual object/group/material labels without guessing parts.
import selects a flat midsurface, maps every triangle to supplied orthotropic
sections, remaps explicit rib/support vertices and locates each bridge station.
export creates an editable, material-labelled midsurface from an existing native
board, including all its sections, ribs, supports and physical bridge positions.
The native grand_piano --board-geometry path assembles and radiates the result.

Units, frame, flatness tolerance, support choice and physical material constants
must be explicit. Solid case/crowned meshes do not become flat plate models.
Visual MTL values are never elastic constants. No download or MTL file is opened.
Outputs must be fresh paths. See grand_piano/MESH_IMPORT.md for the SI format.";

fn read_bounded(path: &str, cap: usize) -> Result<String, String> {
    let file = std::fs::File::open(path).map_err(|e| format!("{path}: {e}"))?;
    let mut bytes = Vec::new();
    file.take(cap as u64 + 1).read_to_end(&mut bytes).map_err(|e| format!("{path}: {e}"))?;
    if bytes.len() > cap { return Err(format!("{path}: exceeds {cap}-byte input cap")); }
    String::from_utf8(bytes).map_err(|_| format!("{path}: expected UTF-8"))
}
fn create_output(path: &str) -> Result<std::fs::File,String> {
    std::fs::OpenOptions::new().write(true).create_new(true).open(path)
        .map_err(|e| format!("{path}: output must be new and writable: {e}"))
}
fn write_output(mut file:std::fs::File, text:&str, path:&str)->Result<(),String> {
    file.write_all(text.as_bytes()).and_then(|()|file.sync_all()).map_err(|e|format!("{path}: {e}"))
}
fn run(args:&[String])->Result<(),String> {
    match args {
        [] => { println!("{USAGE}"); Ok(()) }
        [flag] if flag == "--help" || flag == "-h" => { println!("{USAGE}"); Ok(()) }
        [command,path] if command == "inspect" => {
            let text=read_bounded(path,mesh_import::MAX_OBJ_BYTES)?;
            let doc=fs_io::obj::read_obj_document(&text).map_err(|e|e.to_string())?;
            println!("{} vertices, {} triangles; SOURCE coordinates, units not inferred",
                doc.soup.positions.len(),doc.soup.triangles.len());
            let mut groups=std::collections::BTreeMap::new();
            for r in &doc.regions {
                *groups.entry((r.object.clone(),r.groups.clone(),r.material.clone())).or_insert(0usize)
                    += r.triangles.len();
            }
            for ((object,groups,material),count) in groups {
                println!("{count} triangles: object={object:?} groups={groups:?} material={material:?}");
            }
            println!("MTL references (not opened): {:?}",doc.material_libraries);
            Ok(())
        }
        [command,obj,spec,output] if command == "import" => {
            let imported=mesh_import::import(&read_bounded(obj,mesh_import::MAX_OBJ_BYTES)?,
                &read_bounded(spec,mesh_import::MAX_SPEC_BYTES)?)?;
            write_output(create_output(output)?,&imported.fsb,output)?;
            println!("Imported {} selected vertices, {} triangles; maximum plane projection {:.9e} m",
                imported.source_vertices.len(),imported.triangles,imported.max_projection_m);
            println!("Written {output}; admitted by the existing piano plate reader. Not a global intersection, measured-geometry or acoustic-fidelity certificate.");
            Ok(())
        }
        [command,input,obj,spec] if command == "export" => {
            if obj == spec { return Err("OBJ and specification outputs must differ".into()); }
            let (geometry,materials)=mesh_import::export(&read_bounded(input,mesh_import::MAX_SPEC_BYTES)?)?;
            // Validate the whole pair before opening either destination. Neither
            // path is overwritten. An OS write/open failure can leave a partial
            // new pair: report the error, never advertise transactional export.
            let a=create_output(obj)?; let b=create_output(spec)?;
            write_output(a,&geometry,obj)?; write_output(b,&materials,spec)?;
            println!("Written {obj} and {spec}; editable physical midsurface and section/beam/bridge sidecar.");
            Ok(())
        }
        _ => Err(USAGE.into()),
    }
}
fn main() {
    if let Err(error)=run(&std::env::args().skip(1).collect::<Vec<_>>()) {
        eprintln!("piano_board_import: {error}"); std::process::exit(1);
    }
}

#[cfg(test)]
mod cli_tests {
    use super::*;
    #[test]
    fn invalid_commands_refuse_and_help_needs_no_files() {
        assert!(run(&["import".into(),"missing.obj".into()]).is_err());
        assert!(run(&["--help".into()]).is_ok());
        assert!(run(&["export".into(),"x.fsb".into(),"same".into(),"same".into()]).is_err());
    }
    #[test]
    fn outputs_do_not_overwrite_measurements() {
        let dir=std::env::temp_dir().join(format!("fs-piano-import-{}-{}",std::process::id(),
            std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos()));
        std::fs::create_dir(&dir).unwrap();
        let path=dir.join("measurement.fsb"); let name=path.to_str().unwrap();
        write_output(create_output(name).unwrap(),"keep these measurements",name).unwrap();
        assert!(create_output(name).is_err());
        assert_eq!(std::fs::read_to_string(&path).unwrap(),"keep these measurements");
        assert!(read_bounded(name,4).is_err());
        // This test owns the freshly created temporary directory, not user files.
        std::fs::remove_file(path).unwrap(); std::fs::remove_dir(dir).unwrap();
    }
}
