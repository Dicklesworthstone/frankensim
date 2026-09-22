//! Editable OBJ geometry + physical sections -> the real grand_piano plate/shell path.
//!
//! cargo run -p fs-couple --example piano_board_import -- inspect piano.obj
//! cargo run -p fs-couple --example piano_board_import -- import piano.obj board.fspi board.fsb
//! cargo run -p fs-couple --example piano_board_import -- import-crowned piano.obj board.fspi board.fss
//!
//! Offline utility: no new physical solver and no triangle mesh in audio.
#![allow(dead_code)] // Shared piano modules also expose runtime-only operations.
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

const USAGE: &str = "piano_board_import inspect INPUT.obj
piano_board_import import INPUT.obj MATERIALS.fspi OUTPUT.fsb
piano_board_import import-crowned INPUT.obj MATERIALS.fspi OUTPUT.fss
piano_board_import export INPUT.fsb OUTPUT.obj OUTPUT.fspi
piano_board_import render-steinway INPUT.fsb|INPUT.fss OUTPUT.wav SECONDS [PERFORMANCE.mid]

inspect reports the actual object/group/material labels without guessing parts.
import selects a flat midsurface, maps every triangle to supplied orthotropic
sections, remaps explicit rib/support vertices and locates each bridge station.
import-crowned retains the selected OBJ's source heights as a 3-D CST/DKT shell.
It requires clamped supports, zero prestress and a shallow graph within 50 mm
of the declared reference plane. Rectangular beam sections are reconstructed
from their supplied area and bending inertia; actual offsets retain coupling.
This is initial crown geometry, NOT a solved downbearing/prestress equilibrium.
export creates an editable, material-labelled midsurface from a FLAT native
board, including all its sections, ribs, supports and physical bridge positions.
render-steinway accepts either native board type with source Model D strings,
felt cards and shanks, then the existing 48 kHz physical-pressure/PCM path.
Optional MIDI uses channel 1 and velocity 127 -> 4.5 m/s; otherwise strike key
69 at 2 m/s. Raw source tensions are preserved. Retention: 400 Hz board, 24
string partials. Crowned radiation projects signed normal volume velocity to
the flat baffle; this is NOT full 3-D radiation or room/lid scattering.
The native grand_piano --board-geometry path remains flat-only for custom cards.

Units, frame, support choice and physical material constants must be explicit.
Solid cabinet meshes do not become soundboard midsurfaces. Visual MTL values
are never elastic constants. No download or MTL file is opened.
Outputs must be fresh paths. See grand_piano/MESH_IMPORT.md and CROWNED_BOARD.md.";

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
        [command,obj,spec,output] if command == "import" || command == "import-crowned" => {
            let obj=read_bounded(obj,mesh_import::MAX_OBJ_BYTES)?;
            let spec=read_bounded(spec,mesh_import::MAX_SPEC_BYTES)?;
            let crowned=command=="import-crowned";
            let imported=if crowned {mesh_import::crowned::import(&obj,&spec)?}
                else {mesh_import::import(&obj,&spec)?};
            write_output(create_output(output)?,&imported.fsb,output)?;
            println!("Imported {} selected vertices, {} triangles; {} {:.9e} m",
                imported.source_vertices.len(),imported.triangles,
                if crowned {"retained maximum |crown height|"} else {"maximum plane projection"},
                imported.max_projection_m);
            println!("Written {output}; {}. Not a global intersection, measured-geometry or acoustic-fidelity certificate.",
                if crowned {"3-D shell with supplied crown, orthotropy and eccentric reinforcement"}
                else {"admitted by the existing piano plate reader"});
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
        [command,input,output,seconds] if command == "render-steinway" => {
            render_file(input,output,seconds,None)
        }
        [command,input,output,seconds,midi] if command == "render-steinway" => {
            render_file(input,output,seconds,Some(midi))
        }
        _ => Err(USAGE.into()),
    }
}
fn render_file(input:&str, output:&str, seconds:&str, midi:Option<&str>)->Result<(),String> {
    let count=mesh_render::frames(seconds)?;
    if std::path::Path::new(output).exists() {
        return Err(format!("{output}: render output must be a fresh path"));
    }
    let text=read_bounded(input,mesh_import::MAX_SPEC_BYTES)?;
    let rendered=mesh_render::render(&text,count,midi)?;
    // Race-safe no-clobber creation after successful preparation and rendering.
    let mut file=create_output(output)?;
    file.write_all(&rendered.wav).and_then(|()|file.sync_all()).map_err(|e|format!("{output}: {e}"))?;
    println!("{}",rendered.report);
    Ok(())
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
        assert!(run(&["import-crowned".into(),"missing.obj".into()]).is_err());
        assert!(run(&["--help".into()]).is_ok());
        assert!(run(&["render-steinway".into(),"missing.fsb".into(),"out.wav".into(),"NaN".into()]).is_err());
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
