//! Actual geometric piano motion -> supplied finite acoustic body -> receiver Pa.
//! No sampled piano, infinite baffle, synthesized cabinet, or second mechanics.
#![allow(dead_code)] // Shared piano modules expose additional front doors.
#[path="grand_piano/geometry.rs"] mod geometry;
#[path="grand_piano/linear.rs"] mod linear;
#[path="grand_piano/board.rs"] mod board;
#[path="grand_piano/board_geometry.rs"] mod board_geometry;
#[path="grand_piano/crowned_board.rs"] mod crowned_board;
#[path="grand_piano/steinway_d.rs"] mod steinway_d;
#[path="grand_piano/steinway_scale.rs"] mod steinway_scale;
#[path="grand_piano/performance.rs"] mod performance;
#[path="grand_piano/felt.rs"] mod felt;
#[path="grand_piano/engine.rs"] mod engine;
#[path="grand_piano/microphone.rs"] mod microphone;
#[path="grand_piano/audio.rs"] mod audio;
#[path="grand_piano/hammer_materials.rs"] mod hammer_materials;
#[path="grand_piano/mesh_import.rs"] mod mesh_import;
#[path="grand_piano/mesh_render.rs"] mod mesh_render;
#[path="grand_piano/exterior_geometry.rs"] mod exterior_geometry;
#[path="grand_piano/exterior_audio.rs"] mod exterior_audio;
use exterior_geometry::{Boundary,Specification,RATE};
use std::{fmt::Write as _,io::{Read,Write}};

const USAGE:&str="piano_exterior response BOARD.fsb|BOARD.fss SCALE.csv|steinway-d BODY.obj ACOUSTICS.fspe OUTPUT.csv
piano_exterior render BOARD.fsb|BOARD.fss SCALE.csv|steinway-d BODY.obj ACOUSTICS.fspe OUTPUT.wav SECONDS [PERFORMANCE.mid]

Use one explicitly supplied closed outward acoustic skin, including both sides
and edges of a finite soundboard. Label every part as moving or rigid in the
acoustic specification. Separate closed rigid lid/cabinet components participate
in the SAME boundary solve, not an extra source or output EQ. Coordinates and
physical surface motion are three dimensional, including the loaded crown.

The scale keyword steinway-d retains all 88 raw source courses and source felt/
shank mechanics; a CSV preserves its own supplied tensions. Use the SAME scale
that produced any settled/downbearing board. No missing geometry is inferred.
response writes complex pressure per mass-normalized modal acceleration in
exp(-i omega t) convention. render fits causal fixed-receiver transfers with
held-out checks, then observes EVERY mechanics substep before PCM encoding.
It uses one score and physical clock for both receivers, no channel normalization.
Default gesture: A4 at 2 m/s; optional MIDI uses the existing importer. Output
is 48 kHz with four mechanics substeps and at most 24 partials per string.
Transfer accuracy is checked only inside the declared sampled band; an attack
contains out-of-band energy, so this is NOT a full-band realism certificate.
Outputs must be fresh paths; render duration must be 0.05..60 seconds.
See grand_piano/EXTERIOR_ACOUSTICS.md for SI rows and approximation boundaries.";

fn read_bounded(path:&str,cap:usize)->Result<String,String> {
    let file=std::fs::File::open(path).map_err(|e|format!("{path}: {e}"))?;
    let mut bytes=Vec::new();file.take(cap as u64+1).read_to_end(&mut bytes).map_err(|e|e.to_string())?;
    if bytes.len()>cap {return Err(format!("{path}: input exceeds {cap} bytes"));}
    String::from_utf8(bytes).map_err(|_|format!("{path}: input must be UTF-8"))
}
fn publish(path:&str,bytes:&[u8])->Result<(),String> {
    let mut file=std::fs::OpenOptions::new().write(true).create_new(true).open(path)
        .map_err(|e|format!("{path}: fresh writable output required: {e}"))?;
    file.write_all(bytes).and_then(|()|file.sync_all()).map_err(|e|e.to_string())
}
fn scale(input:&str)->Result<Vec<geometry::Course>,String> {
    if input=="steinway-d" {steinway_scale::courses()}
    else {geometry::read_scale(&read_bounded(input,512*1024)?)}
}
struct Scene {
    piano:engine::Instrument,
    board:board_geometry::PreparedBoard,
    boundary:Boundary,
    spec:Specification,
}
fn prepare(board_text:&str,courses:Vec<geometry::Course>,obj:&str,spec:Specification)->Result<Scene,String> {
    let keys:Vec<_>=courses.iter().map(|c|c.midi).collect();
    let board=if crowned_board::is_crowned(board_text) {
        crowned_board::CrownedBoard::read(board_text)?.prepare_with_motion(&keys,spec.board_band_hz)?
    } else {board_geometry::BoardGeometry::read(board_text)?.prepare_with_motion(&keys,spec.board_band_hz)?};
    let materials=courses.iter().map(steinway_scale::hammer_material).collect::<Result<Vec<_>,_>>()?;
    let piano=engine::Instrument::new_with_course_shanks(courses,&board.modes,RATE,4,24,true,
        materials,engine::ShankGeometry::published())?;
    let boundary=Boundary::from_obj(obj,&spec,board.motion.as_ref().ok_or("missing full-vector structural motion")?)?
        .loaded(&piano.bank)?;
    Ok(Scene {piano,board,boundary,spec})
}
fn response(scene:&Scene)->Result<String,String> {
    let samples=scene.boundary.sample(&scene.spec)?;
    let mut csv=format!("# finite exterior BEM, exp(-i omega t), Pa per unit mass-normalized modal acceleration\n# source: {}\n# structure: {}\n# panels={}, components={}, min_panels_per_wavelength={}, condition_lower_bound_max={}\n# rigid scatterers, one-way acoustics; no radiation loading, flexible cabinet, room or above-band claim\nfrequency_hz,receiver,input,real_pa_per_acceleration,imag_pa_per_acceleration\n",
        scene.spec.source,scene.board.provenance,scene.boundary.surface.areas().len(),scene.boundary.components,
        samples.minimum_ppw,samples.maximum_condition_lower_bound);
    for (f,w) in samples.omega.iter().enumerate() {for (receiver,inputs) in samples.values.iter().enumerate() {
        for (input,row) in inputs.iter().enumerate() {
            writeln!(csv,"{:.17e},{receiver},{input},{:.17e},{:.17e}",w/std::f64::consts::TAU,row[f].re,row[f].im).unwrap();
        }
    }}
    Ok(csv)
}
fn run(args:&[String])->Result<(),String> {
    match args {
        []=>{println!("{USAGE}");Ok(())},
        [help] if help=="--help" || help=="-h"=>{println!("{USAGE}");Ok(())},
        [command,board,strings,obj,spec,output,tail @ ..] if command=="response" || command=="render"=>{
            let frames=match command.as_str() {
                "response" if tail.is_empty()=>None,
                "render" if (1..=2).contains(&tail.len())=>Some(mesh_render::frames(&tail[0])?),
                _=>return Err(USAGE.into()),
            };
            if std::path::Path::new(output).exists() {return Err("output must be a fresh path".into());}
            let spec=Specification::read(&read_bounded(spec,exterior_geometry::MAX_SPEC_BYTES)?)?;
            let courses=scale(strings)?;let keys:Vec<_>=courses.iter().map(|c|c.midi).collect();
            // Score admission precedes structural/BEM preparation and all writes.
            let score=frames.map(|n|match tail.get(1) {
                Some(path)=>{
                    let midi=performance::midi::load(path,&keys,RATE,n as u64,performance::midi::Mapping::default())?;
                    performance::Performance::from_events(midi.events,&keys,n as u64)
                }
                None=>performance::Performance::demonstration(&keys,RATE,n as u64,Some(69),Some(2.)),
            }).transpose()?;
            let geometry=read_bounded(board,8*1024*1024)?;
            let obj=read_bounded(obj,exterior_geometry::MAX_OBJ_BYTES)?;
            let mut scene=prepare(&geometry,courses,&obj,spec)?;
            if let (Some(n),Some(score))=(frames,score) {
                let samples=scene.boundary.sample(&scene.spec)?;
                let baked=exterior_audio::Baked::from_samples(&samples,scene.spec.fit_order)?;
                let audio=exterior_audio::render(&mut scene.piano,score,n,&baked,scene.spec.full_scale_pa)?;
                publish(output,&audio.wav)?;
                println!("{}\nAcoustic source: {}. Structural source: {}.\nBand {:?} Hz; {} panels, {} closed components, minimum panels/wavelength={}, conditioning lower bound={}. Written {output}.",
                    audio.report,scene.spec.source,scene.board.provenance,scene.spec.band_hz,
                    scene.boundary.surface.areas().len(),scene.boundary.components,samples.minimum_ppw,samples.maximum_condition_lower_bound);
            } else {
                let csv=response(&scene)?;publish(output,csv.as_bytes())?;
                println!("Written {output}: finite-body pressure from actual structural modes and supplied acoustic geometry; no measured-fidelity claim.");
            }
            Ok(())
        }
        _=>Err(USAGE.into()),
    }
}
fn main() {
    if let Err(error)=run(&std::env::args().skip(1).collect::<Vec<_>>()) {
        eprintln!("piano_exterior: {error}");std::process::exit(1);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn malformed_front_doors_do_not_open_or_overwrite_files() {
        assert!(run(&["--help".into()]).is_ok());
        assert!(run(&["response".into(),"missing.fss".into()]).is_err());
        assert!(Specification::read("frankensim-piano-exterior-si-v1\n").is_err());
    }
    fn small_source_scene()->Scene {
        let mut board=String::from("frankensim-board-geometry-si-v1\nsource,estimated,soft-panel acoustic integration NOT Steinway geometry\nsupport,clamped\npretension,0\ndamping,0.01\n");
        for (i,p) in [[0.,0.],[0.1,0.],[0.1,0.1],[0.,0.1],[0.05,0.05]].iter().enumerate() {
            writeln!(board,"node,{i},{},{}",p[0],p[1]).unwrap();
        }
        for (i,t) in [[0,1,4],[1,2,4],[2,3,4],[3,0,4]].iter().enumerate() {
            writeln!(board,"triangle,{i},{},{},{},0.003,450,1e7,8e5,0.3,6e5,0",t[0],t[1],t[2]).unwrap();
        }
        board.push_str("fixed,0\nfixed,1\nfixed,2\nfixed,3\nbridge,69,0,0,0,1\n");
        let spec=exterior_geometry::tests::specification()
            .replace("band-hz,40,400,17","band-hz,40,300,41")
            .replace("board-band-hz,400","board-band-hz,300");
        // Duplicate receivers must share the same mechanics, yet own histories.
        let spec=Specification::read(&format!("{spec}receiver-m,0.05,0.05,1\n")).unwrap();
        let obj=exterior_geometry::tests::box_obj("skin",[0.,0.,-0.0015],[0.1,0.1,0.003]);
        let courses=steinway_scale::courses().unwrap().into_iter().filter(|c|c.midi==69).collect();
        prepare(&board,courses,&obj,spec).unwrap()
    }
    #[test]
    fn actual_source_hammer_board_and_exterior_bem_reach_stereo_pcm_on_one_clock() {
        let mut scene=small_source_scene();let mut manual=small_source_scene();
        let samples=scene.boundary.sample(&scene.spec).unwrap();
        let baked=exterior_audio::Baked::from_samples(&samples,scene.spec.fit_order).unwrap();
        let score=||performance::Performance::read("sample,event,key,value\n0,note_on,69,0.1\n1200,note_off,69,0\n",&[69],2400).unwrap();
        let audio=exterior_audio::render(&mut scene.piano,score(),2400,&baked,2.).unwrap();
        assert!(audio.peak_pa>1e-14);assert_eq!(&audio.wav[..4],b"RIFF");
        assert_eq!(u16::from_le_bytes([audio.wav[22],audio.wav[23]]),2);
        // No separately advanced left/right piano or observation backreaction.
        let mut program=score();
        for n in 0..2400 {program.dispatch(n,&mut manual.piano).unwrap();manual.piano.step().unwrap();}
        assert_eq!(scene.piano.bank.q,manual.piano.bank.q);assert_eq!(scene.piano.bank.v,manual.piano.bank.v);
        assert_eq!(scene.piano.accounting.input_work_j,manual.piano.accounting.input_work_j);
        let residual=scene.piano.accounting.input_work_j-scene.piano.energy_j()-scene.piano.accounting.dissipated_j();
        assert!(residual.abs()<1e-7);
        let data=audio.wav.windows(4).position(|w|w==b"data").unwrap()+8;
        for frame in audio.wav[data..].chunks_exact(4) {assert_eq!(&frame[..2],&frame[2..]);}
    }

}
