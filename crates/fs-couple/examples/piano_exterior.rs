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
#[path="grand_piano/bridge_response.rs"] mod bridge_response;
#[path="grand_piano/exterior_loading.rs"] mod exterior_loading;
#[path="grand_piano/radiation_fit.rs"] mod radiation_fit;
#[path="grand_piano/exterior_playback.rs"] mod playback;
use exterior_geometry::{Boundary,Specification,RATE};
use std::{fmt::Write as _,io::{Read,Write}};

const USAGE:&str="piano_exterior admittance BOARD.fsb|BOARD.fss SCALE.csv|steinway-d BODY.obj ACOUSTICS.fspe DRIVE_KEY OUTPUT.csv
piano_exterior response BOARD.fsb|BOARD.fss SCALE.csv|steinway-d BODY.obj ACOUSTICS.fspe OUTPUT.csv
piano_exterior render BOARD.fsb|BOARD.fss SCALE.csv|steinway-d BODY.obj ACOUSTICS.fspe OUTPUT.wav SECONDS [PERFORMANCE.mid]
piano_exterior render-loaded BOARD.fsb|BOARD.fss SCALE.csv|steinway-d BODY.obj ACOUSTICS.fspe OUTPUT.wav SECONDS [PERFORMANCE.mid]
    [--modes 1..512] [--substeps 1..16]
    [--hammers materials.fsh] [--hammer-footprints faces.fshp]
    [--dampers estimated|pads.fspd]
    [--performance events.csv | --midi performance.mid]
    [--midi-channel 1..16] [--midi-velocity-max-m-s V] [--midi-half-pedal]
    [--note 21..108] [--velocity m/s]
These playback options apply to both render and render-loaded.
response and admittance also accept --modes/--substeps after the output path,
so harmonic comparisons can use the SAME retained string/board system.
admittance alone accepts --lossless-structure to remove the existing wood and
string material damping for a declared conservative-structure comparison.
It retains the complete complex radiation load and all original modes. Near
fixed-interface string poles use a bounded coupled solve, not fabricated loss.
This option is not accepted by response, render or render-loaded.

Use one explicitly supplied closed outward acoustic skin, including both sides
and edges of a finite soundboard. Label every part as moving or rigid in the
acoustic specification. Separate closed rigid lid/cabinet components participate
in the SAME boundary solve, not an extra source or output EQ. Coordinates and
physical surface motion are three dimensional, including the loaded crown.

The scale keyword steinway-d retains all 88 raw source courses and source felt/
shank mechanics; a CSV preserves its own supplied tensions. Use the SAME scale
that produced any settled/downbearing board. No missing geometry is inferred.
admittance solves a unit peak bridge-force experiment with BEM pressure reacting
on ALL retained string/board coordinates. It writes all bridge mobilities,
receiver Pa/N and wood/string/radiation power balance, alongside an explicit
one-way comparison. This harmonic image excludes hammer/key-damper contacts.
render-loaded fits a passive full-matrix BEM load and couples acoustic storage
into every nonlinear hammer/string/board substep. Acoustic loss is separate from
wood and felt loss. It requires 33..257 odd frequency samples and at most 32
complete board modes; failed passive fits REFUSE, never fall back to one-way.
The ordinary render command remains one-way for a controlled comparison.
response writes complex pressure per mass-normalized modal acceleration in
exp(-i omega t) convention. render fits causal fixed-receiver transfers with
held-out checks, then observes EVERY mechanics substep before PCM encoding.
It uses one score and physical clock for both receivers, no channel normalization.
Default gesture: A4 or the nearest available key at 2 m/s; --note/--velocity
select a supplied-key study. --performance plays the existing sample-accurate CSV,
including jack_staccato/jack_legato force pulses in N and all three pedals.
--midi (or the legacy positional path) uses the existing importer. Channel,
velocity-to-hammer mapping and optional continuous CC64 travel are explicit;
ignored messages and synthetic end releases are reported, never hidden.
CSV, MIDI and demonstration overrides are mutually exclusive. Output
is 48 kHz. Defaults retain four mechanics substeps and at most 24 partials per
string; --modes and --substeps expose the existing larger physical budgets.
They do not retune strings, widen the admitted output band, or certify accuracy.
--hammers supplies complete per-key WoolFelt/Prony cards. --hammer-footprints
selects point or finite longitudinal contact sites for every key. --dampers
selects supplied finite pads or explicit estimates. All use the SAME nonlinear
engine, actual loaded bridge basis and radiation feedback, not output filters.
Invalid or missing supplied cards refuse before any structural/BEM preparation.
See grand_piano/EXTERIOR_PLAYBACK.md for the physical controls.
Transfer accuracy is checked only inside the declared sampled band; an attack
contains out-of-band energy, so this is NOT a full-band realism certificate.
Outputs must be fresh paths; render duration must be 0.05..60 seconds.
See grand_piano/EXTERIOR_ACOUSTICS.md for SI rows and RADIATION_FEEDBACK.md
for the passive-fit, second-order splitting and approximation boundaries.";

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
    let controls=playback::Controls::from_texts(&courses,None,None,None)?;
    prepare_controlled(board_text,courses,obj,spec,&playback::Options::default(),controls)
}
fn prepare_controlled(board_text:&str,courses:Vec<geometry::Course>,obj:&str,spec:Specification,
    options:&playback::Options,controls:playback::Controls)->Result<Scene,String> {
    options.validate()?;
    let keys:Vec<_>=courses.iter().map(|c|c.midi).collect();
    let board=if crowned_board::is_crowned(board_text) {
        crowned_board::CrownedBoard::read(board_text)?.prepare_with_motion(&keys,spec.board_band_hz)?
    } else {board_geometry::BoardGeometry::read(board_text)?.prepare_with_motion(&keys,spec.board_band_hz)?};
    let piano=controls.instrument(courses,&board.modes,options)?;
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
/// Only admittance admits the explicit conservative-structure comparison.
/// Preserve option/value boundaries: a flag cannot repair a missing value.
fn admittance_options(args:&[String])->Result<(playback::Options,bool),String> {
    let mut numeric=Vec::new();let mut damping=true;let mut args=args.iter();
    while let Some(flag)=args.next() {
        if flag=="--lossless-structure" {
            if !damping {return Err("duplicate --lossless-structure".into());}
            damping=false;
        } else {
            numeric.push(flag.clone());
            let value=args.next().filter(|s|!s.starts_with("--"))
                .ok_or_else(||format!("missing value for {flag}"))?;
            numeric.push(value.clone());
        }
    }
    Ok((playback::Options::harmonic(&numeric)?,damping))
}
fn admittance(board_text:&str,courses:&[geometry::Course],obj:&str,spec:&Specification,drive:u8)->Result<String,String> {
    admittance_controlled(board_text,courses,obj,spec,drive,&playback::Options::default(),true)
}
fn admittance_controlled(board_text:&str,courses:&[geometry::Course],obj:&str,spec:&Specification,
    drive:u8,options:&playback::Options,damping:bool)->Result<String,String> {
    options.validate()?;
    let keys:Vec<_>=courses.iter().map(|c|c.midi).collect();
    if !keys.contains(&drive) {return Err("admittance drive key is absent from the scale".into());}
    let board=if crowned_board::is_crowned(board_text) {
        crowned_board::CrownedBoard::read(board_text)?.prepare_with_motion(&keys,spec.board_band_hz)?
    } else {board_geometry::BoardGeometry::read(board_text)?.prepare_with_motion(&keys,spec.board_band_hz)?};
    let model=bridge_response::BridgeResponse::new(courses,&board.modes,RATE*options.substeps as u32,
        0.45*f64::from(RATE),options.modes,damping)?;
    let boundary=Boundary::from_obj(obj,spec,board.motion.as_ref().ok_or("missing harmonic surface motion")?)?
        .loaded(model.bank())?;
    let csv=exterior_loading::sweep(&boundary,&model,spec,drive)?;
    Ok(format!("# structure: {}\n# structural loss: {}; radiation loading remains in the coupled columns\n# board modes={}, retained string coordinates={}, omitted high-frequency duplex mode sets={}\n{}",
        board.provenance,if damping {"original wood/string material damping"}else{"explicitly disabled by --lossless-structure"},
        board.modes.len(),model.bank().modes.len(),model.bank().omitted_duplex_modes,csv))
}
/// Admit both acoustic realizations before attaching the load or dispatching
/// any score event. Loaded and one-way output share the original PCM path.
fn bake(scene:&mut Scene,loaded:bool)->Result<(exterior_audio::Baked,exterior_geometry::Samples,String),String> {
    let (load,samples)=if loaded {
        let (fit,samples)=radiation_fit::prepare(&scene.boundary,&scene.spec,scene.piano.bank.rate)?;
        (Some(fit),samples)
    } else {(None,scene.boundary.sample(&scene.spec)?)};
    let baked=exterior_audio::Baked::from_samples(&samples,scene.spec.fit_order)?;
    let report=if let Some(fit)=load {
        scene.piano.configure_radiation(&fit.model)?;
        format!("Passive load: {} acoustic coordinates; complex matrix peak/RMS={:.6}/{:.6}; resistance peak/RMS={:.6}/{:.6}. These sampled bounds do not certify time-step or spatial convergence.",
            fit.model.poles.len(),fit.peak_error,fit.rms_error,fit.resistance_peak_error,fit.resistance_rms_error)
    } else {String::from("One-way reference: no acoustic reaction on the mechanics.")};
    Ok((baked,samples,report))
}
fn run(args:&[String])->Result<(),String> {
    match args {
        []=>{println!("{USAGE}");Ok(())},
        [help] if help=="--help" || help=="-h"=>{println!("{USAGE}");Ok(())},
        [command,board,strings,obj,spec,drive,output,tail @ ..] if command=="admittance"=>{
            let (options,damping)=admittance_options(tail)?;
            let drive:u8=drive.parse().map_err(|_|"admittance requires a MIDI bridge key in 21..108")?;
            if !(21..=108).contains(&drive) {return Err("admittance drive key outside 21..108".into());}
            if std::path::Path::new(output).exists() {return Err("output must be a fresh path".into());}
            let spec=Specification::read(&read_bounded(spec,exterior_geometry::MAX_SPEC_BYTES)?)?;
            let courses=scale(strings)?;
            if !courses.iter().any(|c|c.midi==drive) {return Err("admittance drive key is absent from the scale".into());}
            let board=read_bounded(board,8*1024*1024)?;
            let obj=read_bounded(obj,exterior_geometry::MAX_OBJ_BYTES)?;
            let csv=admittance_controlled(&board,&courses,&obj,&spec,drive,&options,damping)?;
            publish(output,csv.as_bytes())?;
            println!("Written {output}: radiation-loaded bridge mobility and pressure per 1 N peak, all retained strings and physical loss channels. No time-domain feedback or measured-fidelity claim.");
            Ok(())
        }
        [command,board,strings,obj,spec,output,tail @ ..] if command=="response" || command=="render" || command=="render-loaded"=>{
            let (frames,options)=match command.as_str() {
                "response"=>(None,playback::Options::harmonic(tail)?),
                "render"|"render-loaded" if !tail.is_empty()=>
                    (Some(mesh_render::frames(&tail[0])?),playback::Options::parse(&tail[1..])?),
                _=>return Err(USAGE.into()),
            };
            if std::path::Path::new(output).exists() {return Err("output must be a fresh path".into());}
            let spec=Specification::read(&read_bounded(spec,exterior_geometry::MAX_SPEC_BYTES)?)?;
            if command=="render-loaded" && spec.frequencies<33 {
                return Err("render-loaded requires at least 33 odd-grid frequencies".into());
            }
            let courses=scale(strings)?;let keys:Vec<_>=courses.iter().map(|c|c.midi).collect();
            // Score admission precedes structural/BEM preparation and all writes.
            let score=frames.map(|n|options.score(&keys,n as u64)).transpose()?;
            let controls=playback::Controls::load(&options,&courses)?;
            let geometry=read_bounded(board,8*1024*1024)?;
            let obj=read_bounded(obj,exterior_geometry::MAX_OBJ_BYTES)?;
            let mut scene=prepare_controlled(&geometry,courses,&obj,spec,&options,controls)?;
            if let (Some(n),Some(score))=(frames,score) {
                let (baked,samples,load_report)=bake(&mut scene,command=="render-loaded")?;
                let score_report=score.report;
                let audio=exterior_audio::render(&mut scene.piano,score.performance,n,&baked,scene.spec.full_scale_pa)?;
                publish(output,&audio.wav)?;
                let physical_report=options.report(&scene.piano);
                println!("{}\n{score_report}\n{physical_report}\n{load_report}\nAcoustic source: {}. Structural source: {}.\nBand {:?} Hz; {} panels, {} closed components, minimum panels/wavelength={}, conditioning lower bound={}. Written {output}.",
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
        for key in ["0","109","NaN"] {
            let args=["admittance","missing.fss","missing.csv","missing.obj","missing.fspe",key,"unused.csv"].map(str::to_owned);
            assert!(run(&args).is_err());
        }
    }
    pub(super) fn small_source_inputs()->(String,Vec<geometry::Course>,String,Specification) {
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
        (board,courses,obj,spec)
    }
    pub(super) fn small_source_scene()->Scene {
        let (board,courses,obj,spec)=small_source_inputs();
        prepare(&board,courses,&obj,spec).unwrap()
    }
    #[test]
    fn actual_source_hammer_board_and_exterior_bem_reach_stereo_pcm_on_one_clock() {
        let mut scene=small_source_scene();let mut manual=small_source_scene();
        let samples=scene.boundary.sample(&scene.spec).unwrap();
        let baked=exterior_audio::Baked::from_samples(&samples,scene.spec.fit_order).unwrap();
        let score=||performance::Performance::read("sample,event,key,value\n0,note_on,69,0.5\n1200,note_off,69,0\n",&[69],2400).unwrap();
        let audio=exterior_audio::render(&mut scene.piano,score(),2400,&baked,2.).unwrap();
        assert!(scene.piano.accounting.felt_loss_j>0.,"the source hammer must actually contact the string");
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

    #[test]
    fn supplied_structure_scale_and_body_produce_loaded_bridge_csv() {
        let (board,courses,obj,spec)=small_source_inputs();
        let csv=admittance(&board,&courses,&obj,&spec,69).unwrap();
        let rows:Vec<_>=csv.lines().filter(|l|!l.starts_with('#')).collect();
        assert_eq!(rows.len(),1+spec.frequencies*(courses.len()+spec.receivers.len()));
        assert!(rows[0].contains("radiation_w"));
        let mut changed=false;
        for line in &rows[1..] {
            let cells:Vec<_>=line.split(',').collect();assert_eq!(cells.len(),15);
            let v:Vec<f64>=cells[3..].iter().map(|s|s.parse::<f64>().unwrap()).collect();
            assert!(v.iter().all(|v|v.is_finite()));
            if cells[1]=="bridge" && (v[0]-v[2]).hypot(v[1]-v[3])>1e-8*v[2].hypot(v[3]) {changed=true;}
            assert!(v[7]>=-1e-10); // radiation W: it must not inject mechanical power.
        }
        assert!(changed,"fluid reaction must alter mobility, not merely the printed pressure");
        assert!(admittance(&board,&courses,&obj,&spec,60).is_err());
    }
}

#[cfg(test)]
#[path="grand_piano/radiation_render_tests.rs"]
mod radiation_render_tests;

#[cfg(test)]
#[path="grand_piano/lossless_admittance_tests.rs"]
mod lossless_admittance_tests;
