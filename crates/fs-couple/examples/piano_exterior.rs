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
#[path="grand_piano/string_polarization.rs"] mod string_polarization;
#[path="grand_piano/cavity.rs"] mod cavity;
#[path="grand_piano/mesh_import.rs"] mod mesh_import;
#[path="grand_piano/mesh_render.rs"] mod mesh_render;
#[path="grand_piano/exterior_geometry.rs"] mod exterior_geometry;
#[path="grand_piano/exterior_audio.rs"] mod exterior_audio;
#[path="grand_piano/bridge_response.rs"] mod bridge_response;
#[path="grand_piano/exterior_loading.rs"] mod exterior_loading;
#[path="grand_piano/radiation_fit.rs"] mod radiation_fit;
#[path="grand_piano/exterior_playback.rs"] mod playback;
#[path="grand_piano/section_skin_cli.rs"] mod section_skin;
use exterior_geometry::{Boundary,Specification,RATE,rigid::Assembly};
use std::{fmt::Write as _,io::{Read,Write}};

const USAGE:&str="piano_exterior export-skin BOARD.fsb|BOARD.fss SCALE.csv|steinway-d ACOUSTICS.fspe OUTPUT.obj [--continuous-thickness]
piano_exterior export-rigid ASSEMBLY.fspr OUTPUT.obj
piano_exterior admittance BOARD.fsb|BOARD.fss SCALE.csv|steinway-d BODY.obj ACOUSTICS.fspe DRIVE_KEY OUTPUT.csv
piano_exterior response BOARD.fsb|BOARD.fss SCALE.csv|steinway-d BODY.obj ACOUSTICS.fspe OUTPUT.csv
piano_exterior render BOARD.fsb|BOARD.fss SCALE.csv|steinway-d BODY.obj ACOUSTICS.fspe OUTPUT.wav SECONDS [PERFORMANCE.mid]
piano_exterior render-loaded BOARD.fsb|BOARD.fss SCALE.csv|steinway-d BODY.obj ACOUSTICS.fspe OUTPUT.wav SECONDS [PERFORMANCE.mid]
    [--modes 1..512] [--substeps 1..16] [--rigid-assembly ASSEMBLY.fspr]
    [--equilibrate-board-mass] [--consistent-board-mass | --edge-cubic-board-mass]
    [--board-reduction max_modes,keep_low_modes,Hz,...]
    [--hammers materials.fsh] [--hammer-footprints faces.fshp]
    [--rt0425-hammer-stiffness] [--rt0425-hammer-dissipation]
    [--rt0425-string-damping]
    [--dampers estimated|pads.fspd]
    [--string-polarization bridge-frames.fspp]
    [--cavity enclosure.fspc]
    [--performance events.csv | --midi performance.mid]
    [--midi-channel 1..16] [--midi-velocity-max-m-s V] [--midi-half-pedal]
    [--note 21..108] [--velocity m/s]
These playback options apply to both render and render-loaded.
--string-polarization uses both transverse directions of each physical string,
projected from the same full-vector board modes. Every key must supply its
bridge site, 3-D arm, orthonormal string/hammer frame and lateral damper ratio.
Missing motion or a primary projection inconsistent with the board refuses;
no lateral coupling is guessed. See grand_piano/STRING_POLARIZATION.md.
--cavity adds reciprocal sealed-cavity compression and standing waves below
the flat board, using supplied dimensions, air state and momentum damping.
It composes with both render and render-loaded, including stereo receivers;
see grand_piano/CAVITY.md. Supply a closed enclosure OBJ with the exposed board
face moving and enclosure walls rigid. Bare board-skin selections and moving
underside panels refuse: that face is already inside the cavity. The enclosing
volume must agree with the cavity card; it is not inferred from the exterior OBJ.
The RT-0425 hammer flags require the steinway-d scale and no supplied hammer
cards. Dissipation also requires RT-0425 stiffness. They select the already
implemented per-string K_H and published R_H contact laws, not output EQ.
The RT-0425 string flag requires the steinway-d scale and projects published
per-key R_u and eta_u onto the existing reduced string modes; it is opt-in.
response and admittance also accept --modes/--substeps, --rigid-assembly,
--equilibrate-board-mass, either flat-board mass option, --board-reduction,
--string-polarization and --rt0425-string-damping after the output path. Admittance retains both
transverse string directions and the selected intrinsic loss law from playback.
Its applied bridge force and reported bridge velocity remain in the primary
hammer direction; secondary strings react through the same board. The response
command is pressure per supplied modal acceleration, so intrinsic damping does
not turn that acoustic transfer into a force-driven structural response.
Admittance also accepts --cavity, preserving the same geometry-derived air
compression, standing-wave inertia and explicit momentum loss as playback.
All acoustic coordinates remain in the coupled solve, including at lossless
fixed-wall resonances and coincident string resonances. Selected sweeps append
a cavity_w power column; one_way columns still retain the selected cavity and
omit only exterior radiation reaction. The prescribed-motion response command
rejects --cavity because it has no mechanical force/reaction solve.
Mass equilibration is an opt-in numerical solve for a flat geometric board;
it leaves geometry, materials, mode cap and original residual admission unchanged.
--consistent-board-mass selects consistent P1 panel inertia. --edge-cubic-board-mass
selects the existing cubic transverse field for inertia, bridge coupling and
acoustic surface motion, including its analytic physical rotations. Both retain
the same eigensolve throughout playback and can compose with string polarization.
--board-reduction admits up to 512 certified source modes within the supplied
FSPE board-band-hz, then retains at most 128 coordinates using bridge
static/harmonic response directions and the requested exact low modes. Both
directions of supplied --string-polarization frames guide the basis. Its
1..16 increasing target frequencies must lie inside that source band. Full
projected wood damping and the transformed bridge/acoustic motion are shared
by playback and admittance. Set max_modes <=32 for render-loaded, whose passive
radiation fit keeps its existing independent budget. Snapshot displacement
projection error is not a transfer, acoustic or mesh-convergence certificate.
See grand_piano/BOARD_REDUCTION.md for the explicit format and scope.
Reduction also accepts crowned shells and uses the equilibrium tangent modes
when downbearing is supplied. Mass options remain flat-only; higher-band
convergence remains an independent requirement. The optional FSB stiffener-mass row accepts lumped,
consistent-hermite or consistent-eccentric. The last adds supplied bending rotary
and offset centroid inertia to Hermite translation; no torsional polar inertia is inferred.
admittance alone accepts --lossless-structure to remove the existing wood and
string material damping for a declared conservative-structure comparison.
It retains the complete complex radiation load and all selected coordinates. Near
fixed-interface string poles use a bounded coupled solve, not fabricated loss.
This option is not accepted by response, render or render-loaded.
It also excludes --rt0425-string-damping, which explicitly selects a lossy law.

BODY.obj can instead be the explicit keyword board-skin: derive both faces,
rim/hole walls and section-thickness steps from the SAME prepared physical board.
Both board-skin and board-skin-continuous require SI identity transform and only
moving,soundboard_skin. The continuous variant EXPLICITLY reconstructs a nodal
thickness field while preserving total section volume; per-facet changes are
reported. It is suitable for a smooth taper represented by cellwise sections.
Neither variant invents a cabinet/lid or changes the structural cards.
export-skin (optionally --continuous-thickness) writes that equilibrium surface as OBJ,
without a BEM solve; it can exceed the renderer's separate 2048-panel limit.
The 1 nm height grid and volume discrepancy are reported. See SECTION_SKIN.md.

--rigid-assembly appends selected, posed OBJ lid/cabinet parts to EITHER native
board-skin variant or a supplied BODY.obj. Exact source labels, source units,
SI hinge axes/pivots and translations are explicit. Rigid parts join the SAME
BEM with zero source velocity; pressure AND radiation reaction reflect the pose.
Assets load once before modal preparation. No auto-remeshing, material guesses,
missing-file fallback, dynamic hinges or extra board mass. The combined scene
still obeys the 2048-panel budget. export-rigid writes the posed SI parts for
inspection without a board eigensolve or BEM. See RIGID_ASSEMBLY.md.

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
exp(-i omega t) convention. response and admittance admit 1..64 receivers,
sharing the same source solve; render and render-loaded require one or two.
Larger audio requests refuse before structural or acoustic preparation.
render fits causal fixed-receiver transfers with
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
selects point, uniform span or authored crown profile for every key. Profile
sites supply offsets, recession, local thickness and fractions of the original
area; their gaps control engagement and their felt/Prony histories are independent.
Published R_H requires the original uniform felt thickness. See HAMMER_FOOTPRINTS.md.
--dampers selects supplied finite pads or explicit estimates. All use the SAME nonlinear
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
    prepare_controlled_body(board_text,courses,Some(obj),spec,options,controls,false)
}
/// One structural solve owns bridge, string-frame and acoustic motion fields.
fn prepare_board_motion(board_text:&str,keys:&[u8],band_hz:f64,options:&playback::Options)
    ->Result<board_geometry::PreparedBoard,String> {
    prepare_board_motion_with_source_ports(board_text,keys,band_hz,options,None)
}
fn prepare_board_motion_with_source_ports(board_text:&str,keys:&[u8],band_hz:f64,
    options:&playback::Options,frames:Option<&[board_geometry::motion::SourceBridgeFrame]>)
    ->Result<board_geometry::PreparedBoard,String> {
    options.validate()?;
    if crowned_board::is_crowned(board_text) {
        if options.equilibrate_board_mass || options.consistent_board_mass || options.edge_cubic_board_mass {
            return Err("flat-board mass controls require a flat geometric board".into());
        }
        let geometry=crowned_board::CrownedBoard::read(board_text)?;
        if let Some(reduction)=&options.board_reduction {
            geometry.prepare_reduced_with_ports(keys,band_hz,true,reduction,frames)
        } else {geometry.prepare_with_motion(keys,band_hz)}
    } else {
        let geometry=board_geometry::BoardGeometry::read(board_text)?;
        if let Some(reduction)=&options.board_reduction {
            return geometry.prepare_reduced_with_ports(keys,band_hz,true,options.equilibrate_board_mass,
                options.consistent_board_mass,options.edge_cubic_board_mass,reduction,frames);
        }
        if options.edge_cubic_board_mass {
            geometry.prepare_with_motion_edge_cubic_transverse_mass(keys,band_hz,options.equilibrate_board_mass)
        } else if options.consistent_board_mass {
            geometry.prepare_with_motion_consistent_transverse_mass(keys,band_hz,options.equilibrate_board_mass)
        } else if options.equilibrate_board_mass {
            geometry.prepare_with_motion_mass_equilibrated(keys,band_hz)
        } else {geometry.prepare_with_motion(keys,band_hz)}
    }
}
fn prepare_controlled_body(board_text:&str,courses:Vec<geometry::Course>,obj:Option<&str>,mut spec:Specification,
    options:&playback::Options,controls:playback::Controls,continuous:bool)->Result<Scene,String> {
    options.validate()?;
    if obj.is_none() && controls.has_cavity() {
        return Err("a sealed cavity requires a supplied outer enclosure OBJ; bare board-skin exposes the interior board face".into());
    }
    if obj.is_none() {spec.require_board_skin()?;}
    let rigid=options.rigid_assembly.as_deref().map(Assembly::load).transpose()?;
    let keys:Vec<_>=courses.iter().map(|c|c.midi).collect();
    let source_ports=controls.source_ports(&courses)?;
    let board=prepare_board_motion_with_source_ports(board_text,&keys,spec.board_band_hz,
        options,source_ports.as_deref())?;
    let (mut piano,cavity_report)=controls.instrument_with_cavity_report(courses,&board.modes,board.motion.as_ref(),options)?;
    if let Some(report)=cavity_report {spec.source.push_str(&format!("; {report}"));}
    if let Some(c)=board.physical_damping.as_deref() {piano.configure_bare_board_damping(c)?;}
    let (bare,description)=section_skin::boundary(obj,board_text,&spec,
        board.motion.as_ref().ok_or("missing full-vector structural motion")?,continuous)?;
    if piano.has_cavity() {validate_cavity_boundary(&bare)?;}
    let bare=if let Some(rigid)=&rigid {
        let combined=rigid.attach(bare)?;
        spec.source.push_str(&format!("; {}",rigid.report()));
        combined
    } else {bare};
    let boundary=bare.loaded(&piano.bank)?;
    spec.source.push_str(&format!("; {description}"));
    Ok(Scene {piano,board,boundary,spec})
}
/// The sealed chart is below a planar +z board. Its underside cannot also be
/// an exterior moving source. Full enclosure/volume agreement remains explicit
/// input geometry; this rejects the known duplicated-fluid configuration.
fn validate_cavity_boundary(boundary:&Boundary)->Result<(),String> {
    for (panel,normal) in boundary.surface.normals().iter().enumerate() {
        if normal[2] < -1e-10 && boundary.weights.iter().any(|row|row[panel]!=0.0) {
            return Err("sealed cavity exterior cannot expose a moving board underside; supply the outer enclosure with rigid bottom and walls".into());
        }
    }
    Ok(())
}
/// Preserve the source FE certificates separately from the reduced pencil.
/// These comments describe the prepared basis, not an acoustic error estimate.
fn board_reduction_report(board:&board_geometry::PreparedBoard)->String {
    let Some(report)=&board.reduction else {return String::new();};
    let mut out=format!("# soundboard reduction: source_modes={}, retained_modes={}, protected_low_modes={}; bridge static/harmonic targets={:?} Hz (primary and any supplied secondary directions)\n# nonzero_snapshots={}, maximum_relative_displacement_projection_error={:.17e}; no transfer, acoustic or mesh-convergence certificate\n",
        report.source_modes,board.modes.len(),report.protected_low_modes,report.sample_hz,
        report.snapshot_count,report.max_relative_snapshot_error);
    for (i,(lo,hi)) in report.source_frequency_intervals_hz.iter().enumerate() {
        writeln!(out,"# source FE mode {i}: [{lo:.17e},{hi:.17e}] Hz").unwrap();
    }
    for (i,(lo,hi)) in board.frequency_intervals_hz.iter().enumerate() {
        let scope=if i<report.protected_low_modes {"source FE"}else{"projected-pencil Ritz"};
        writeln!(out,"# retained board mode {i} ({scope}): [{lo:.17e},{hi:.17e}] Hz").unwrap();
    }
    out
}
fn response(scene:&Scene)->Result<String,String> {
    if scene.piano.has_cavity() {
        return Err("prescribed-motion response has no cavity reaction; use admittance or rendering".into());
    }
    let samples=scene.boundary.sample(&scene.spec)?;
    let mut csv=format!("# finite exterior BEM, exp(-i omega t), Pa per unit mass-normalized modal acceleration\n# source: {}\n# structure: {}\n# retained string coordinates={}, two transverse directions={}; acoustic motion-to-pressure map, not force-driven structural response\n# panels={}, components={}, min_panels_per_wavelength={}, condition_lower_bound_max={}\n# rigid scatterers, one-way acoustics; no radiation loading, flexible cabinet, room or above-band claim\nfrequency_hz,receiver,input,real_pa_per_acceleration,imag_pa_per_acceleration\n",
        scene.spec.source,scene.board.provenance,scene.piano.bank.modes.len(),scene.piano.bank.has_secondary_polarization(),
        scene.boundary.surface.areas().len(),scene.boundary.components,
        samples.minimum_ppw,samples.maximum_condition_lower_bound);
    for (f,w) in samples.omega.iter().enumerate() {for (receiver,inputs) in samples.values.iter().enumerate() {
        for (input,row) in inputs.iter().enumerate() {
            writeln!(csv,"{:.17e},{receiver},{input},{:.17e},{:.17e}",w/std::f64::consts::TAU,row[f].re,row[f].im).unwrap();
        }
    }}
    Ok(format!("{}{csv}",board_reduction_report(&scene.board)))
}
/// Only admittance admits the explicit conservative-structure comparison.
/// Preserve option/value boundaries: a flag cannot repair a missing value.
fn admittance_options(args:&[String])->Result<(playback::Options,bool),String> {
    let mut numeric=Vec::new();let mut damping=true;let mut args=args.iter();
    while let Some(flag)=args.next() {
        if flag=="--lossless-structure" {
            if !damping {return Err("duplicate --lossless-structure".into());}
            damping=false;
        } else if matches!(flag.as_str(),"--equilibrate-board-mass"|"--consistent-board-mass"|"--edge-cubic-board-mass"|"--rt0425-string-damping") {
            numeric.push(flag.clone());
        } else {
            numeric.push(flag.clone());
            let value=args.next().filter(|s|!s.starts_with("--"))
                .ok_or_else(||format!("missing value for {flag}"))?;
            numeric.push(value.clone());
        }
    }
    let options=playback::Options::harmonic(&numeric)?;
    if !damping && options.rt0425_string_damping {
        return Err("--lossless-structure excludes --rt0425-string-damping".into());
    }
    Ok((options,damping))
}
fn admittance(board_text:&str,courses:&[geometry::Course],obj:&str,spec:&Specification,drive:u8)->Result<String,String> {
    admittance_controlled(board_text,courses,obj,spec,drive,&playback::Options::default(),true)
}
fn admittance_controlled(board_text:&str,courses:&[geometry::Course],obj:&str,spec:&Specification,
    drive:u8,options:&playback::Options,damping:bool)->Result<String,String> {
    admittance_controlled_body(board_text,courses,Some(obj),spec,drive,options,damping,false)
}
fn admittance_controlled_body(board_text:&str,courses:&[geometry::Course],obj:Option<&str>,spec:&Specification,
    drive:u8,options:&playback::Options,damping:bool,continuous:bool)->Result<String,String> {
    options.validate_harmonic()?;
    if !damping && options.rt0425_string_damping {
        return Err("--lossless-structure excludes --rt0425-string-damping".into());
    }
    if obj.is_none() && options.cavity.is_some() {
        return Err("a sealed cavity requires a supplied outer enclosure OBJ; bare board-skin exposes the interior board face".into());
    }
    if obj.is_none() {spec.require_board_skin()?;}
    let rigid=options.rigid_assembly.as_deref().map(Assembly::load).transpose()?;
    let keys:Vec<_>=courses.iter().map(|c|c.midi).collect();
    if !keys.contains(&drive) {return Err("admittance drive key is absent from the scale".into());}
    // Admit the complete supplied frame card before the board eigensolve.
    // Its projection must use the SAME retained motion as played preparation.
    let polarization=options.string_polarization.as_deref().map(|path|
        string_polarization::Specification::load(path,courses)).transpose()?;
    let cavity=options.cavity.as_deref().map(cavity::Specification::load).transpose()?;
    let source_ports=polarization.as_ref().map(|frames|frames.source_ports(courses)).transpose()?;
    let board=prepare_board_motion_with_source_ports(board_text,&keys,spec.board_band_hz,
        options,source_ports.as_deref())?;
    let projected=polarization.as_ref().map(|frames|
        frames.project(courses,&board.modes,board.motion.as_ref())).transpose()?;
    let cavity=cavity.as_ref().map(|cavity|cavity.project(&board)).transpose()?;
    let mut model=bridge_response::BridgeResponse::new_with_string_damping(courses,&board.modes,
        RATE*options.substeps as u32,0.45*f64::from(RATE),options.modes,damping,
        projected.as_ref().map(|frames|frames.secondary().0),options.rt0425_string_damping)?;
    if let Some(c)=board.physical_damping.as_deref() {model.configure_bare_board_damping(c)?;}
    if let Some(cavity)=&cavity {model.configure_cavity(&cavity.loaded(model.bank())?)?;}
    let (bare,description)=section_skin::boundary(obj,board_text,spec,
        board.motion.as_ref().ok_or("missing harmonic surface motion")?,continuous)?;
    if model.has_cavity() {validate_cavity_boundary(&bare)?;}
    let (bare,description)=if let Some(rigid)=&rigid {
        (rigid.attach(bare)?,format!("{description}; {}",rigid.report()))
    } else {(bare,description)};
    let boundary=bare.loaded(model.bank())?;
    let csv=exterior_loading::sweep(&boundary,&model,spec,drive)?;
    let cavity_report=cavity.as_ref().map_or_else(String::new,|cavity|
        format!("# {}; retained in both coupled and one-way columns\n",cavity.report()));
    Ok(format!("{}{cavity_report}# acoustic geometry: {description}\n# structure: {}\n# structural loss: {}; radiation loading remains in the coupled columns\n# intrinsic string loss: {}\n# two transverse directions={}; bridge force and reported velocity use the primary hammer direction\n# board modes={}, retained string coordinates={}, omitted high-frequency duplex mode sets={}\n{}",
        board_reduction_report(&board),board.provenance,if damping {"physical wood damping and selected intrinsic string law"}else{"explicitly disabled by --lossless-structure"},
        if !damping {"disabled"}else if options.rt0425_string_damping {"RT-0425 per-key R_u and eta_u"}else{"estimated common law"},
        model.bank().has_secondary_polarization(),
        board.modes.len(),model.bank().modes.len(),model.bank().omitted_duplex_modes,csv))
}
/// Admit both acoustic realizations before attaching the load or dispatching
/// any score event. Loaded and one-way output share the original PCM path.
fn bake(scene:&mut Scene,loaded:bool)->Result<(exterior_audio::Baked,exterior_geometry::Samples,String),String> {
    scene.spec.require_audio_receivers()?;
    let (load,samples)=if loaded {
        let (fit,samples)=radiation_fit::prepare(&scene.boundary,&scene.spec,scene.piano.bank.rate)?;
        (Some(fit),samples)
    } else {(None,scene.boundary.sample(&scene.spec)?)};
    let baked=exterior_audio::Baked::from_samples(&samples,scene.spec.fit_order)?;
    let report=if let Some(fit)=load {
        scene.piano.configure_radiation(&fit.model)?;
        format!("Passive load: {} acoustic coordinates; complex matrix peak/RMS={:.6}/{:.6}; resistance peak/RMS={:.6}/{:.6}. These sampled bounds do not certify time-step or spatial convergence.",
            fit.model.poles.len(),fit.peak_error,fit.rms_error,fit.resistance_peak_error,fit.resistance_rms_error)
    } else {String::from("One-way exterior reference: no exterior radiation reaction on the mechanics.")};
    Ok((baked,samples,report))
}
fn run(args:&[String])->Result<(),String> {
    match args {
        []=>{println!("{USAGE}");Ok(())},
        [help] if help=="--help" || help=="-h"=>{println!("{USAGE}");Ok(())},
        [command,input,output] if command=="export-rigid"=>{
            if std::path::Path::new(output).exists() {return Err("output must be a fresh path".into());}
            let assembly=Assembly::load(input)?;
            publish(output,assembly.obj().as_bytes())?;
            println!("Written {output}: {} Geometry only, no BEM or render certificate.",assembly.report());
            Ok(())
        },
        [command,board,strings,spec,output,tail @ ..] if command=="export-skin"=>{
            let continuous=match tail {
                []=>false,[flag] if flag=="--continuous-thickness"=>true,
                _=>return Err("export-skin accepts only optional --continuous-thickness".into()),
            };
            section_skin::export(board,strings,spec,output,continuous)
        },
        [command,board,strings,obj,spec,drive,output,tail @ ..] if command=="admittance"=>{
            let (options,damping)=admittance_options(tail)?;
            if options.rt0425_string_damping && strings!="steinway-d" {
                return Err("RT-0425 source laws require the steinway-d source scale".into());
            }
            let drive:u8=drive.parse().map_err(|_|"admittance requires a MIDI bridge key in 21..108")?;
            if !(21..=108).contains(&drive) {return Err("admittance drive key outside 21..108".into());}
            if std::path::Path::new(output).exists() {return Err("output must be a fresh path".into());}
            let spec=Specification::read(&read_bounded(spec,exterior_geometry::MAX_SPEC_BYTES)?)?;
            let courses=scale(strings)?;
            if !courses.iter().any(|c|c.midi==drive) {return Err("admittance drive key is absent from the scale".into());}
            let board=read_bounded(board,8*1024*1024)?;
            let (obj,continuous)=section_skin::read_body(obj)?;
            let csv=admittance_controlled_body(&board,&courses,obj.as_deref(),&spec,drive,&options,damping,continuous)?;
            publish(output,csv.as_bytes())?;
            println!("Written {output}: radiation-loaded bridge mobility and pressure per 1 N peak, all retained strings and physical loss channels. No time-domain feedback or measured-fidelity claim.");
            Ok(())
        }
        [command,board,strings,obj,spec,output,tail @ ..] if command=="response" || command=="render" || command=="render-loaded"=>{
            let (frames,options)=match command.as_str() {
                "response"=>{
                    let options=playback::Options::harmonic(tail)?;
                    if options.cavity.is_some() {
                        return Err("prescribed-motion response has no cavity reaction; use admittance or rendering".into());
                    }
                    (None,options)
                },
                "render"|"render-loaded" if !tail.is_empty()=>
                    (Some(mesh_render::frames(&tail[0])?),playback::Options::parse(&tail[1..])?),
                _=>return Err(USAGE.into()),
            };
            if (options.rt0425_hammer_stiffness || options.rt0425_string_damping)
                && strings!="steinway-d" {
                return Err("RT-0425 source laws require the steinway-d source scale".into());
            }
            if std::path::Path::new(output).exists() {return Err("output must be a fresh path".into());}
            let spec=Specification::read(&read_bounded(spec,exterior_geometry::MAX_SPEC_BYTES)?)?;
            if frames.is_some() {spec.require_audio_receivers()?;}
            if command=="render-loaded" && spec.frequencies<33 {
                return Err("render-loaded requires at least 33 odd-grid frequencies".into());
            }
            let courses=scale(strings)?;let keys:Vec<_>=courses.iter().map(|c|c.midi).collect();
            // Score admission precedes structural/BEM preparation and all writes.
            let score=frames.map(|n|options.score(&keys,n as u64)).transpose()?;
            let controls=playback::Controls::load(&options,&courses)?;
            let geometry=read_bounded(board,8*1024*1024)?;
            let (obj,continuous)=section_skin::read_body(obj)?;
            let mut scene=prepare_controlled_body(&geometry,courses,obj.as_deref(),spec,&options,controls,continuous)?;
            if let (Some(n),Some(score))=(frames,score) {
                let (baked,samples,load_report)=bake(&mut scene,command=="render-loaded")?;
                let score_report=score.report;
                let audio=exterior_audio::render(&mut scene.piano,score.performance,n,&baked,scene.spec.full_scale_pa)?;
                publish(output,&audio.wav)?;
                let physical_report=options.report(&scene.piano);
                print!("{}",board_reduction_report(&scene.board));
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
    #[test]
    fn harmonic_arrays_refuse_playback_before_opening_structure_or_preparing_acoustics() {
        let stamp=std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
        let dir=std::env::temp_dir().join(format!("fs-piano-array-{}-{stamp}",std::process::id()));
        std::fs::create_dir(&dir).unwrap();let name=|file:&str|dir.join(file).to_str().unwrap().to_owned();
        let spec=format!("{}receiver-m,0.06,0.05,1\nreceiver-m,0.07,0.05,1\n",exterior_geometry::tests::specification());
        std::fs::OpenOptions::new().create_new(true).write(true).open(name("array.fspe")).unwrap()
            .write_all(spec.as_bytes()).unwrap();
        for command in ["render","render-loaded"] {
            let output=name(&format!("{command}.wav"));
            let args=vec![command.into(),name("missing-board.fsb"),name("missing-scale.csv"),
                name("missing-body.obj"),name("array.fspe"),output.clone(),"0.05".into()];
            assert!(run(&args).unwrap_err().contains("one or two receivers"));
            assert!(!std::path::Path::new(&output).exists());
        }
        let mut scene=small_source_scene();scene.spec.receivers.push([0.06,0.05,1.]);
        for loaded in [false,true] {
            assert!(bake(&mut scene,loaded).err().unwrap().contains("one or two receivers"));
            assert!(!scene.piano.has_radiation());
            assert!(scene.piano.bank.q.iter().chain(&scene.piano.bank.v).all(|x|*x==0.));
        }
    }
    #[test]
    fn published_laws_require_the_source_scale_at_the_cli_boundary() {
        for flag in ["--rt0425-hammer-stiffness", "--rt0425-string-damping"] {
            let args=["render","missing.fsb","other-scale.csv","missing.obj",
                "missing.fspe","unused.wav","0.05",flag].map(str::to_owned);
            assert!(run(&args).unwrap_err().contains("steinway-d source scale"));
        }
        assert!(playback::Options::harmonic(&["--rt0425-hammer-stiffness".into()]).is_err());
    }
    #[test]
    fn admittance_accepts_the_same_opt_in_board_scaling_as_response() {
        let (options,damping)=admittance_options(&[
            "--equilibrate-board-mass".into(), "--modes".into(), "24".into(),
            "--lossless-structure".into()]).unwrap();
        assert!(options.equilibrate_board_mass);
        assert_eq!(options.modes,24);
        assert!(!damping);
        assert!(admittance_options(&[
            "--equilibrate-board-mass".into(), "--equilibrate-board-mass".into()]).is_err());
        for flag in ["--consistent-board-mass", "--edge-cubic-board-mass"] {
            let (options,damping)=admittance_options(&[
                flag.into(), "--equilibrate-board-mass".into(), "--lossless-structure".into()]).unwrap();
            assert!(options.equilibrate_board_mass);
            assert_eq!(options.edge_cubic_board_mass,flag=="--edge-cubic-board-mass");
            assert_eq!(options.consistent_board_mass,flag=="--consistent-board-mass");
            assert!(!damping);
            assert!(admittance_options(&["--modes".into(),flag.into()]).is_err());
        }
        assert!(admittance_options(&["--consistent-board-mass".into(),"--edge-cubic-board-mass".into()]).is_err());
    }
    #[test]
    fn exterior_inertia_selection_retains_the_selected_motion_and_refuses_crowned_substitution() {
        let (board,_,_,spec)=small_source_inputs();
        for flag in ["--consistent-board-mass", "--edge-cubic-board-mass"] {
            let options=playback::Options::parse(&[flag.into(),"--equilibrate-board-mass".into()]).unwrap();
            let prepared=prepare_board_motion(&board,&[69],spec.board_band_hz,&options).unwrap();
            assert_eq!(prepared.motion.as_ref().unwrap().is_edge_cubic(),flag=="--edge-cubic-board-mass");
            assert!(prepared.provenance.contains("mass-diagonal solver equilibration"));
            assert!(prepare_board_motion(crowned_board::HEADER,&[69],spec.board_band_hz,&options)
                .err().unwrap().contains("flat geometric board"));
        }
    }
    #[test]
    fn opt_in_board_scaling_reaches_the_exterior_scene_and_is_reported() {
        let (board,courses,obj,spec)=small_source_inputs();
        let options=playback::Options::parse(&["--equilibrate-board-mass".into()]).unwrap();
        let controls=playback::Controls::from_texts(&courses,None,None,None).unwrap();
        let scene=prepare_controlled_body(&board,courses,Some(&obj),spec,&options,controls,false).unwrap();
        assert!(scene.board.provenance.contains("mass-diagonal solver equilibration"));
        assert!(!scene.board.modes.is_empty());
        assert!(!scene.boundary.surface.areas().is_empty());
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
    fn projected_vector_strings_reach_the_same_finite_body_renderer() {
        for edge_cubic_board_mass in [false,true] {
        let (board,courses,obj,spec)=small_source_inputs();
        let board=board.replace("node,4,0.05,0.05", "node,4,0.043,0.054");
        let frames=format!("{}\nsource,estimated,test bridge height and frame\ncourse,69,0,0,0,1,0,0,0.02,0,1,0,0,0,1,0.3\n",
            string_polarization::HEADER);
        let options=playback::Options {edge_cubic_board_mass,..playback::Options::default()};
        let controls=playback::Controls::from_texts(&courses,None,None,Some("estimated")).unwrap()
            .with_string_polarization(&frames,&courses).unwrap();
        let mut scene=prepare_controlled_body(&board,courses,Some(&obj),spec,&options,controls,false).unwrap();
        assert!(scene.piano.bank.has_secondary_polarization());
        assert!(scene.piano.bank.strings.iter().filter(|s|s.polarization==1)
            .any(|s|s.bridge.iter().any(|g|g.abs()>1e-8)));
        let (baked,_,_)=bake(&mut scene,false).unwrap();
        let score=performance::Performance::read("sample,event,key,value\n0,note_on,69,0.5\n1200,note_off,69,0\n",&[69],2400).unwrap();
        let audio=exterior_audio::render(&mut scene.piano,score,2400,&baked,2.).unwrap();
        assert_eq!(u16::from_le_bytes([audio.wav[22],audio.wav[23]]),2);
        assert!(audio.peak_pa>1e-14);
        assert!(scene.piano.accounting.felt_loss_j>0. && scene.piano.accounting.damper_loss_j>0.);
        assert!(scene.piano.bank.strings.iter().filter(|s|s.polarization==1).any(|s|
            scene.piano.bank.q[s.modes.clone()].iter().any(|q|q.abs()>1e-14)));
        assert!((scene.piano.accounting.input_work_j-scene.piano.energy_j()
            -scene.piano.accounting.dissipated_j()).abs()<1e-7);
        }
    }

    #[test]
    fn published_hammer_contact_reaches_finite_body_pcm_and_passive_loss() {
        let options=playback::Options::parse(&[
            "--rt0425-hammer-stiffness".into(), "--rt0425-hammer-dissipation".into()]).unwrap();
        for loaded in [false,true] {
            let (board,courses,obj,spec)=small_source_inputs();
            let (_,_,_,baseline_spec)=small_source_inputs();
            let controls=playback::Controls::load(&options,&courses).unwrap();
            let mut selected=prepare_controlled_body(&board,courses.clone(),Some(&obj),spec,
                &options,controls,false).unwrap();
            let mut baseline=prepare(&board,courses,&obj,baseline_spec).unwrap();
            let (baked,_,_)=bake(&mut selected,loaded).unwrap();
            bake(&mut baseline,loaded).unwrap();
            let score=||performance::Performance::read(
                "sample,event,key,value\n0,note_on,69,2.48\n1200,note_off,69,0\n",
                &[69],2400).unwrap();
            let source=exterior_audio::render(&mut selected.piano,score(),2400,&baked,2.).unwrap();
            let old=exterior_audio::render(&mut baseline.piano,score(),2400,&baked,2.).unwrap();
            assert_ne!(source.wav,old.wav);
            assert!(selected.piano.accounting.felt_relaxation_loss_j>0.0);
            assert!(source.peak_pa>0.0);
            assert_eq!(selected.piano.has_radiation(),loaded);
            let defect=selected.piano.accounting.input_work_j-selected.piano.energy_j()
                -selected.piano.accounting.dissipated_j();
            assert!(defect.abs()<1e-7);
        }
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
#[path="grand_piano/cavity_admittance_tests.rs"]
mod cavity_admittance_tests;

#[cfg(test)]
#[path="grand_piano/cavity_exterior_tests.rs"]
mod cavity_exterior_tests;

#[cfg(test)]
#[path="grand_piano/radiation_render_tests.rs"]
mod radiation_render_tests;

#[cfg(test)]
#[path="grand_piano/lossless_admittance_tests.rs"]
mod lossless_admittance_tests;

#[cfg(test)]
#[path="grand_piano/rigid_assembly_render_tests.rs"]
mod rigid_assembly_render_tests;

#[cfg(test)]
#[path="grand_piano/harmonic_controls_render_tests.rs"]
mod harmonic_controls_render_tests;

#[cfg(test)]
#[path="grand_piano/board_reduction_exterior_tests.rs"]
mod board_reduction_exterior_tests;
