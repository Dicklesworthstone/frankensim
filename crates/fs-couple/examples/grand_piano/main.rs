//! Geometry-first grand-piano composition of existing FrankenSim owners.
//! `cargo run --release -p fs-couple --example grand_piano -- --preset steinway-d --render piano.wav`
//! Source-derived geometry is approximate, not a measured specimen digital twin.
mod geometry;
mod linear;
mod board;
mod board_geometry;
mod crowned_board;
mod steinway_d;
mod steinway_scale;
mod performance;
use performance::midi;
mod felt;
mod engine;
mod microphone;
mod pressure_basis;
mod audio;
mod hammer_materials;
mod string_polarization;

const USAGE: &str = "grand_piano [--render piano.wav] [--scale strings.csv]
    [--preset steinway-d] [--board board.csv | --board-geometry panel.fsb|panel.fss]
    [--rt0425-bridge-contacts] [--rt0425-hammer-stiffness]
    [--rt0425-hammer-dissipation] [--rt0425-string-damping]
    [--equilibrate-board-mass] [--consistent-board-mass | --edge-cubic-board-mass]
    [--hammers materials.fsh] [--hammer-footprints faces.fshp]
    [--dampers estimated | pads.fspd] [--string-stretching axial.fspx]
    [--string-polarization bridge-frames.fspp]
    [--concert-pitch 430..450 | --raw-tensions]
    [--mesh-divisions 4..32] [--dump-geometry panel.fsb] [--dump-obj soundboard.obj]
    [--board-band-hz Hz] [--board-reduction max_modes,keep_low_modes,Hz,...]
    [--performance events.csv] [--observer-gain Pa/(m^3/s)]
    [--midi performance.mid] [--midi-channel 1..16]
    [--midi-velocity-max-m-s V] [--midi-half-pedal]
    [--microphone x_m,y_m,z_m] [--microphone-right x_m,y_m,z_m] [--diagnostic-volume]
    [--acoustic-refinement-levels 0..3]
    [--note 21..108] [--velocity m/s] [--duration seconds]
    [--bridge-trace-csv paired.csv]
    [--modal-pressure-csv modal.csv]
    [--receiver-pressure-csv receivers.csv]
    [--pressure-basis-json basis.json]
    [--sample-rate Hz] [--substeps 1..16] [--modes 1..512]
    [--pcm-full-scale-pa positive-Pa]
    [--dump-scale strings.csv] [--dump-board board.csv]
--preset steinway-d reconstructs the published 17-rib Model D drawing, with
spruce panel, sugar-pine ribs, maple bridges, cut-off bar and 88 bridge stations.
--board-geometry may replace that board while retaining the preset strings,
felt cards and shank mechanics. Its native header selects flat FSB or crowned
3-D shell FSS. An invalid/missing supplied board never falls back to the preset.
Mesh-generation/export controls cannot accompany a supplied board override.
The preset cannot be combined with --board modal CSV.
--rt0425-bridge-contacts projects the source's 84 published coupling points
onto this approximate board's bridges; four end keys remain extrapolated.
--rt0425-hammer-stiffness applies the source K_H to each unison string.
--rt0425-hammer-dissipation requires that stiffness and uses the published
per-key R_H d(e^p)/dt instead of estimated crush and Prony relaxation.
--rt0425-string-damping uses the published per-key R_u and eta_u for the
preset scale's intrinsic string losses instead of the estimated common law.
--equilibrate-board-mass solves a flat geometric board in mass-diagonal
coordinates, then certifies modes in the original SI coordinates. It is an
opt-in numerical trial, not a different soundboard or a tuned piano preset.
--consistent-board-mass integrates the flat panel's P1 transverse inertia
exactly; slope inertia remains lumped and beam inertia is lumped by default. It is an opt-in numerical
trial, not a measured Model D material correction.
--edge-cubic-board-mass integrates a declared cubic panel displacement field
and applies that same field at bridge and acoustic surface samples. This is
an opt-in numerical trial; slope inertia remains lumped and beam inertia is lumped by default.
The optional FSB stiffener-mass row accepts lumped, consistent-hermite or
consistent-eccentric. The last adds the supplied bending rotary and offset
centroid inertia to Hermite translation; it does not infer torsional polar inertia.
These opt-in corrections have not passed a perceptual similarity gate.
--acoustic-refinement-levels uniformly subdivides flat P1 radiating triangles
for Rayleigh integration only. It preserves the structural mesh, modes and
bridge mechanics. Level 0 is unchanged; the 120000-point receiver budget still
applies. Crowned/cubic fields and diagnostic-volume output are excluded.
It uses Chabassier/Durufle's wrapped-string MODEL table (84 notes plus four
estimated extensions), separate per-key hammer force and relaxation cards, and
published shank geometry reduced to rigid rotation plus one bending coordinate.
The shank is not a fitted oscillator: mass, compliance and jack projection come
from its dimensions/material. This is a linearized reduction, not the full action.
Lengths, effective winding mass and EI remain fixed. Preset tensions are tuned
to A4=440 Hz by default to compensate rounded source values; --raw-tensions
preserves the published table. --concert-pitch tunes first partials by changing
physical tension, not oscillator frequencies. This is not a stretch-tuning fit.
An explicit --scale always supplies the geometry/masses; absent --concert-pitch,
its tensions are preserved even with a preset. Preset hammer voicing still applies
unless --hammers supplies a complete per-key WoolFelt/Prony material file.
--hammers requires --render and a card for every key in the supplied scale, not
only the struck keys. Missing, duplicate or invalid cards refuse without fallback.
These cards change physical contact forces and relaxation, not an output EQ.
The preset shank and the scale's hammer mass/patch geometry remain unchanged.
See HAMMERS.md for the SI format; importing values does not certify measurements.
--hammer-footprints supplies point, uniform span or authored crown profile for
EVERY key. Profiles supply one to four ordered sites with longitudinal offsets,
face recession, local felt thickness and positive fractions of the original area.
Sites engage according to their gaps and retain independent felt/Prony histories
while sharing the hammer. Point/span thickness still comes from the scale.
It requires --render and works with supplied materials, MIDI and existing pedals.
Published R_H requires the original uniform thickness. Varying thickness applies
the selected felt/Prony law locally and cannot silently reuse that source rate.
See HAMMER_FOOTPRINTS.md for SI rows and the parallel-column contact model.
--dampers selects finite-footprint viscous pads instead of the default point
damper. 'estimated' declares approximate spans and drag; a file must cover
every scale key with a pad or explicit free row. It requires --render and
uses existing MIDI/CSV key, sustain and sostenuto controls. See DAMPERS.md.
This is spatial drag, not falling-pad or hysteretic felt contact mechanics.
--string-polarization adds both transverse directions of each physical string.
Its complete per-key file supplies bridge sites, 3-D arms, string/hammer axes
and lateral damper ratios. The same board solve projects both bridge rows;
the primary row must agree with the board's existing bridge geometry.
It requires a geometric --render, using the selected P1, edge-cubic or crowned
motion field. Modal CSV lacks that motion. No lateral coupling or drag is guessed.
See STRING_POLARIZATION.md for the physical input format and scope.
--string-stretching supplies linear or geometric-extension selection for EVERY
scale key. A stretch row supplies axial rigidity EA in N and a moderate-slope
bound; these are not inferred from tension, EI or winding mass. It requires
--render, retains the same strings/bridge/hammer/pedal states, and solves
extension and felt contact in the same mechanical tick. All-linear selection
preserves the original image. No pitch automation, clipping, modal retuning or
output processing is substituted; see STRING_STRETCHING.md for the SI format.
--dump-geometry/--dump-obj export the board model; export alone skips eigenanalysis.
Thickness taper, material constants and key assignment include explicit estimates.
Default per-key WoolFelt loading envelopes are source-derived; crush/unloading
parameters, felt patch geometry and Prony time constants remain estimates.
--note performs a single-key study; otherwise the demo also plays a chord of
available keys. --velocity overrides the three demo hammer launch speeds.
Velocity is POST-ESCAPEMENT hammer velocity, not MIDI velocity or key motion.
Geometric boards assemble a flat orthotropic plate or a supplied crowned shell.
Crowned structures retain 3-D motion, grain and eccentric ribs/bridges; their
radiation is a projected flat-baffle approximation, not 3-D exterior BEM.
Crown is reference geometry, not solved downbearing; see CROWNED_BOARD.md.
The explicit frequency band ordinarily admits at most 128 board modes.
--board-reduction explicitly permits up to 512 source modes on a flat or crowned board,
then retains at most max_modes (1..128), including keep_low_modes exact low
modes. Static and damped responses at 1..16 increasing target frequencies guide
the remaining basis at every admitted primary bridge. Targets must lie within
--board-band-hz; that source band is never widened implicitly. Full projected
wood damping, bridge motion and the radiating field stay in the same basis.
For supplied downbearing, reduction uses the solved equilibrium tangent modes.
This requires --render and excludes --dump-board: modal CSV cannot preserve the
dense material operator. Flat-board mass controls remain flat-only. See BOARD_REDUCTION.md.
--modes independently admits up to 512 string partials. These budgets and the
reported snapshot projection error do not certify transfer, acoustic or mesh
convergence, perceptual similarity, or real-time performance.
--performance uses sample,event,key,value CSV instead of the demo and cannot be
combined with --note or --velocity. note_on values are hammer velocity in m/s.
With the preset, jack_staccato and jack_legato instead take peak force in N at
the physical jack station, with 7/100 ms pulses and 1.5 mm let-off. For example:
0,jack_staccato,27,70
48000,jack_legato,69,30
Use note_off events to release keys before repeating them. Shank damping and
the backcheck remain estimates; jack timing is resolved at the mechanical rate.
--midi replaces the demo/--performance with a format-0/1 Standard MIDI File.
It selects channel 1 by default. Velocity 127 maps to 4.5 m/s by default;
--midi-velocity-max-m-s explicitly changes this linear, uncalibrated hammer
launch mapping, not audio gain. CC64/66/67 drive the existing three pedals;
--midi-half-pedal opts into CC64/127 travel instead of the default switch.
The complete score must fit --duration, leaving room for acoustic ringdown.
Unsupported physical controls refuse or are reported as ignored; see MIDI.md.
Geometric boards default to a spatial Rayleigh half-space pressure microphone
at (0.675,1,1) metres in the mesh coordinate system. --microphone moves it.
--microphone-right adds a second spatial receiver and writes left/right stereo
from one mechanical performance. The original/default microphone is left;
no panning, duplicated mono, peak normalization or second physics run. See STEREO.md.
This assumes an infinite baffle, with no lid/room scattering or air backreaction.
Pressure uses every mechanics substep and causal anti-alias filtering before
output-rate propagation. At 4x oversampling the filter adds 44 audio samples of
latency, in addition to acoustic travel time. Histories persist across blocks.
--bridge-trace-csv requires --note or an isolated --performance and geometric pressure. It
writes that key's modeled vertical bridge velocity, a centered-difference
acceleration, and the left pressure sample on the same output clock. The
acceleration is a diagnostic derivative of output-rate velocity, not a sensor
model; pressure retains its filter and travel-time delay. Performance traces
require exactly one scheduled hammer launch or jack pulse, with releases only
for that key. --note remains a sustain/restrike demo, not a long isolated strike.
--modal-pressure-csv requires the same single-note geometric render. It writes
each loaded-board mode's delayed pressure at each receiver on the WAV clock;
the signed modal sum reconstructs the pressure, including cancellation. Modes
are indexed in the loaded basis, not identified as bare-board eigenfrequencies.
--receiver-pressure-csv writes unquantized total Pa at each physical receiver
on the WAV clock, for isolated notes or music, without allocating modal traces.
--pressure-basis-json exports the exact bare-from-loaded microphone map and
loaded diagonal reference frequencies for that render. Its columns match
modal pressure indices; these frequencies are not full coupled instrument poles.
--diagnostic-volume retains the old volume-velocity observer; --observer-gain
applies only to that diagnostic, not to physical microphone pressure.
--pcm-full-scale-pa declares the pressure mapped to PCM full scale (default 2 Pa).
All output paths must be fresh, with existing parent directories. Equivalent
paths through dot components or symlinked parents refuse before preparation.
An over-range render refuses before writing a clipped WAV; this option changes
encoding gain, not the mechanics, microphone position, or acoustic calibration.";

#[derive(Debug)]
struct Options {
    render: Option<String>, scale: Option<String>, board: Option<String>,
    board_geometry: Option<String>, performance: Option<String>, preset: Option<String>,
    hammers: Option<String>, hammer_footprints: Option<String>, dampers: Option<String>,
    string_stretching: Option<String>, string_polarization: Option<String>,
    midi: Option<String>, midi_mapping: midi::Mapping,
    concert_pitch: Option<f64>, raw_tensions: bool,
    rt0425_bridge_contacts: bool, rt0425_hammer_stiffness: bool,
    rt0425_hammer_dissipation: bool, rt0425_string_damping: bool,
    equilibrate_board_mass: bool, consistent_board_mass: bool, edge_cubic_board_mass: bool,
    mesh_divisions: usize, dump_geometry: Option<String>, dump_obj: Option<String>,
    board_band_hz: f64, observer_gain: f64,
    board_reduction: Option<board_geometry::ritz::RitzOptions>,
    microphone: Option<[f64; 3]>, microphone_right: Option<[f64; 3]>, diagnostic_volume: bool,
    acoustic_refinement_levels: usize,
    dump_scale: Option<String>, dump_board: Option<String>,
    bridge_trace_csv: Option<String>, modal_pressure_csv: Option<String>, receiver_pressure_csv: Option<String>,
    pressure_basis_json: Option<String>,
    note: Option<u8>, velocity: Option<f64>, duration: f64,
    sample_rate: u32, substeps: usize, modes: usize, pcm_full_scale_pa: f64, help: bool,
}
impl Default for Options {
    fn default() -> Self {
        Self { render: None, scale: None, board: None, board_geometry: None,
            performance: None, preset: None, hammers: None, hammer_footprints: None, dampers: None, concert_pitch: None, raw_tensions: false,
            rt0425_bridge_contacts: false, rt0425_hammer_stiffness: false,
            rt0425_hammer_dissipation: false, rt0425_string_damping: false,
            equilibrate_board_mass: false, consistent_board_mass: false, edge_cubic_board_mass: false,
            string_stretching: None, string_polarization: None,
            midi: None, midi_mapping: midi::Mapping::default(),
            mesh_divisions: 8, dump_geometry: None, dump_obj: None,
            board_band_hz: 400.0, observer_gain: 10_000.0, dump_scale: None,
            board_reduction: None,
            microphone: None, microphone_right: None, diagnostic_volume: false,
            acoustic_refinement_levels: 0,
            dump_board: None, bridge_trace_csv: None, modal_pressure_csv: None, receiver_pressure_csv: None,
            pressure_basis_json: None,
            note: None, velocity: None, duration: 6.0,
            sample_rate: 48_000, substeps: 4, modes: 24, pcm_full_scale_pa: 2.0, help: false }
    }
}
impl Options {
    fn parse(args: &[String]) -> Result<Self, String> {
        let mut options = Self::default();
        let mut seen = std::collections::BTreeSet::new();
        let mut args = args.iter();
        while let Some(flag) = args.next() {
            if flag == "--help" || flag == "-h" { options.help = true; continue; }
            if !seen.insert(flag.as_str()) { return Err(format!("duplicate option {flag}")); }
            if flag == "--diagnostic-volume" { options.diagnostic_volume = true; continue; }
            if flag == "--raw-tensions" { options.raw_tensions = true; continue; }
            if flag == "--rt0425-bridge-contacts" { options.rt0425_bridge_contacts = true; continue; }
            if flag == "--rt0425-hammer-stiffness" { options.rt0425_hammer_stiffness = true; continue; }
            if flag == "--rt0425-hammer-dissipation" { options.rt0425_hammer_dissipation = true; continue; }
            if flag == "--rt0425-string-damping" { options.rt0425_string_damping = true; continue; }
            if flag == "--equilibrate-board-mass" { options.equilibrate_board_mass = true; continue; }
            if flag == "--consistent-board-mass" { options.consistent_board_mass = true; continue; }
            if flag == "--edge-cubic-board-mass" { options.edge_cubic_board_mass = true; continue; }
            if flag == "--midi-half-pedal" { options.midi_mapping.continuous_sustain = true; continue; }
            let value = args.next().ok_or_else(|| format!("missing value for {flag}"))?;
            let invalid = || format!("invalid value for {flag}: {value}");
            match flag.as_str() {
                "--render" => options.render = Some(value.clone()),
                "--scale" => options.scale = Some(value.clone()),
                "--board" => options.board = Some(value.clone()),
                "--board-geometry" => options.board_geometry = Some(value.clone()),
                "--preset" => options.preset = Some(value.clone()),
                "--hammers" => options.hammers = Some(value.clone()),
                "--hammer-footprints" => options.hammer_footprints = Some(value.clone()),
                "--dampers" => options.dampers = Some(value.clone()),
                "--string-stretching" => options.string_stretching = Some(value.clone()),
                "--string-polarization" => options.string_polarization = Some(value.clone()),
                "--concert-pitch" => options.concert_pitch = Some(value.parse().map_err(|_| invalid())?),
                "--mesh-divisions" => options.mesh_divisions = value.parse().map_err(|_| invalid())?,
                "--dump-geometry" => options.dump_geometry = Some(value.clone()),
                "--dump-obj" => options.dump_obj = Some(value.clone()),
                "--performance" => options.performance = Some(value.clone()),
                "--midi" => options.midi = Some(value.clone()),
                "--midi-channel" => {
                    options.midi_mapping.channel = value.parse::<u8>().map_err(|_| invalid())?
                        .checked_sub(1).ok_or_else(invalid)?;
                }
                "--midi-velocity-max-m-s" => options.midi_mapping.maximum_velocity_m_s =
                    value.parse().map_err(|_| invalid())?,
                "--microphone" | "--microphone-right" => {
                    let values = value.split(',').map(str::parse::<f64>).collect::<Result<Vec<_>,_>>()
                        .map_err(|_| invalid())?;
                    if values.len() != 3 { return Err(invalid()); }
                    let position=Some([values[0], values[1], values[2]]);
                    if flag=="--microphone" {options.microphone=position;} else {options.microphone_right=position;}
                }
                "--board-band-hz" => options.board_band_hz = value.parse().map_err(|_| invalid())?,
                "--board-reduction" => options.board_reduction = Some(board_geometry::ritz::RitzOptions::parse(value)?),
                "--acoustic-refinement-levels" => options.acoustic_refinement_levels = value.parse().map_err(|_| invalid())?,
                "--observer-gain" => options.observer_gain = value.parse().map_err(|_| invalid())?,
                "--dump-scale" => options.dump_scale = Some(value.clone()),
                "--dump-board" => options.dump_board = Some(value.clone()),
                "--bridge-trace-csv" => options.bridge_trace_csv = Some(value.clone()),
                "--modal-pressure-csv" => options.modal_pressure_csv = Some(value.clone()),
                "--receiver-pressure-csv" => options.receiver_pressure_csv = Some(value.clone()),
                "--pressure-basis-json" => options.pressure_basis_json = Some(value.clone()),
                "--note" => options.note = Some(value.parse().map_err(|_| invalid())?),
                "--velocity" => options.velocity = Some(value.parse().map_err(|_| invalid())?),
                "--duration" => options.duration = value.parse().map_err(|_| invalid())?,
                "--sample-rate" => options.sample_rate = value.parse().map_err(|_| invalid())?,
                "--substeps" => options.substeps = value.parse().map_err(|_| invalid())?,
                "--modes" => options.modes = value.parse().map_err(|_| invalid())?,
                "--pcm-full-scale-pa" => options.pcm_full_scale_pa = value.parse().map_err(|_| invalid())?,
                _ => return Err(format!("unknown option {flag}\n{USAGE}")),
            }
        }
        if options.note.is_some_and(|n| !(21..=108).contains(&n))
            || options.velocity.is_some_and(|v| !v.is_finite() || v <= 0.0 || v > 8.0)
            || !options.duration.is_finite() || !(0.001..=120.0).contains(&options.duration)
            || !(8_000..=192_000).contains(&options.sample_rate)
            || !(1..=16).contains(&options.substeps) || !(1..=linear::MAX_STRING_MODES).contains(&options.modes) {
            return Err("render control outside its finite admitted range".into());
        }
        if options.concert_pitch.is_some_and(|f| !f.is_finite() || !(430.0..=450.0).contains(&f))
            || (options.raw_tensions && (options.preset.is_none() || options.concert_pitch.is_some())) {
            return Err("concert pitch must be 430..450 Hz; --raw-tensions requires a preset and excludes --concert-pitch".into());
        }
        if options.board.is_some() && (options.board_geometry.is_some() || options.preset.is_some()) {
            return Err("modal board CSV excludes a preset and geometric board; --board-geometry may override a preset board".into());
        }
        if options.preset.as_deref().is_some_and(|p| p != "steinway-d") {
            return Err("unknown piano preset; available: steinway-d".into());
        }
        if (options.rt0425_bridge_contacts && !options.uses_preset_board())
            || (options.rt0425_hammer_stiffness && (options.preset.is_none()
                || options.hammers.is_some() || options.render.is_none()))
            || (options.rt0425_hammer_dissipation && !options.rt0425_hammer_stiffness) {
            return Err("RT-0425 contacts need the preset board; hammer stiffness needs a preset render without --hammers; hammer dissipation also requires source stiffness".into());
        }
        if options.rt0425_string_damping && (options.preset.is_none()
            || options.scale.is_some() || options.render.is_none()) {
            return Err("RT-0425 string damping requires a preset render without a supplied scale".into());
        }
        if options.hammers.is_some() && options.render.is_none() {
            return Err("--hammers requires --render; material input is not an export-only option".into());
        }
        if options.hammer_footprints.as_ref().is_some_and(|s| s.trim().is_empty() || options.render.is_none()) {
            return Err("--hammer-footprints requires --render and a complete nonempty specification path".into());
        }
        if options.dampers.as_ref().is_some_and(|s| s.trim().is_empty() || options.render.is_none()) {
            return Err("--dampers requires --render and either estimated or a nonempty specification path".into());
        }
        if options.string_stretching.as_ref().is_some_and(|s|
            s.trim().is_empty() || s.starts_with("--") || options.render.is_none()) {
            return Err("--string-stretching requires --render and a complete nonempty specification path".into());
        }
        if !(4..=32).contains(&options.mesh_divisions)
            || (!options.uses_preset_board() && (seen.contains("--mesh-divisions")
                || options.dump_geometry.is_some() || options.dump_obj.is_some())) {
            return Err("mesh/export controls require --preset steinway-d without a supplied board override; divisions must be 4..32".into());
        }
        if seen.contains("--board-band-hz") && options.board_geometry.is_none() && options.preset.is_none() {
            return Err("--board-band-hz requires --board-geometry or --preset".into());
        }
        if !options.board_band_hz.is_finite() || options.board_band_hz <= 0.0
            || options.board_band_hz >= 0.45 * f64::from(options.sample_rate)
            || !options.observer_gain.is_finite() || options.observer_gain <= 0.0 {
            return Err("invalid board frequency band or diagnostic observer gain".into());
        }
        if !options.pcm_full_scale_pa.is_finite() || options.pcm_full_scale_pa <= 0.0
            || (seen.contains("--pcm-full-scale-pa") && options.render.is_none()) {
            return Err("--pcm-full-scale-pa requires --render and a positive finite pressure".into());
        }
        if options.performance.is_some()
            && (options.render.is_none() || options.note.is_some() || options.velocity.is_some()) {
            return Err("--performance requires --render and replaces --note/--velocity demo controls".into());
        }
        if options.midi.is_some() && (options.render.is_none() || options.performance.is_some()
            || options.note.is_some() || options.velocity.is_some()) {
            return Err("--midi requires --render and excludes --performance/--note/--velocity".into());
        }
        let mapping = options.midi_mapping;
        if mapping.channel > 15 || !mapping.maximum_velocity_m_s.is_finite()
            || mapping.maximum_velocity_m_s <= 0.0 || mapping.maximum_velocity_m_s > 8.0
            || mapping.maximum_velocity_m_s / 127.0 == 0.0
            || (options.midi.is_none() && ["--midi-channel", "--midi-velocity-max-m-s", "--midi-half-pedal"]
                .iter().any(|flag| seen.contains(flag))) {
            return Err("MIDI controls require --midi, channel 1..16 and finite maximum hammer velocity in (0,8] m/s".into());
        }
        let geometric = options.preset.is_some() || options.board_geometry.is_some();
        if let Some(reduction) = &options.board_reduction {
            if !geometric || options.render.is_none() || options.dump_board.is_some() {
                return Err("--board-reduction requires a geometric --render and excludes --dump-board; modal CSV cannot preserve full damping".into());
            }
            if reduction.sample_hz.iter().any(|hz| *hz > options.board_band_hz) {
                return Err("--board-reduction target frequencies must lie within --board-band-hz".into());
            }
        }
        if options.string_polarization.as_ref().is_some_and(|s| s.trim().is_empty()
            || s.starts_with("--") || options.render.is_none() || !geometric) {
            return Err("--string-polarization requires a complete nonempty specification and a geometric render with full-vector motion; modal CSV is unsupported".into());
        }
        if options.acoustic_refinement_levels > 3 || (seen.contains("--acoustic-refinement-levels")
            && (!geometric || options.render.is_none() || options.diagnostic_volume
                || options.edge_cubic_board_mass)) {
            return Err("--acoustic-refinement-levels requires a flat P1 geometric pressure render and levels 0..3".into());
        }
        if options.equilibrate_board_mass && !geometric {
            return Err("--equilibrate-board-mass requires a flat geometric board".into());
        }
        if options.consistent_board_mass && (!geometric || (options.render.is_none()
            && options.dump_board.is_none())) {
            return Err("--consistent-board-mass requires a flat geometric board solve for --render or --dump-board".into());
        }
        if options.edge_cubic_board_mass && (!geometric || options.consistent_board_mass
            || (options.render.is_none() && options.dump_board.is_none())) {
            return Err("--edge-cubic-board-mass requires a flat geometric board solve for --render or --dump-board and excludes --consistent-board-mass".into());
        }
        if [options.microphone,options.microphone_right].iter().flatten()
            .any(|p| p.iter().any(|x|!x.is_finite()) || p[2]<0.05)
            || ((options.microphone.is_some() || options.microphone_right.is_some())
                && (!geometric || options.render.is_none() || options.diagnostic_volume)) {
            return Err("microphones need a geometric render, finite x,y,z with z>=0.05, and no --diagnostic-volume".into());
        }
        if geometric && !options.diagnostic_volume && seen.contains("--observer-gain") {
            return Err("--observer-gain requires --diagnostic-volume for a geometric board".into());
        }
        if options.bridge_trace_csv.is_some()
            && (options.render.is_none() || (options.note.is_none() && options.performance.is_none())
                || !geometric || options.diagnostic_volume) {
            return Err("--bridge-trace-csv requires --render, --note or an isolated --performance, and geometric pressure without --diagnostic-volume".into());
        }
        if options.modal_pressure_csv.is_some()
            && (options.render.is_none() || (options.note.is_none() && options.performance.is_none())
                || !geometric || options.diagnostic_volume) {
            return Err("--modal-pressure-csv requires --render, --note or an isolated --performance, and geometric pressure without --diagnostic-volume".into());
        }
        if options.receiver_pressure_csv.is_some()
            && (options.render.is_none() || !geometric || options.diagnostic_volume) {
            return Err("--receiver-pressure-csv requires a geometric pressure render without --diagnostic-volume".into());
        }
        if options.pressure_basis_json.is_some()
            && (options.render.is_none() || !geometric || options.diagnostic_volume) {
            return Err("--pressure-basis-json requires a geometric pressure render without --diagnostic-volume".into());
        }
        // Do not overwrite the very measurements that a render was asked to use.
        let inputs = [options.scale.as_ref(), options.board.as_ref(),
            options.board_geometry.as_ref(), options.performance.as_ref(), options.hammers.as_ref(), options.hammer_footprints.as_ref(), options.midi.as_ref(),
            options.dampers.as_ref().filter(|s| s.as_str() != "estimated"),
            options.string_stretching.as_ref(), options.string_polarization.as_ref()];
        let outputs = [options.render.as_ref(), options.dump_scale.as_ref(), options.dump_board.as_ref(),
            options.dump_geometry.as_ref(), options.dump_obj.as_ref(), options.bridge_trace_csv.as_ref(),
            options.modal_pressure_csv.as_ref(), options.receiver_pressure_csv.as_ref(),
            options.pressure_basis_json.as_ref()];
        for (i, output) in outputs.iter().enumerate() {
            if let Some(path) = output {
                if path.is_empty() || inputs.iter().flatten().any(|input| input == path)
                    || outputs[..i].iter().flatten().any(|previous| previous == path) {
                    return Err("output paths must be distinct from inputs and each other".into());
                }
            }
        }
        let mut identities = std::collections::BTreeSet::new();
        for path in outputs.iter().flatten() {
            if !identities.insert(fresh_output_identity(path)?) {
                return Err("output paths resolve to the same file".into());
            }
        }
        Ok(options)
    }
    fn tuning_hz(&self) -> Option<f64> {
        self.concert_pitch.or_else(||
            (self.preset.is_some() && self.scale.is_none() && !self.raw_tensions).then_some(440.0))
    }
    fn uses_preset_board(&self) -> bool {
        self.preset.is_some() && self.board_geometry.is_none()
    }
}

// Resolve the actual parent before preparation; the leaf must not exist,
// including dangling symlinks. Final create_new admission still owns races.
fn fresh_output_identity(path: &str) -> Result<std::path::PathBuf, String> {
    if matches!(path.rsplit(std::path::is_separator).next(), None | Some("" | "." | "..")) {
        return Err(format!("{path}: output needs a regular-file leaf"));
    }
    let output = std::path::Path::new(path);
    match output.symlink_metadata() {
        Ok(_) => return Err(format!("{path}: output must be a fresh path")),
        Err(error) if error.kind() == std::io::ErrorKind::NotFound => {}
        Err(error) => return Err(format!("{path}: {error}")),
    }
    let name = output.file_name().ok_or_else(|| format!("{path}: output needs a filename"))?;
    let parent = output.parent().filter(|p| !p.as_os_str().is_empty())
        .unwrap_or_else(|| std::path::Path::new("."));
    let parent = parent.canonicalize().map_err(|e| format!("{path}: output parent: {e}"))?;
    if !parent.is_dir() {
        return Err(format!("{path}: output parent must be an existing directory"));
    }
    Ok(parent.join(name))
}
fn fresh_output(path: &str) -> Result<std::fs::File, String> {
    std::fs::OpenOptions::new().write(true).create_new(true).open(path)
        .map_err(|e| format!("{path}: fresh writable output required: {e}"))
}
fn write_fresh_output(path: &str, bytes: &[u8]) -> Result<(), String> {
    use std::io::Write;
    fresh_output(path)?.write_all(bytes).map_err(|e| format!("{path}: {e}"))
}

fn load_scale(text: Option<&str>) -> Result<Vec<geometry::Course>, String> {
    match text { Some(text) => geometry::read_scale(text), None => geometry::demonstration_scale() }
}
fn selected_scale(text: Option<&str>, options: &Options) -> Result<Vec<geometry::Course>, String> {
    let mut scale = if text.is_none() && options.preset.is_some() {
        steinway_scale::courses()?
    } else { load_scale(text)? };
    if let Some(reference) = options.tuning_hz() {
        for c in &mut scale {
            let target = reference * fs_math::det::pow(2.0, (f64::from(c.midi) - 69.0) / 12.0);
            let cents = 1200.0 * fs_math::det::ln(target / c.partial_hz(1, c.tension_n))
                / std::f64::consts::LN_2;
            c.tension_n = c.tension_at_cents(cents)
                .map_err(|e| format!("key {} tension retuning: {e}", c.midi))?;
        }
    }
    Ok(scale)
}
fn load_string_stretching(scale: &[geometry::Course], options: &Options)
    -> Result<Option<linear::string_stretching::Specification>, String> {
    options.string_stretching.as_deref().map(|path|
        linear::string_stretching::Specification::load(path, scale)).transpose()
}
fn prepare_instrument(scale: Vec<geometry::Course>, modes: &[linear::BoardMode],
    options: &Options) -> Result<engine::Instrument, String> {
    let stretching = load_string_stretching(&scale, options)?;
    prepare_instrument_with_string_material(scale, modes, options, stretching.as_ref())
}
fn prepare_instrument_with_string_material(scale: Vec<geometry::Course>, modes: &[linear::BoardMode],
    options: &Options, stretching: Option<&linear::string_stretching::Specification>)
    -> Result<engine::Instrument, String> {
    prepare_instrument_with_physical_controls(scale, modes, options, stretching, None)
}
fn prepare_instrument_with_physical_controls(scale: Vec<geometry::Course>, modes: &[linear::BoardMode],
    options: &Options, stretching: Option<&linear::string_stretching::Specification>,
    polarization: Option<&string_polarization::Prepared>) -> Result<engine::Instrument, String> {
    prepare_instrument_with_board_damping(scale, modes, options, stretching, polarization, None)
}
fn prepare_instrument_with_board_damping(scale: Vec<geometry::Course>, modes: &[linear::BoardMode],
    options: &Options, stretching: Option<&linear::string_stretching::Specification>,
    polarization: Option<&string_polarization::Prepared>, physical_damping: Option<&[f64]>)
    -> Result<engine::Instrument, String> {
    if options.board_reduction.is_some() && physical_damping.is_none() {
        return Err("reduced soundboard requires its full projected material damping; no diagonal fallback".into());
    }
    let dampers = match options.dampers.as_deref() {
        None => None,
        Some("estimated") => Some(linear::dampers::Specification::estimated(&scale)?),
        Some(path) => Some(linear::dampers::Specification::load(path, &scale)?),
    };
    let text = options.hammers.as_ref().map(|path| std::fs::read_to_string(path)
        .map_err(|e| format!("{path}: {e}"))).transpose()?;
    let mut piano = prepare_instrument_with_admitted_materials(scale, modes, options, text.as_deref(), stretching, polarization)?;
    if let Some(c) = physical_damping { piano.configure_bare_board_damping(c)?; }
    if let Some(spec) = &dampers { piano.configure_dampers(spec)?; }
    Ok(piano)
}
fn prepare_instrument_with_hammers(scale: Vec<geometry::Course>, modes: &[linear::BoardMode],
    options: &Options, text: Option<&str>) -> Result<engine::Instrument, String> {
    let stretching = load_string_stretching(&scale, options)?;
    prepare_instrument_with_admitted_materials(scale, modes, options, text, stretching.as_ref(), None)
}
fn prepare_instrument_with_admitted_materials(scale: Vec<geometry::Course>, modes: &[linear::BoardMode],
    options: &Options, text: Option<&str>, stretching: Option<&linear::string_stretching::Specification>,
    polarization: Option<&string_polarization::Prepared>)
    -> Result<engine::Instrument, String> {
    if options.string_stretching.is_some() && stretching.is_none() {
        return Err("supplied string stretching was not admitted; no linear fallback".into());
    }
    if options.string_polarization.is_some() && polarization.is_none() {
        return Err("supplied string polarization was not projected from the played board; no one-plane fallback".into());
    }
    let secondary = polarization.map(string_polarization::Prepared::secondary);
    let imported = text.map(|text| hammer_materials::read(text,
        &scale.iter().map(|c| c.midi).collect::<Vec<_>>())).transpose()?;
    let footprints=options.hammer_footprints.as_deref().map(|path|
        linear::hammer_footprint::Specification::load(path,&scale)).transpose()?;
    let source_rates=options.rt0425_hammer_dissipation.then(||
        scale.iter().map(|c|steinway_scale::hammer_relaxation_rt0425(c.midi))
            .collect::<Result<Vec<_>,_>>()).transpose()?;
    let mut piano = if options.preset.is_some() {
        let materials = match imported {
            Some(materials) => materials,
            None => scale.iter().map(if options.rt0425_hammer_dissipation {
                steinway_scale::hammer_material_rt0425_damped
            } else if options.rt0425_hammer_stiffness {
                steinway_scale::hammer_material_rt0425
            } else { steinway_scale::hammer_material }).collect::<Result<Vec<_>,_>>()?,
        };
        engine::Instrument::new_with_string_damping(scale, modes, options.sample_rate,
            options.substeps, options.modes, true, materials, Some(engine::ShankGeometry::published()),
            footprints.as_ref(), secondary, options.rt0425_string_damping)
    } else if let Some(materials) = imported {
        engine::Instrument::new_with_transverse_contact_geometry(scale, modes, options.sample_rate,
            options.substeps, options.modes, true, materials, None, footprints.as_ref(), secondary)
    } else if secondary.is_some() {
        engine::Instrument::new_with_demonstration_geometry(scale, modes, options.sample_rate,
            options.substeps, options.modes, true, footprints.as_ref(), secondary)
    } else if let Some(spec)=&footprints {
        engine::Instrument::new_with_footprints(scale,modes,options.sample_rate,options.substeps,
            options.modes,true,spec)
    } else {
        engine::Instrument::new(scale, modes, options.sample_rate, options.substeps, options.modes, true)
    }?;
    if let Some(rates)=source_rates {piano.configure_source_hammer_dissipation(&rates)?;}
    if let Some(spec) = stretching { piano.configure_string_stretching(spec)?; }
    Ok(piano)
}
fn load_board(text: Option<&str>, scale: &[geometry::Course]) -> Result<Vec<linear::BoardMode>, String> {
    match text {
        Some(text) => board::read(text, &scale.iter().map(|c| c.midi).collect::<Vec<_>>()),
        None => Ok(board::demonstration()),
    }
}
fn prepare_geometric_board(text: &str, keys: &[u8], band_hz: f64,
    equilibrate_mass: bool, consistent_mass: bool, edge_cubic_mass: bool, acoustic_refinement_levels: usize)
    -> Result<board_geometry::PreparedBoard, String> {
    prepare_geometric_board_motion(text, keys, band_hz, equilibrate_mass, consistent_mass,
        edge_cubic_mass, acoustic_refinement_levels, false)
}
#[allow(clippy::too_many_arguments)]
fn prepare_geometric_board_motion(text: &str, keys: &[u8], band_hz: f64,
    equilibrate_mass: bool, consistent_mass: bool, edge_cubic_mass: bool,
    acoustic_refinement_levels: usize, retain_motion: bool) -> Result<board_geometry::PreparedBoard, String> {
    prepare_geometric_board_with_reduction(text, keys, band_hz, equilibrate_mass,
        consistent_mass, edge_cubic_mass, acoustic_refinement_levels, retain_motion, None)
}
#[allow(clippy::too_many_arguments)]
fn prepare_geometric_board_with_reduction(text: &str, keys: &[u8], band_hz: f64,
    equilibrate_mass: bool, consistent_mass: bool, edge_cubic_mass: bool,
    acoustic_refinement_levels: usize, retain_motion: bool,
    reduction: Option<&board_geometry::ritz::RitzOptions>) -> Result<board_geometry::PreparedBoard, String> {
    if let Some(reduction) = reduction {
        if crowned_board::is_crowned(text) {
            if equilibrate_mass || consistent_mass || edge_cubic_mass || acoustic_refinement_levels != 0 {
                return Err("flat-board mass controls require a flat geometric board".into());
            }
            return crowned_board::CrownedBoard::read(text)?
                .prepare_reduced(keys, band_hz, retain_motion, reduction);
        }
        return board_geometry::BoardGeometry::read(text)?
            .with_acoustic_refinement(acoustic_refinement_levels)?
            .prepare_reduced(keys, band_hz, retain_motion, equilibrate_mass,
                consistent_mass, edge_cubic_mass, reduction);
    }
    if crowned_board::is_crowned(text) {
        if equilibrate_mass || consistent_mass || edge_cubic_mass || acoustic_refinement_levels != 0 {
            return Err("flat-board mass controls require a flat geometric board".into());
        }
        let geometry = crowned_board::CrownedBoard::read(text)?;
        if retain_motion { geometry.prepare_with_motion(keys, band_hz) }
        else { geometry.prepare(keys, band_hz) }
    } else {
        let geometry=board_geometry::BoardGeometry::read(text)?
            .with_acoustic_refinement(acoustic_refinement_levels)?;
        if retain_motion && edge_cubic_mass { geometry.prepare_with_motion_edge_cubic_transverse_mass(keys,band_hz,equilibrate_mass) }
        else if retain_motion && consistent_mass { geometry.prepare_with_motion_consistent_transverse_mass(keys,band_hz,equilibrate_mass) }
        else if retain_motion && equilibrate_mass { geometry.prepare_with_motion_mass_equilibrated(keys,band_hz) }
        else if retain_motion { geometry.prepare_with_motion(keys,band_hz) }
        else if edge_cubic_mass { geometry.prepare_edge_cubic_transverse_mass(keys,band_hz,equilibrate_mass) }
        else if consistent_mass { geometry.prepare_consistent_transverse_mass(keys,band_hz,equilibrate_mass) }
        else if equilibrate_mass { geometry.prepare_mass_equilibrated(keys,band_hz) }
        else { geometry.prepare(keys,band_hz) }
    }
}
/// Export only the admitted keys. An absent measurement must not become an
/// apparently measured zero bridge coefficient when the table is re-imported.
fn write_board_for_scale(modes: &[linear::BoardMode], scale: &[geometry::Course]) -> String {
    let all = board::write(modes);
    let mut out = String::new();
    for row in all.lines() {
        let admitted = row == board::HEADER || row.split(',').nth(4)
            .and_then(|key| key.parse::<u8>().ok())
            .is_some_and(|key| scale.iter().any(|c| c.midi == key));
        if admitted { out.push_str(row); out.push('\n'); }
    }
    out
}
fn study_key(scale: &[geometry::Course], requested: Option<u8>) -> Result<u8, String> {
    if let Some(key) = requested {
        return scale.iter().any(|c| c.midi == key).then_some(key)
            .ok_or_else(|| format!("requested key {key} is absent from the input scale"));
    }
    scale.iter().min_by_key(|c| c.midi.abs_diff(69)).map(|c| c.midi)
        .ok_or_else(|| "empty input scale".into())
}

fn render(path: &str, scale: Vec<geometry::Course>, modes: &[linear::BoardMode],
    surface: Option<&[board_geometry::SurfaceSample]>, options: &Options) -> Result<(), String> {
    let stretching = load_string_stretching(&scale, options)?;
    render_with_string_material(path, scale, modes, surface, options, stretching.as_ref())
}
fn check_pcm_headroom(peak_pa: f64, full_scale_pa: f64) -> Result<(), String> {
    if !peak_pa.is_finite() {
        return Err("non-finite rendered pressure cannot be encoded".into());
    }
    if peak_pa >= full_scale_pa {
        return Err(format!("rendered peak {peak_pa:.6} Pa reaches PCM full scale {full_scale_pa:.6} Pa; rerun with --pcm-full-scale-pa greater than {peak_pa:.6} to avoid clipping"));
    }
    Ok(())
}
fn bridge_acceleration(velocity: &[f64], sample: usize, rate: u32) -> f64 {
    let left = sample.saturating_sub(1);
    let right = (sample + 1).min(velocity.len() - 1);
    (velocity[right] - velocity[left]) * f64::from(rate) / (right - left) as f64
}
fn write_bridge_trace(path: &str, velocity: &[f64], pressure: &[f64],
    rate: u32, channels: usize) -> Result<(), String> {
    use std::io::Write;
    let file = fresh_output(path)?;
    let mut out = std::io::BufWriter::new(file);
    writeln!(out, "sample,time_s,bridge_velocity_m_s,bridge_acceleration_m_s2,pressure_left_pa")
        .map_err(|e| format!("{path}: {e}"))?;
    for (i, &v) in velocity.iter().enumerate() {
        writeln!(out, "{i},{:.17e},{:.17e},{:.17e},{:.17e}",
            i as f64 / f64::from(rate), v, bridge_acceleration(velocity, i, rate),
            pressure[i * channels]).map_err(|e| format!("{path}: {e}"))?;
    }
    out.flush().map_err(|e| format!("{path}: {e}"))
}
fn write_modal_pressure(path: &str, modal: &[f64], pressure: &[f64],
    rate: u32, channels: usize, modes: usize) -> Result<(), String> {
    use std::io::Write;
    let file = fresh_output(path)?;
    let mut out = std::io::BufWriter::new(file);
    write!(out, "sample,time_s").map_err(|e| format!("{path}: {e}"))?;
    for channel in 0..channels {
        let side = if channel == 0 {"left"} else {"right"};
        write!(out, ",pressure_{side}_pa").map_err(|e| format!("{path}: {e}"))?;
        for mode in 0..modes {
            write!(out, ",{side}_mode_{mode}_pa").map_err(|e| format!("{path}: {e}"))?;
        }
    }
    writeln!(out).map_err(|e| format!("{path}: {e}"))?;
    for (sample, frame) in pressure.chunks_exact(channels).enumerate() {
        write!(out, "{sample},{:.17e}", sample as f64 / f64::from(rate))
            .map_err(|e| format!("{path}: {e}"))?;
        for channel in 0..channels {
            write!(out, ",{:.17e}", frame[channel]).map_err(|e| format!("{path}: {e}"))?;
            let start = (sample * channels + channel) * modes;
            let components = &modal[start..start + modes];
            let reconstructed: f64 = components.iter().sum();
            if (reconstructed - frame[channel]).abs() > 1e-10 * frame[channel].abs().max(1.0) {
                return Err(format!("modal pressure fails to reconstruct channel {channel} sample {sample}"));
            }
            for value in components {
                write!(out, ",{value:.17e}").map_err(|e| format!("{path}: {e}"))?;
            }
        }
        writeln!(out).map_err(|e| format!("{path}: {e}"))?;
    }
    out.flush().map_err(|e| format!("{path}: {e}"))
}
fn write_receiver_pressure(out: &mut impl std::io::Write, pressure: &[f64],
    rate: u32, channels: usize) -> Result<(), String> {
    if !(1..=2).contains(&channels) || rate == 0 || !pressure.len().is_multiple_of(channels)
        || pressure.iter().any(|p|!p.is_finite()) {
        return Err("receiver pressure needs complete finite mono/stereo frames and a positive clock".into());
    }
    write!(out,"sample,time_s,pressure_left_pa").map_err(|e|e.to_string())?;
    if channels==2 {write!(out,",pressure_right_pa").map_err(|e|e.to_string())?;}
    writeln!(out).map_err(|e|e.to_string())?;
    for (sample,frame) in pressure.chunks_exact(channels).enumerate() {
        write!(out,"{sample},{:.17e}",sample as f64/f64::from(rate)).map_err(|e|e.to_string())?;
        for value in frame {write!(out,",{value:.17e}").map_err(|e|e.to_string())?;}
        writeln!(out).map_err(|e|e.to_string())?;
    }
    out.flush().map_err(|e|e.to_string())
}
fn render_with_string_material(path: &str, scale: Vec<geometry::Course>, modes: &[linear::BoardMode],
    surface: Option<&[board_geometry::SurfaceSample]>, options: &Options,
    stretching: Option<&linear::string_stretching::Specification>) -> Result<(), String> {
    render_with_physical_controls(path, scale, modes, surface, options, stretching, None)
}
fn render_with_physical_controls(path: &str, scale: Vec<geometry::Course>, modes: &[linear::BoardMode],
    surface: Option<&[board_geometry::SurfaceSample]>, options: &Options,
    stretching: Option<&linear::string_stretching::Specification>,
    polarization: Option<&string_polarization::Prepared>) -> Result<(), String> {
    render_with_board_damping(path, scale, modes, surface, options, stretching, polarization, None)
}
#[allow(clippy::too_many_arguments)]
fn render_with_board_damping(path: &str, scale: Vec<geometry::Course>, modes: &[linear::BoardMode],
    surface: Option<&[board_geometry::SurfaceSample]>, options: &Options,
    stretching: Option<&linear::string_stretching::Specification>,
    polarization: Option<&string_polarization::Prepared>, physical_damping: Option<&[f64]>)
    -> Result<(), String> {
    let keys: Vec<u8> = scale.iter().map(|c| c.midi).collect();
    let rate = options.sample_rate;
    let count = (options.duration * f64::from(rate)).round() as u32;
    let score = if let Some(path) = &options.midi {
        let parsed = midi::load(path, &keys, rate, u64::from(count), options.midi_mapping)?;
        let report = &parsed.report;
        println!("MIDI: {} tracks, channel {}, {} hammer launches, end sample {}, {} end releases; {} other-channel messages, {} unsupported channel messages and {} SysEx events ignored.",
            report.tracks, options.midi_mapping.channel + 1, report.selected_note_ons,
            report.end_sample, report.end_releases, report.other_channel_messages,
            report.ignored_channel_messages, report.ignored_sysex_events);
        println!("Uncalibrated MIDI mapping: velocity 127 -> {} m/s, linear hammer launch; sustain mapping {}. No output gain, sample bank or pitch-wheel substitution.",
            options.midi_mapping.maximum_velocity_m_s,
            if options.midi_mapping.continuous_sustain { "CC64/127 travel" } else { "switch at 64" });
        performance::Performance::from_events(parsed.events, &keys, u64::from(count))?
    } else { match &options.performance {
        Some(path) => performance::Performance::read(
            &std::fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?,
            &keys, u64::from(count))?,
        None => performance::Performance::demonstration(&keys, rate, u64::from(count),
            options.note, options.velocity)?,
    }};
    let key = if options.performance.is_some()
        && (options.bridge_trace_csv.is_some() || options.modal_pressure_csv.is_some()) {
        score.single_excitation_key()?
    } else { study_key(&scale, options.note)? };
    let observed_course = options.bridge_trace_csv.as_ref().map(|_| {
        scale.iter().position(|c| c.midi == key)
            .ok_or_else(|| format!("bridge observation key {key} is absent"))
    }).transpose()?;
    let piano = prepare_instrument_with_board_damping(scale, modes, options, stretching, polarization, physical_damping)?;
    let bridge_row = observed_course.map(|course| {
        piano.bank.strings.iter().find(|s|
            s.course == course && s.member == 0 && s.polarization == 0 && !s.duplex)
            .map(|s| s.bridge.clone())
            .ok_or_else(|| format!("no speaking vertical bridge port for key {key}"))
    }).transpose()?;
    debug_assert_eq!(piano.sample_rate(), rate);
    if let Some(polarization) = polarization {
        println!("Two transverse directions per string; bridge motion and lateral damper ratios supplied by {}. One hammer/contact area and shared mechanics clock; see STRING_POLARIZATION.md.", polarization.source);
    }
    if let Some(path) = &options.string_stretching {
        let count = (0..piano.bank.strings.len())
            .filter(|&i| piano.bank.string_stretching_observation(i).is_some()).count();
        println!("String extension: {path}; {count} geometric speaking/duplex channels; same-tick contact and reciprocal bridge work. Zero means explicit all-linear selection. EA is supplied, not inferred; slope/iteration limits refuse rather than clamp. No real-time or calibration claim; see STRING_STRETCHING.md.");
    }
    if let Some(path)=&options.hammer_footprints {
        println!("Hammer contact geometry: {path}; {} independent felt sites, unchanged total course area and hammer mass. No inferred/calibrated face width; see HAMMER_FOOTPRINTS.md.",
            piano.hammer_contact_count());
    }
    if let Some((strings,cells)) = piano.damper_resolution() {
        println!("Spatial viscous dampers: {} speaking-string pads, {} quadrature stations; source {}. No measured pad/action or real-time claim; see DAMPERS.md.",
            strings, cells, options.dampers.as_deref().unwrap_or("supplied"));
    }
    let surface = if options.diagnostic_volume { None } else { surface };
    let left=options.microphone.unwrap_or([0.675,1.0,1.0]);
    let mut stream = match options.microphone_right {
        Some(right)=>audio::AudioStream::new_stereo(piano,score,
            surface.ok_or("stereo pressure requires a geometric soundboard surface")?,
            [left,right],fs_bem::helmholtz::Medium::air())?,
        None=>audio::AudioStream::new(piano,score,surface,left,
            fs_bem::helmholtz::Medium::air(),options.observer_gain)?,
    };
    if let Some(mic) = stream.microphone() {
        println!("Rayleigh pressure microphone at {:?} m; propagation {:?} samples plus {} anti-alias delay samples; {} modal multiplies/output sample.",
            mic.position_m, mic.delay_samples, mic.filter_delay_samples(), mic.multiply_adds_per_sample());
    }
    if let Some([_,right])=stream.stereo_microphones() {
        println!("Right Rayleigh microphone at {:?} m; propagation {:?} samples plus {} anti-alias delay samples. Independent spatial pressure; shared mechanics.",
            right.position_m,right.delay_samples,right.filter_delay_samples());
    }
    let channels=stream.channels();
    let mut pressure = vec![0.0; count as usize*channels];
    let mut bridge_velocity = bridge_row.as_ref().map(|_| Vec::with_capacity(count as usize));
    let modal_modes = stream.instrument().bank.board_count;
    let mut modal_pressure = if options.modal_pressure_csv.is_some() {
        let values = (count as usize).checked_mul(channels)
            .and_then(|n| n.checked_mul(modal_modes))
            .filter(|&n| n <= 4_000_000)
            .ok_or("modal pressure trace exceeds the 4-million-value diagnostic budget")?;
        Some(vec![0.0; values])
    } else { None };
    let mut trace_frame = 0;
    let start = std::time::Instant::now();
    for block in pressure.chunks_mut(256*channels) {
        if bridge_velocity.is_some() || modal_pressure.is_some() {
            for frame in block.chunks_mut(channels) {
                stream.render_interleaved_block(frame).map_err(|e| e.to_string())?;
                if let (Some(row), Some(velocity)) = (&bridge_row, &mut bridge_velocity) {
                    let bank = &stream.instrument().bank;
                    velocity.push(row.iter().zip(&bank.v[bank.modes.len()..])
                        .map(|(shape, speed)| shape * speed).sum::<f64>());
                }
                if let Some(modal) = &mut modal_pressure {
                    let offset = trace_frame * channels * modal_modes;
                    if let Some([left, right]) = stream.stereo_microphones() {
                        left.mode_pressures(&mut modal[offset..offset + modal_modes])?;
                        right.mode_pressures(&mut modal[offset + modal_modes..offset + 2*modal_modes])?;
                    } else {
                        stream.microphone().ok_or("modal trace needs a physical microphone")?
                            .mode_pressures(&mut modal[offset..offset + modal_modes])?;
                    }
                }
                trace_frame += 1;
            }
        } else {
            stream.render_interleaved_block(block).map_err(|e| e.to_string())?;
        }
    }
    let elapsed = start.elapsed().as_secs_f64();
    let peak = pressure.iter().fold(0.0_f64, |a, p| a.max(p.abs()));
    let seconds = f64::from(count) / f64::from(rate);
    check_pcm_headroom(peak, options.pcm_full_scale_pa)?;
    let (wav, clips) = fs_couple::pcm_wav::encode_pcm16_wav_interleaved(&pressure, rate, channels as u16, options.pcm_full_scale_pa).map_err(|e| e.to_string())?;
    write_fresh_output(path, &wav)?;
    if let (Some(csv), Some(velocity)) = (&options.bridge_trace_csv, &bridge_velocity) {
        write_bridge_trace(csv, velocity, &pressure, rate, channels)?;
        println!("Modeled vertical bridge motion at key {key}: {csv}; centered output-rate acceleration, left pressure, same sample indices. Pressure includes propagation and anti-alias delay.");
    }
    if let (Some(csv), Some(modal)) = (&options.modal_pressure_csv, &modal_pressure) {
        write_modal_pressure(csv, modal, &pressure, rate, channels, modal_modes)?;
        println!("Loaded-board modal pressure at key {key}: {csv}; signed contributions at each physical receiver on the WAV clock. Basis indices are not bare-board eigenfrequencies.");
    }
    if let Some(csv)=&options.receiver_pressure_csv {
        let file=fresh_output(csv)?;
        write_receiver_pressure(&mut std::io::BufWriter::new(file),&pressure,rate,channels)
            .map_err(|e|format!("{csv}: {e}"))?;
        println!("Total receiver pressure: {csv}; unquantized Pa, same frames and channels as WAV, no modal trace allocation.");
    }
    if let Some(json) = &options.pressure_basis_json {
        let file = fresh_output(json)?;
        pressure_basis::write(&mut std::io::BufWriter::new(file), &stream.instrument().bank)
            .map_err(|e| format!("{json}: {e}"))?;
        println!("Loaded pressure basis: {json}; exact projection map and diagonal reference frequencies, not full instrument poles.");
    }
    if stream.microphone().is_some() {
        println!("Computed half-space pressure in Pa; PCM full scale {} Pa, no peak normalization. Infinite baffle; no room/lid scattering, radiation loading or measured-SPL calibration.", options.pcm_full_scale_pa);
    } else {
        println!("Diagnostic volume-velocity observer, gain {} Pa/(m^3/s); no peak normalization or calibrated SPL claim.", options.observer_gain);
    }
    let piano = stream.instrument();
    println!("{} string modes; {} board modes; {} above-band duplex segments omitted from dynamic retention (static attachment retained).",
        piano.bank.modes.len(), piano.bank.board_count, piano.bank.omitted_duplex_modes);
    println!("{seconds:.6} s, {channels} channel(s) rendered in {elapsed:.6} s; wall/audio ratio {:.4}; peak {peak:.6} Pa-equivalent; {clips} PCM clips.", elapsed / seconds);
    println!("Input {:.9} J; stored {:.9} J; component losses {:.9} J; closure {:.3e} J; worst substep defect {:.3e} J.",
        piano.accounting.input_work_j, piano.energy_j(), piano.accounting.dissipated_j(),
        piano.accounting.input_work_j - piano.energy_j() - piano.accounting.dissipated_j(), piano.accounting.max_balance_error_j);
    println!("Felt loss {:.9} J, including {:.9} J time-dependent relaxation; shank damping {:.9} J.",
        piano.accounting.felt_loss_j, piano.accounting.felt_relaxation_loss_j, piano.accounting.shank_loss_j);
    if piano.damper_resolution().is_some() {
        println!("Spatial damper loss {:.9} J, already included in component losses; no output-envelope damping.",
            piano.accounting.damper_loss_j);
    }
    Ok(())
}

fn run() -> Result<(), String> {
    let options = Options::parse(&std::env::args().skip(1).collect::<Vec<_>>())?;
    if options.help { println!("{USAGE}"); return Ok(()); }
    let read = |path: &String| std::fs::read_to_string(path).map_err(|e| format!("{path}: {e}"));
    let scale_text = options.scale.as_ref().map(read).transpose()?;
    let board_text = options.board.as_ref().map(read).transpose()?;
    let scale = selected_scale(scale_text.as_deref(), &options)?;
    // Freeze supplied string material before geometry work or any exports.
    // Playback receives this admitted value, never a later read of the path.
    let stretching = load_string_stretching(&scale, &options)?;
    let polarization = options.string_polarization.as_deref().map(|path|
        string_polarization::Specification::load(path, &scale)).transpose()?;
    let scale_source = options.scale.as_deref().unwrap_or(if options.preset.is_some() {
        "RT-0425 Appendix A wrapped-string MODEL: 84 published courses plus four estimated extensions"
    } else { "ESTIMATED demonstration" });
    let tuning_source = options.tuning_hz().map_or_else(|| "input tensions preserved".to_owned(),
        |f| format!("tensions adjusted to A4={f} Hz first-partial equal temperament; L, mass and EI preserved"));
    let preset = options.uses_preset_board().then(|| if options.rt0425_bridge_contacts {
        steinway_d::build_with_rt0425_contacts(options.mesh_divisions)
    } else { steinway_d::build(options.mesh_divisions) }).transpose()?;
    if let Some(preset) = &preset {
        if let Some(path) = &options.dump_geometry { write_fresh_output(path, preset.geometry.as_bytes())?; }
        if let Some(path) = &options.dump_obj { write_fresh_output(path, preset.obj.as_bytes())?; }
        if options.render.is_none() && options.dump_board.is_none() && options.dump_scale.is_none()
            && (options.dump_geometry.is_some() || options.dump_obj.is_some()) {
            println!("Exported source-derived Model D: tapered panel, 17 ribs, maple bridges, cut-off bar and 88 bridge stations. No eigenanalysis was needed.");
            return Ok(());
        }
    }
    let geometry_text = match (&preset, &options.board_geometry) {
        (_, Some(path)) => Some(read(path)?),
        (Some(p), _) => Some(p.geometry.clone()),
        _ => None,
    };
    let (modes, board_source, surface, polarization, physical_damping) = if let Some(text) = &geometry_text {
        let start = std::time::Instant::now();
        let prepared = prepare_geometric_board_with_reduction(text, &scale.iter().map(|c| c.midi).collect::<Vec<_>>(),
            options.board_band_hz, options.equilibrate_board_mass, options.consistent_board_mass,
            options.edge_cubic_board_mass, options.acoustic_refinement_levels, polarization.is_some(),
            options.board_reduction.as_ref())?;
        let projected = polarization.as_ref().map(|spec|
            spec.project(&scale, &prepared.modes, prepared.motion.as_ref())).transpose()?;
        let model_name = if crowned_board::is_crowned(text) { "Crowned shell" } else { "Flat plate" };
        println!("{model_name} {:.6} m^2, {:.6} kg (panel+ribs/bridges), {} free DOFs, {} retained modes; source band (0,{}] Hz; preparation {:.6} s.",
            prepared.area_m2, prepared.mass_kg, prepared.free_dofs, prepared.modes.len(),
            options.board_band_hz, start.elapsed().as_secs_f64());
        if let Some(report) = &prepared.reduction {
            println!("Bridge-informed reduction: {} certified source modes -> {} retained coordinates; {} exact low modes; {} normalized static/harmonic snapshots at {:?} Hz; maximum displacement projection residual {:.9e}. This residual does not certify transfer, acoustic or mesh convergence.",
                report.source_modes, prepared.modes.len(), report.protected_low_modes,
                report.snapshot_count, report.sample_hz, report.max_relative_snapshot_error);
            for (i, interval) in report.source_frequency_intervals_hz.iter().enumerate() {
                println!("source FE mode {i}: [{:.9}, {:.9}] Hz", interval.0, interval.1);
            }
        }
        for (i, interval) in prepared.frequency_intervals_hz.iter().enumerate() {
            let scope = if prepared.reduction.as_ref().is_some_and(|r| i >= r.protected_low_modes) {
                "projected-pencil Ritz"
            } else { "source FE" };
            println!("board mode {i} ({scope}): [{:.9}, {:.9}] Hz", interval.0, interval.1);
        }
        (prepared.modes, format!("GEOMETRY-DERIVED {model_name}; {}; rim compliance not modeled", prepared.provenance),
            Some(prepared.surface), projected, prepared.physical_damping)
    } else {
        (load_board(board_text.as_deref(), &scale)?, options.board.as_deref()
            .unwrap_or("AUTHORED illustrative modes; not measured Steinway geometry").to_owned(), None, None, None)
    };
    if modes.iter().any(|m| m.frequency_hz >= 0.45 * f64::from(options.sample_rate)) {
        return Err("soundboard mode at/above output retention ceiling; use an explicitly reduced board".into());
    }
    study_key(&scale, options.note)?;
    println!("String scale: {scale_source}; {tuning_source}.");
    println!("Soundboard: {board_source}.");
    if let Some(path) = &options.hammers {
        println!("Hammer materials supplied by {path}; no automatic coupon-fit or calibration claim. Scale mass/patch geometry unchanged.");
    } else if options.preset.is_some() {
        println!("Per-key source-derived hammer loading envelopes{}; {}.",
            if options.rt0425_hammer_stiffness { " with RT-0425 stiffness per unison string" } else { " with legacy stiffness divided across the unison" },
            if options.rt0425_hammer_dissipation { "RT-0425 R_H power-law rate loss, no estimated crush or Prony terms" }
            else { "estimated unloading/crush and tangent-scaled Prony relaxation, not RT-0425 R_H" });
    }
    if options.preset.is_some() {
        println!("Published shank geometry -> rigid rotation + bending; reciprocal jack port and 1.5 mm let-off. Linearized action fragment; damping/backcheck estimated.");
        println!("String damping: {}.", if options.rt0425_string_damping {
            "RT-0425 per-key R_u and eta_u projected onto scalar stiff-string modes; four end keys extrapolated"
        } else { "estimated common rate and authored Maxwell bending spectrum" });
    }
    println!("Source authority belongs to the inputs, not the model name; imported files are not independently certified measurements.");
    if let Some(path) = &options.dump_scale {
        write_fresh_output(path, format!("# Source: {scale_source}; {tuning_source}\n{}", geometry::write_scale(&scale)).as_bytes())?;
    }
    if let Some(path) = &options.dump_board {
        write_fresh_output(path, format!("# Source: {board_source}\n{}", write_board_for_scale(&modes, &scale)).as_bytes())?;
    }
    if let Some(path) = &options.render {
        return render_with_board_damping(path, scale, &modes, surface.as_deref(), &options,
            stretching.as_ref(), polarization.as_ref(), physical_damping.as_deref());
    }
    println!("Model D published envelope: {} x {} m; board {} -> {} m (center -> edge).",
        geometry::D_LENGTH_M, geometry::D_WIDTH_M, geometry::D_BOARD_CENTER_M, geometry::D_BOARD_EDGE_M);
    println!("{} courses; {} speaking strings; {} board modes.", scale.len(), scale.iter().map(|c| c.unison).sum::<usize>(), modes.len());
    for c in scale { println!("key {}: L={:.4} m, mu={:.6} kg/m, T={:.2} N, f1={:.3} Hz", c.midi, c.length_m, c.linear_density_kg_m, c.tension_n, c.partial_hz(1, c.tension_n)); }
    Ok(())
}
fn main() {
    if let Err(error) = run() { eprintln!("grand_piano: {error}"); std::process::exit(1); }
}

#[cfg(test)]
#[path = "board_reduction_render_tests.rs"]
mod board_reduction_render_tests;

#[cfg(test)]
#[path = "string_polarization_render_tests.rs"]
mod string_polarization_render_tests;

#[cfg(test)]
#[path = "string_stretching_render_tests.rs"]
mod string_stretching_render_tests;

#[cfg(test)]
#[path = "hammer_footprint_render_tests.rs"]
mod footprint_render_tests;

#[cfg(test)]
#[path = "hammer_profile_render_tests.rs"]
mod hammer_profile_render_tests;

#[cfg(test)]
#[path = "crowned_render_tests.rs"]
mod crowned_render_tests;

#[cfg(test)]
mod render_tests {
    use super::*;
    fn options(args: &[&str]) -> Result<Options, String> {
        Options::parse(&args.iter().map(|s| (*s).to_owned()).collect::<Vec<_>>())
    }
    fn output_fixture() -> std::path::PathBuf {
        let stamp = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH)
            .unwrap().as_nanos();
        let dir = std::env::temp_dir().join(format!("fs-piano-outputs-{}-{stamp}", std::process::id()));
        std::fs::create_dir(&dir).unwrap();
        std::fs::create_dir(dir.join("child")).unwrap();
        dir
    }
    fn export_options(wav: &std::path::Path, csv: &std::path::Path) -> Result<Options, String> {
        options(&["--preset", "steinway-d", "--render", wav.to_str().unwrap(),
            "--receiver-pressure-csv", csv.to_str().unwrap()])
    }
    #[test]
    fn pressure_basis_export_requires_physical_render_and_fresh_distinct_paths() {
        let dir = output_fixture(); let wav = dir.join("piano.wav");
        let basis = dir.join("basis.json");
        let common = ["--preset", "steinway-d", "--render", wav.to_str().unwrap()];
        let args: Vec<_> = common.into_iter().chain(["--pressure-basis-json",
            basis.to_str().unwrap()]).collect();
        let accepted = options(&args).unwrap();
        assert_eq!(accepted.pressure_basis_json.as_deref(), basis.to_str());
        for invalid in [vec!["--pressure-basis-json", basis.to_str().unwrap()],
            vec!["--render", wav.to_str().unwrap(), "--pressure-basis-json", basis.to_str().unwrap()],
            vec!["--preset", "steinway-d", "--render", wav.to_str().unwrap(),
                "--diagnostic-volume", "--pressure-basis-json", basis.to_str().unwrap()],
            vec!["--preset", "steinway-d", "--render", wav.to_str().unwrap(),
                "--pressure-basis-json", wav.to_str().unwrap()]] {
            assert!(options(&invalid).is_err());
        }
        let alias = dir.join("child/../piano.wav");
        assert!(options(&common.into_iter().chain(["--pressure-basis-json",
            alias.to_str().unwrap()]).collect::<Vec<_>>()).is_err());
        write_fresh_output(basis.to_str().unwrap(), b"existing basis").unwrap();
        assert!(options(&args).is_err());
        assert_eq!(std::fs::read(basis).unwrap(), b"existing basis");
        assert!(!wav.exists());
    }
    #[test]
    fn export_paths_refuse_existing_entries_and_dot_aliases_before_preparation() {
        let dir = output_fixture(); let wav = dir.join("fresh.wav");
        assert!(export_options(&wav, &dir.join("fresh.csv")).is_ok());
        for alias in [dir.join("./fresh.wav"), dir.join("child/../fresh.wav")] {
            assert!(export_options(&wav, &alias).unwrap_err().contains("same file"));
        }
        let source = dir.join("events.csv");
        write_fresh_output(source.to_str().unwrap(), b"source input remains intact").unwrap();
        let args = ["--preset", "steinway-d", "--render", wav.to_str().unwrap(),
            "--performance", source.to_str().unwrap(), "--receiver-pressure-csv",
            dir.join("./events.csv").to_str().unwrap()].map(str::to_owned);
        assert!(Options::parse(&args).unwrap_err().contains("fresh path"));
        assert!(export_options(&wav, &dir.join("child")).is_err());
        assert!(export_options(&wav, &dir.join("missing/pressure.csv")).is_err());
        assert!(export_options(&wav, &source.join("pressure.csv")).is_err());
        let trailing = format!("{}{}", dir.join("absent.csv").display(), std::path::MAIN_SEPARATOR);
        assert!(export_options(&wav, std::path::Path::new(&trailing)).is_err());
        for leaf in [format!("{trailing}."), format!("{trailing}..") ] {
            assert!(export_options(&wav, std::path::Path::new(&leaf)).is_err());
        }
        assert_eq!(std::fs::read(source).unwrap(), b"source input remains intact");
        assert!(!wav.exists());
        // Retain the test's newly created files; never overwrite or delete input.
    }
    #[test]
    fn publication_preserves_an_entry_created_after_option_admission() {
        let dir = output_fixture(); let wav = dir.join("fresh.wav");
        assert!(export_options(&wav, &dir.join("fresh.csv")).is_ok());
        write_fresh_output(wav.to_str().unwrap(), b"first writer").unwrap();
        assert!(write_fresh_output(wav.to_str().unwrap(), b"second writer").is_err());
        assert_eq!(std::fs::read(&wav).unwrap(), b"first writer");
        assert!(fresh_output(wav.to_str().unwrap()).is_err());
    }
    #[cfg(unix)]
    #[test]
    fn export_paths_refuse_symlinked_parents_dangling_links_and_hard_links() {
        let dir = output_fixture(); let child = dir.join("child"); let link = dir.join("link");
        std::os::unix::fs::symlink(&child, &link).unwrap();
        assert!(export_options(&child.join("same.wav"), &link.join("same.wav"))
            .unwrap_err().contains("same file"));
        let dangling = dir.join("dangling.csv");
        std::os::unix::fs::symlink(dir.join("absent.csv"), &dangling).unwrap();
        assert!(export_options(&dir.join("fresh.wav"), &dangling).is_err());
        let source = dir.join("events.csv");
        write_fresh_output(source.to_str().unwrap(), b"hard-linked source").unwrap();
        let alias = dir.join("hard.csv"); std::fs::hard_link(&source, &alias).unwrap();
        assert!(export_options(&dir.join("fresh.wav"), &alias).is_err());
        assert_eq!(std::fs::read(source).unwrap(), b"hard-linked source");
        assert!(dangling.symlink_metadata().is_ok());
        assert!(!dir.join("absent.csv").exists());
    }
    #[test]
    fn rt0425_string_damping_requires_the_preset_scale_and_a_render() {
        assert!(options(&["--preset", "steinway-d", "--render", "a.wav",
            "--rt0425-string-damping"]).unwrap().rt0425_string_damping);
        assert!(options(&["--preset", "steinway-d", "--rt0425-string-damping"]).is_err());
        assert!(options(&["--render", "a.wav", "--rt0425-string-damping"]).is_err());
        assert!(options(&["--preset", "steinway-d", "--scale", "custom.csv",
            "--render", "a.wav", "--rt0425-string-damping"]).is_err());
    }
    #[test]
    fn pcm_pressure_scale_is_explicit_and_overrange_refuses_before_write() {
        assert_eq!(options(&[]).unwrap().pcm_full_scale_pa, 2.0);
        let o = options(&["--render", "piano.wav", "--pcm-full-scale-pa", "4"]).unwrap();
        assert_eq!(o.pcm_full_scale_pa, 4.0);
        assert!(check_pcm_headroom(3.193, o.pcm_full_scale_pa).is_ok());
        let refusal = check_pcm_headroom(3.193, 2.0).unwrap_err();
        assert!(refusal.contains("--pcm-full-scale-pa greater than 3.193000"));
        assert!(check_pcm_headroom(2.0, 2.0).is_err());
        assert!(check_pcm_headroom(f64::NAN, 4.0).is_err());
        for args in [vec!["--pcm-full-scale-pa", "4"],
            vec!["--render", "piano.wav", "--pcm-full-scale-pa", "0"],
            vec!["--render", "piano.wav", "--pcm-full-scale-pa", "NaN"]] {
            assert!(options(&args).is_err(), "accepted {args:?}");
        }
    }
    #[test]
    fn bridge_trace_requires_one_key_and_physical_pressure() {
        let accepted = options(&["--preset", "steinway-d", "--render", "a.wav",
            "--note", "69", "--bridge-trace-csv", "paired.csv"]).unwrap();
        assert_eq!(accepted.bridge_trace_csv.as_deref(), Some("paired.csv"));
        let isolated = options(&["--preset", "steinway-d", "--render", "a.wav",
            "--performance", "held.performance", "--bridge-trace-csv", "paired.csv"]).unwrap();
        assert_eq!(isolated.performance.as_deref(), Some("held.performance"));
        for args in [
            vec!["--preset", "steinway-d", "--note", "69", "--bridge-trace-csv", "paired.csv"],
            vec!["--preset", "steinway-d", "--render", "a.wav", "--bridge-trace-csv", "paired.csv"],
            vec!["--render", "a.wav", "--note", "69", "--bridge-trace-csv", "paired.csv"],
            vec!["--preset", "steinway-d", "--render", "a.wav", "--note", "69",
                "--diagnostic-volume", "--bridge-trace-csv", "paired.csv"],
            vec!["--preset", "steinway-d", "--render", "a.wav", "--note", "69",
                "--bridge-trace-csv", "a.wav"],
        ] { assert!(options(&args).is_err(), "accepted {args:?}"); }
        let linear = [0.0, 0.01, 0.02, 0.03];
        for i in 0..linear.len() {
            assert!((bridge_acceleration(&linear, i, 100) - 1.0).abs() < 1e-12);
        }
    }
    #[test]
    fn modal_pressure_trace_requires_single_key_and_distinct_output() {
        let accepted = options(&["--preset", "steinway-d", "--render", "a.wav",
            "--note", "69", "--modal-pressure-csv", "modal.csv",
            "--bridge-trace-csv", "bridge.csv", "--microphone-right", "0.775,1,1"])
            .unwrap();
        assert_eq!(accepted.modal_pressure_csv.as_deref(), Some("modal.csv"));
        assert!(options(&["--preset", "steinway-d", "--render", "a.wav",
            "--performance", "held.performance", "--modal-pressure-csv", "modal.csv"]).is_ok());
        for args in [
            vec!["--preset", "steinway-d", "--note", "69", "--modal-pressure-csv", "modal.csv"],
            vec!["--preset", "steinway-d", "--render", "a.wav", "--modal-pressure-csv", "modal.csv"],
            vec!["--render", "a.wav", "--note", "69", "--modal-pressure-csv", "modal.csv"],
            vec!["--preset", "steinway-d", "--render", "a.wav", "--note", "69",
                "--diagnostic-volume", "--modal-pressure-csv", "modal.csv"],
            vec!["--preset", "steinway-d", "--render", "a.wav", "--note", "69",
                "--modal-pressure-csv", "a.wav"],
        ] { assert!(options(&args).is_err(), "accepted {args:?}"); }
    }
    #[test]
    fn damper_selection_composes_with_physical_and_midi_inputs_and_protects_its_source() {
        let o = options(&["--preset", "steinway-d", "--midi", "score.mid", "--render", "piano.wav",
            "--dampers", "pads.fspd", "--hammers", "felt.fsh"]).unwrap();
        assert_eq!(o.dampers.as_deref(), Some("pads.fspd"));
        assert!(options(&["--render", "piano.wav", "--dampers", "estimated"]).is_ok());
        assert_eq!(options(&[]).unwrap().dampers, None);
        for args in [vec!["--dampers"], vec!["--dampers", "estimated"],
            vec!["--render", "piano.wav", "--dampers", ""],
            vec!["--render", "pads.fspd", "--dampers", "pads.fspd"],
            vec!["--render", "piano.wav", "--dampers", "pads.fspd", "--dump-scale", "pads.fspd"],
            vec!["--render", "piano.wav", "--dampers", "estimated", "--dampers", "other.fspd"]] {
            assert!(options(&args).is_err(), "accepted {args:?}");
        }
    }
    #[test]
    fn rt0425_corrections_are_explicit_preset_inputs() {
        let o = options(&["--preset", "steinway-d", "--render", "piano.wav",
            "--rt0425-bridge-contacts", "--rt0425-hammer-stiffness"]).unwrap();
        assert!(o.rt0425_bridge_contacts && o.rt0425_hammer_stiffness);
        assert!(!options(&["--preset", "steinway-d", "--render", "piano.wav"])
            .unwrap().rt0425_bridge_contacts);
        for args in [
            vec!["--rt0425-bridge-contacts"],
            vec!["--rt0425-hammer-stiffness"],
            vec!["--preset", "steinway-d", "--rt0425-hammer-stiffness"],
            vec!["--preset", "steinway-d", "--board-geometry", "board.fsb", "--rt0425-bridge-contacts"],
            vec!["--preset", "steinway-d", "--hammers", "felt.fsh", "--render", "piano.wav", "--rt0425-hammer-stiffness"],
        ] { assert!(options(&args).is_err(), "accepted {args:?}"); }
        let source = options(&["--preset", "steinway-d", "--render", "piano.wav",
            "--rt0425-hammer-stiffness", "--rt0425-hammer-dissipation"]).unwrap();
        assert!(source.rt0425_hammer_dissipation);
        assert!(options(&["--preset", "steinway-d", "--render", "piano.wav",
            "--rt0425-hammer-dissipation"]).is_err());
    }
    #[test]
    fn stereo_receiver_is_explicit_physical_input_and_preserves_existing_controls() {
        let o=options(&["--preset","steinway-d","--render","stereo.wav",
            "--microphone","0.2,0.8,1","--microphone-right","1.2,0.8,1",
            "--midi","score.mid","--dampers","estimated"]).unwrap();
        assert_eq!(o.microphone,Some([0.2,0.8,1.0]));
        assert_eq!(o.microphone_right,Some([1.2,0.8,1.0]));
        assert_eq!(o.midi.as_deref(),Some("score.mid"));
        assert!(options(&["--preset","steinway-d","--render","s.wav",
            "--microphone-right","1,1,1"]).unwrap().microphone.is_none());
        for args in [vec!["--microphone-right","1,1,1"],
            vec!["--render","s.wav","--microphone-right","1,1,1"],
            vec!["--preset","steinway-d","--render","s.wav","--microphone-right","1,1,NaN"],
            vec!["--preset","steinway-d","--render","s.wav","--microphone-right","1,1,0"],
            vec!["--preset","steinway-d","--render","s.wav","--microphone-right","1,1"],
            vec!["--preset","steinway-d","--render","s.wav","--microphone-right","1,1,1","--diagnostic-volume"],
            vec!["--preset","steinway-d","--render","s.wav","--microphone-right","1,1,1","--microphone-right","2,1,1"]] {
            assert!(options(&args).is_err(),"accepted {args:?}");
        }
    }
    #[test]
    fn acoustic_refinement_is_an_explicit_bounded_pressure_control() {
        let ordinary=options(&["--preset","steinway-d","--render","p.wav"]).unwrap();
        assert_eq!(ordinary.acoustic_refinement_levels,0);
        let refined=options(&["--preset","steinway-d","--render","p.wav",
            "--acoustic-refinement-levels","2"]).unwrap();
        assert_eq!(refined.acoustic_refinement_levels,2);
        for args in [vec!["--acoustic-refinement-levels","1"],
            vec!["--render","p.wav","--acoustic-refinement-levels","1"],
            vec!["--preset","steinway-d","--render","p.wav","--acoustic-refinement-levels","4"],
            vec!["--preset","steinway-d","--render","p.wav","--acoustic-refinement-levels","1","--diagnostic-volume"],
            vec!["--preset","steinway-d","--render","p.wav","--acoustic-refinement-levels","1","--edge-cubic-board-mass"],
            vec!["--preset","steinway-d","--render","p.wav","--acoustic-refinement-levels","1","--acoustic-refinement-levels","2"]] {
            assert!(options(&args).is_err(),"accepted {args:?}");
        }
    }
    #[test]
    fn compact_receiver_csv_retains_frame_channel_order_and_physical_units() {
        for channels in [1,2] {
            let values=[0.125,-0.25,0.5,-1.0]; let mut output=Vec::new();
            write_receiver_pressure(&mut output,&values,48_000,channels).unwrap();
            let text=String::from_utf8(output).unwrap(); let mut lines=text.lines();
            assert_eq!(lines.next().unwrap(),if channels==1 {
                "sample,time_s,pressure_left_pa"
            } else {"sample,time_s,pressure_left_pa,pressure_right_pa"});
            for (sample,line) in lines.enumerate() {
                let fields=line.split(',').collect::<Vec<_>>();
                assert_eq!(fields[0].parse::<usize>().unwrap(),sample);
                assert_eq!(fields[1].parse::<f64>().unwrap(),sample as f64/48_000.0);
                for channel in 0..channels {
                    assert_eq!(fields[2+channel].parse::<f64>().unwrap(),values[sample*channels+channel]);
                }
            }
        }
        for (values,rate,channels) in [(&[0.0][..],48_000,2),(&[f64::NAN][..],48_000,1),
            (&[0.0][..],0,1),(&[0.0][..],48_000,0)] {
            let mut output=Vec::new();
            assert!(write_receiver_pressure(&mut output,values,rate,channels).is_err());
            assert!(output.is_empty());
        }
        assert!(options(&["--preset","steinway-d","--midi","score.mid","--render","p.wav",
            "--receiver-pressure-csv","receivers.csv"]).is_ok());
        for args in [vec!["--receiver-pressure-csv","receivers.csv"],
            vec!["--render","p.wav","--receiver-pressure-csv","receivers.csv"],
            vec!["--preset","steinway-d","--render","p.wav","--diagnostic-volume","--receiver-pressure-csv","receivers.csv"],
            vec!["--preset","steinway-d","--render","p.wav","--receiver-pressure-csv","p.wav"],
            vec!["--preset","steinway-d","--midi","score.mid","--render","p.wav","--receiver-pressure-csv","score.mid"]] {
            assert!(options(&args).is_err(),"accepted {args:?}");
        }
    }
    #[test]
    fn midi_options_compose_with_physical_inputs_without_overwriting_the_score() {
        let o = options(&["--preset", "steinway-d", "--midi", "score.mid", "--render", "piano.wav",
            "--midi-channel", "16", "--midi-velocity-max-m-s", "2", "--midi-half-pedal",
            "--hammers", "felt.fsh"]).unwrap();
        assert_eq!(o.midi.as_deref(), Some("score.mid"));
        assert_eq!(o.midi_mapping.channel, 15);
        assert_eq!(o.midi_mapping.maximum_velocity_m_s, 2.0);
        assert!(o.midi_mapping.continuous_sustain);
        for args in [vec!["--midi", "score.mid"], vec!["--midi-channel", "1"],
            vec!["--midi-half-pedal"], vec!["--midi-velocity-max-m-s", "2"],
            vec!["--midi", "score.mid", "--render", "score.mid"],
            vec!["--midi", "score.mid", "--render", "p.wav", "--performance", "p.csv"],
            vec!["--midi", "score.mid", "--render", "p.wav", "--note", "69"],
            vec!["--midi", "score.mid", "--render", "p.wav", "--velocity", "2"],
            vec!["--midi", "score.mid", "--render", "p.wav", "--midi-channel", "0"],
            vec!["--midi", "score.mid", "--render", "p.wav", "--midi-channel", "17"],
            vec!["--midi", "score.mid", "--render", "p.wav", "--midi-velocity-max-m-s", "NaN"],
            vec!["--midi", "score.mid", "--render", "p.wav", "--midi-velocity-max-m-s", "9"]] {
            assert!(options(&args).is_err(), "accepted {args:?}");
        }
    }
    #[test]
    fn midi_and_si_csv_drive_identical_preset_contact_pedals_and_energy() {
        // PPQN 480, default tempo: four ticks -> sample 200 at 48 kHz.
        // End-of-track at tick 12 releases sustain at sample 600, not PCM.
        let bytes = b"MThd\0\0\0\x06\0\0\0\x01\x01\xe0MTrk\0\0\0\x10\
            \0\xb0\x40\x7f\0\x90\x45\x7f\x04\x90\x45\0\x08\xff\x2f\0";
        let mut o = options(&["--preset", "steinway-d", "--render", "piano.wav"]).unwrap(); o.modes = 12;
        for spatial in [false,true] {
        o.dampers = spatial.then(|| "estimated".to_owned());
        let course = selected_scale(None, &o).unwrap()[48];
        let mut a = prepare_instrument(vec![course], &board::demonstration(), &o).unwrap();
        let mut b = prepare_instrument(vec![course], &board::demonstration(), &o).unwrap();
        assert_eq!(a.damper_resolution().is_some(), spatial);
        let decoded = midi::read(bytes, &[69], 48_000, 1500,
            midi::Mapping { maximum_velocity_m_s: 2.0, ..midi::Mapping::default() }).unwrap();
        assert_eq!(decoded.report.end_releases, 1);
        let mut imported = performance::Performance::from_events(decoded.events, &[69], 1500).unwrap();
        let mut csv = performance::Performance::read("sample,event,key,value\n0,sustain,0,1\n\
            0,note_on,69,2\n200,note_off,69,0\n600,sustain,0,0\n", &[69], 1500).unwrap();
        for sample in 0..1500 {
            imported.dispatch(sample, &mut a).unwrap(); csv.dispatch(sample, &mut b).unwrap();
            assert_eq!(a.step().unwrap().to_bits(), b.step().unwrap().to_bits());
        }
        assert!(a.accounting.felt_loss_j > 0.0);
        assert_eq!(a.accounting.input_work_j.to_bits(), b.accounting.input_work_j.to_bits());
        assert!((a.accounting.input_work_j - a.energy_j() - a.accounting.dissipated_j()).abs() < 1e-7);
        assert!(a.accounting.damper_loss_j > 0.0);
        }
    }
    #[test]
    fn rendering_and_both_physical_imports_are_composable() {
        let o = options(&["--scale", "strings.csv", "--render", "piano.wav", "--board", "board.csv",
            "--sample-rate", "44100", "--substeps", "2", "--modes", "40", "--note", "21"]).unwrap();
        assert_eq!(o.render.as_deref(), Some("piano.wav"));
        assert_eq!(o.scale.as_deref(), Some("strings.csv"));
        assert_eq!(o.board.as_deref(), Some("board.csv"));
        assert_eq!((o.sample_rate, o.substeps, o.modes, o.note), (44_100, 2, 40, Some(21)));
    }
    #[test]
    fn imported_subset_reaches_the_prepared_instrument_without_substitution() {
        let mut course = geometry::demonstration_scale().unwrap()[0];
        course.length_m *= 0.97;
        let scale = load_scale(Some(&geometry::write_scale(&[course]))).unwrap();
        assert_eq!(scale, vec![course]);
        assert_eq!(study_key(&scale, None).unwrap(), 21);
        assert!(study_key(&scale, Some(69)).is_err());
        let mut input = board::demonstration();
        input.truncate(1);
        input[0].frequency_hz = 123.0;
        let modes = load_board(Some(&board::write(&input)), &scale).unwrap();
        assert_eq!(modes[0].frequency_hz, 123.0);
        let mut piano = engine::Instrument::new(scale, &modes, 8_000, 4, 4, true).unwrap();
        assert_eq!(piano.bank.board_count, 1);
        piano.note_on(21, 1.0).unwrap();
        for _ in 0..16 { assert!(piano.step().unwrap().is_finite()); }
    }
    #[test]
    fn invalid_inputs_never_fall_back_to_estimates() {
        assert!(load_scale(Some("bad CSV")).is_err());
        let scale = geometry::demonstration_scale().unwrap();
        assert!(load_board(Some(board::HEADER), &scale).is_err());
        for args in [vec!["--velocity", "NaN"], vec!["--duration", "inf"], vec!["--modes", "0"],
            vec!["--scale", "input.csv", "--render", "input.csv"], vec!["--render"],
            vec!["--note", "69", "--note", "70"], vec!["--substeps", "17"]] {
            assert!(options(&args).is_err(), "accepted {args:?}");
        }
    }
    #[test]
    fn geometric_board_and_performance_options_compose_without_replacing_existing_controls() {
        let o = options(&["--scale", "strings.csv", "--board-geometry", "panel.fsb",
            "--board-band-hz", "300", "--performance", "score.csv", "--render", "piano.wav",
            "--sample-rate", "44100", "--duration", "2", "--modes", "40"]).unwrap();
        assert_eq!(o.board_geometry.as_deref(), Some("panel.fsb"));
        assert_eq!(o.performance.as_deref(), Some("score.csv"));
        assert_eq!(o.board_band_hz, 300.0);
        assert_eq!((o.sample_rate, o.duration, o.modes), (44_100, 2.0, 40));
        for args in [vec!["--board", "a.csv", "--board-geometry", "b.fsb"],
            vec!["--board-band-hz", "300"], vec!["--performance", "score.csv"],
            vec!["--render", "a.wav", "--performance", "score.csv", "--note", "69"],
            vec!["--render", "a.wav", "--performance", "score.csv", "--velocity", "1"],
            vec!["--render", "score.csv", "--performance", "score.csv"],
            vec!["--observer-gain", "NaN"],
            vec!["--board-geometry", "a.fsb", "--board-band-hz", "21600"]] {
            assert!(options(&args).is_err(), "accepted {args:?}");
        }
    }
    #[test]
    fn board_export_preserves_missing_measurements_and_reimports_the_admitted_subset() {
        let scale = vec![geometry::demonstration_scale().unwrap()[48]];
        let modes = board::demonstration();
        let text = write_board_for_scale(&modes, &scale);
        assert_eq!(text.lines().count(), modes.len() + 1);
        let roundtrip = board::read(&text, &[69]).unwrap();
        for (a, b) in roundtrip.iter().zip(&modes) {
            assert_eq!(a.frequency_hz, b.frequency_hz);
            assert_eq!(a.bridge[48], b.bridge[48]);
            assert_eq!(a.volume, b.volume);
        }
        assert!(board::read(&text, &[60, 69]).is_err());
    }
    #[test]
    fn source_preset_renders_and_exports_without_conflicting_board_inputs() {
        let o = options(&["--preset", "steinway-d", "--render", "d.wav", "--board-band-hz", "350",
            "--dump-geometry", "d.fsb", "--dump-obj", "d.obj", "--mesh-divisions", "12",
            "--equilibrate-board-mass", "--consistent-board-mass"]).unwrap();
        assert_eq!(o.preset.as_deref(), Some("steinway-d"));
        assert_eq!(o.mesh_divisions, 12);
        assert!(o.equilibrate_board_mass);
        assert!(o.consistent_board_mass);
        assert!(options(&["--preset", "steinway-d", "--consistent-board-mass"]).is_err());
        assert!(options(&["--render", "d.wav", "--consistent-board-mass"]).is_err());
        let cubic = options(&["--preset", "steinway-d", "--render", "d.wav",
            "--edge-cubic-board-mass"]).unwrap();
        assert!(cubic.edge_cubic_board_mass);
        for args in [
            vec!["--preset", "steinway-d", "--edge-cubic-board-mass"],
            vec!["--render", "d.wav", "--edge-cubic-board-mass"],
            vec!["--preset", "steinway-d", "--render", "d.wav",
                "--edge-cubic-board-mass", "--consistent-board-mass"],
        ] { assert!(options(&args).is_err()); }
        for args in [vec!["--preset", "unknown"],vec!["--preset", "steinway-d", "--board", "b.csv"],
            vec!["--preset", "steinway-d", "--board-geometry", "b.fsb", "--dump-obj", "d.obj"],vec!["--dump-obj", "d.obj"],
            vec!["--preset", "steinway-d", "--dump-obj", "d.obj", "--render", "d.obj"],
            vec!["--preset", "steinway-d", "--mesh-divisions", "3"],
            vec!["--preset", "steinway-d", "--mesh-divisions", "33"],
            vec!["--equilibrate-board-mass"]] {assert!(options(&args).is_err());}
    }
    #[test]
    fn physical_microphone_and_diagnostic_observer_are_not_silently_mixed() {
        let o=options(&["--preset","steinway-d","--render","d.wav","--microphone","0.5,1.2,0.8"]).unwrap();
        assert_eq!(o.microphone,Some([0.5,1.2,0.8]));assert!(!o.diagnostic_volume);
        let o=options(&["--preset","steinway-d","--render","d.wav","--diagnostic-volume","--observer-gain","500"]).unwrap();
        assert!(o.diagnostic_volume);
        for args in [vec!["--preset","steinway-d","--observer-gain","1000"],
            vec!["--render","d.wav","--microphone","0,0,1"],
            vec!["--preset","steinway-d","--render","d.wav","--microphone","0,0,NaN"],
            vec!["--preset","steinway-d","--render","d.wav","--microphone","0,0,-1"],
            vec!["--preset","steinway-d","--render","d.wav","--microphone","0,1"],
            vec!["--preset","steinway-d","--render","d.wav","--microphone","0,0,1","--diagnostic-volume"]] {
            assert!(options(&args).is_err());
        }
    }
    #[test]
    fn mass_equilibration_admits_the_refined_source_board() {
        let source = steinway_d::build(20).unwrap();
        let keys: Vec<u8> = (21..=108).collect();
        let prepared = prepare_geometric_board(&source.geometry, &keys, 1200.0, true, false, false, 0)
            .expect("mass-equilibrated mesh-20 source board");
        assert_eq!(prepared.free_dofs, 4236);
        assert_eq!(prepared.modes.len(), 43);
        assert!(prepared.frequency_intervals_hz.iter().all(|(lo, hi)|
            lo.is_finite() && hi.is_finite() && *lo > 0.0 && *lo <= *hi && *hi <= 1200.0));
    }
    #[test]
    fn consistent_panel_mass_reaches_the_real_board_solve_without_changing_geometry() {
        let source = steinway_d::build(8).unwrap();
        let keys: Vec<u8> = (21..=108).collect();
        let lumped = prepare_geometric_board(&source.geometry, &keys, 400.0, true, false, false, 0).unwrap();
        let consistent = prepare_geometric_board(&source.geometry, &keys, 400.0, true, true, false, 0).unwrap();
        let cubic = prepare_geometric_board(&source.geometry, &keys, 400.0, true, false, true, 0).unwrap();
        assert_eq!(lumped.area_m2, consistent.area_m2);
        assert_eq!(lumped.mass_kg, consistent.mass_kg);
        assert_eq!(lumped.free_dofs, consistent.free_dofs);
        assert_eq!(lumped.area_m2, cubic.area_m2);
        assert_eq!(lumped.mass_kg, cubic.mass_kg);
        assert_eq!(lumped.free_dofs, cubic.free_dofs);
        assert!(cubic.provenance.contains("cubic edge-compatible panel mass"));
        assert!(!cubic.modes.is_empty());
        assert_ne!(lumped.modes[0].frequency_hz, cubic.modes[0].frequency_hz);
        assert_eq!(cubic.surface.len(), lumped.surface.len());
        let sampled_area = cubic.surface.iter().map(|sample| sample.area_m2).sum::<f64>();
        assert!((sampled_area - cubic.area_m2).abs() < 1e-11);
        for (index, mode) in cubic.modes.iter().enumerate() {
            let sampled_volume = cubic.surface.iter()
                .map(|sample| sample.area_m2 * sample.mode_shape[index]).sum::<f64>();
            assert!((sampled_volume - mode.volume).abs() < 1e-11,
                "mode {index} acoustic surface and exact volume disagree");
        }
        assert!(consistent.provenance.contains("exact P1 transverse panel mass"));
        assert!(lumped.provenance.contains("mass-diagonal solver equilibration"));
        assert!(!lumped.provenance.contains("exact P1 transverse panel mass"));
        assert!(!lumped.modes.is_empty() && !consistent.modes.is_empty());
        assert_ne!(lumped.modes[0].frequency_hz, consistent.modes[0].frequency_hz);
    }
    #[test]
    fn preset_tunes_tension_not_geometry_and_raw_source_is_available() {
        let o=options(&["--preset","steinway-d"]).unwrap();
        let tuned=selected_scale(None,&o).unwrap();let source=steinway_scale::courses().unwrap();
        for (a,b) in tuned.iter().zip(&source) {
            assert_eq!(a.length_m,b.length_m);assert_eq!(a.linear_density_kg_m,b.linear_density_kg_m);
            assert_eq!(a.flexural_rigidity_nm2,b.flexural_rigidity_nm2);
            let target=440.0*2.0f64.powf((f64::from(a.midi)-69.0)/12.0);
            assert!((a.partial_hz(1,a.tension_n)/target-1.0).abs()<1e-10);
        }
        let raw=options(&["--preset","steinway-d","--raw-tensions"]).unwrap();
        assert_eq!(selected_scale(None,&raw).unwrap(),source);
        let imported=options(&["--preset","steinway-d","--scale","measured.csv"]).unwrap();
        assert_eq!(selected_scale(Some(&geometry::write_scale(&source)),&imported).unwrap(),source);
        for args in [vec!["--concert-pitch","NaN"],vec!["--concert-pitch","400"],
            vec!["--raw-tensions"],vec!["--preset","steinway-d","--raw-tensions","--concert-pitch","442"]] {
            assert!(options(&args).is_err());
        }
    }
    #[test]
    fn source_hammer_cards_are_used_by_the_render_preparation_path() {
        let mut o=options(&["--preset","steinway-d"]).unwrap();o.modes=12;
        let scale=selected_scale(None,&o).unwrap();
        let mut piano=prepare_instrument(vec![scale[48]],&board::demonstration(),&o).unwrap();
        piano.note_on(69,2.0).unwrap();
        for _ in 0..1500 {assert!(piano.step().unwrap().is_finite());}
        assert!(piano.accounting.felt_loss_j>0.0);
        assert!(piano.accounting.felt_relaxation_loss_j>0.0);
        assert!((piano.accounting.input_work_j-piano.energy_j()-piano.accounting.dissipated_j()).abs()<1e-7);
    }
    #[test]
    fn jack_performance_reaches_preset_mechanics_without_velocity_substitution() {
        let mut o=options(&["--preset","steinway-d"]).unwrap();o.modes=12;
        let c=selected_scale(None,&o).unwrap()[48];
        let mut piano=prepare_instrument(vec![c],&board::demonstration(),&o).unwrap();
        let mut score=performance::Performance::read("sample,event,key,value\n0,jack_legato,69,30\n",&[69],2400).unwrap();
        for sample in 0..2400 {score.dispatch(sample,&mut piano).unwrap();piano.step().unwrap();}
        assert!(piano.accounting.shank_loss_j>0.0);assert!(piano.accounting.felt_loss_j>0.0);
        assert!((piano.accounting.input_work_j-piano.energy_j()-piano.accounting.dissipated_j()).abs()<1e-7);
    }
    #[test]
    fn supplied_hammers_compose_with_presets_and_protect_the_input_file() {
        let o = options(&["--preset", "steinway-d", "--hammers", "felt.fsh", "--render", "d.wav"]).unwrap();
        assert_eq!(o.hammers.as_deref(), Some("felt.fsh"));
        assert!(options(&["--hammers", "felt.fsh"]).is_err());
        assert!(options(&["--hammers", "felt.fsh", "--render", "felt.fsh"]).is_err());
        let course = selected_scale(None, &o).unwrap()[48];
        let text = "frankensim-hammer-materials-v1\nfelt,69,400000,0.2,2.5,3.2,0.25,0.8,2500000\n";
        let mut piano = prepare_instrument_with_hammers(vec![course], &board::demonstration(),
            &o, Some(text)).unwrap();
        // An imported elastic card must not resurrect the preset Prony tails,
        // while the preset's geometry-derived jack/shank remains available.
        piano.jack_on(69, 30.0, 0.1).unwrap();
        for _ in 0..1500 { assert!(piano.step().unwrap().is_finite()); }
        assert_eq!(piano.accounting.felt_relaxation_loss_j, 0.0);
        assert!(piano.accounting.shank_loss_j > 0.0);
        assert!(prepare_instrument_with_hammers(vec![course], &board::demonstration(),
            &o, Some("invalid material")).is_err());
    }
}

#[cfg(test)]
#[path = "string_render_tests.rs"]
mod string_render_tests;
