//! Geometry-first grand-piano composition of existing FrankenSim owners.
//! `cargo run --release -p fs-couple --example grand_piano -- --preset steinway-d --render piano.wav`
//! Source-derived geometry is approximate, not a measured specimen digital twin.
mod geometry;
pub(crate) mod linear;
mod board;
mod board_geometry;
mod crowned_board;
mod steinway_d;
mod steinway_scale;
// --felt-pack-hammers. The default felt cards remain the steinway_scale ones.
#[allow(dead_code)]
mod hammer_felt;
mod stretch_tuning;
mod performance;
use performance::midi;
mod felt;
pub(crate) mod engine;
mod microphone;
mod pressure_basis;
pub(crate) mod audio;
mod wav_output;
mod hammer_materials;
mod string_sections;
mod string_polarization;
mod cavity;
#[path = "../../src/bin/music_render/json.rs"]
mod json;
use json::json_string;
// The shared native gas-source loader: --gas-component/--gas-water/--gas-phase.
#[allow(dead_code)]
#[path = "../../src/bin/music_render/thermochemical.rs"]
mod thermochemical;
#[cfg(test)]
#[allow(dead_code)]
#[path = "../../src/render/schedule/reed/field_input/tests/support.rs"]
mod thermochemical_support;
mod rayleigh_loading;
mod ambient;

const USAGE: &str = "grand_piano [--render piano.wav] [--scale strings.csv]
    [--preset steinway-d] [--board board.csv | --board-geometry panel.fsb|panel.fss]
    [--rt0425-bridge-contacts | --estimated-bridge-contacts] [--rt0425-hammer-stiffness]
    [--rt0425-hammer-dissipation] [--felt-pack-hammers] [--felt-pack hammers.fsmatpk] [--legacy-felt-cards]
    [--rt0425-string-damping]
    [--board-panel-pack panel.fsmatpk] [--board-rib-pack ribs.fsmatpk] [--board-bridge-pack bridges.fsmatpk]
    [--equilibrate-board-mass] [--consistent-board-mass | --edge-cubic-board-mass]
    [--hammers materials.fsh] [--string-sections bonded.fsss]
    [--hammer-footprints faces.fshp]
    [--dampers estimated | point | pads.fspd] [--fast-dampers] [--string-stretching axial.fspx]
    [--string-polarization bridge-frames.fspp]
    [--cavity enclosure.fspc] [--gas-record air.gas]
    [--gas-component SPECIES.pack NASA9.pack TRANSPORT.pack ... --gas-water none|ID --gas-phase frozen|liquid-water]
    [--rayleigh-loading sides,min_hz,max_hz,samples]
    [--radiation-ports 1..128]
    [--concert-pitch 430..450 | --stretch-tuning 430..450 | --raw-tensions]
    [--mesh-divisions 4..32] [--dump-geometry panel.fsb] [--dump-obj soundboard.obj]
    [--board-band-hz Hz] [--board-reduction SPEC]
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
    [--sample-rate Hz] [--substeps 1..16] [--modes auto|1..512]
    [--pcm-full-scale-pa positive-Pa] [--wav-format pcm16|pcm24|float32]
    [--dump-scale strings.csv] [--dump-board board.csv]
--preset steinway-d reconstructs the published 17-rib Model D drawing, with
spruce panel, sugar-pine ribs, maple bridges, cut-off bar and 88 bridge stations.
--board-geometry may replace that board while retaining the preset strings,
felt cards and shank mechanics. Its native header selects flat FSB or crowned
3-D shell FSS. An invalid/missing supplied board never falls back to the preset.
Mesh-generation/export controls cannot accompany a supplied board override.
The preset cannot be combined with --board modal CSV.
The preset's strings meet its bridges at the source's 84 published coupling
points, projected onto this approximate board's bridges; four end keys remain
extrapolated (--rt0425-bridge-contacts names this default). The former estimate,
--estimated-bridge-contacts, spaces stations evenly along each bridge: with the
published string lengths it puts most agraffes over the soundboard itself and
reverses the bass order, so it only reproduces earlier renders.
--rt0425-hammer-stiffness applies the source K_H to each unison string.
--rt0425-hammer-dissipation requires that stiffness and uses the published
per-key R_H d(e^p)/dt instead of estimated crush and Prony relaxation.
By default the preset plays the felt pack's hammers, as the browser does:
Stulov's 88 measured Abel hammers, each the whole hammer's K u^p + R d(u^p)/dt
shared by its unison, with the pack's hammer mass (--felt-pack-hammers names
this default). The RT-0425 hammer flags, --hammers or --legacy-felt-cards (the
former default card, RT-0425's K_H divided across the unison with estimated
crush and Prony terms) choose another card instead.
--felt-pack FILE supplies that felt pack (default: the built-in Stulov Abel set).
--board-panel-pack, --board-rib-pack and --board-bridge-pack FILE supply the preset
board's fs-matdb packs (default: the built-in FPL Sitka spruce, sugar pine and
hard maple); the cut-off bar takes the rib pack. The board's wood moisture
follows the selected air (--gas-record or the committed humid record).
A `thermochemical ...` --gas-record (Ne, SF6, He or any packed species)
resolves only through the --gas-component source packs it names, with explicit
--gas-water and --gas-phase; without them it refuses. Nothing is defaulted.
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
--acoustic-refinement-levels uniformly subdivides flat radiating triangles for
Rayleigh integration, evaluating the same P1 or edge-cubic displacement field
at every new point. It preserves the structural mesh, modes and bridge mechanics.
Level 0 is unchanged; the 120000-point receiver budget still applies.
Crowned fields and diagnostic-volume output are excluded.
It uses Chabassier/Durufle's wrapped-string MODEL table (84 notes plus four
estimated extensions), separate per-key hammer force and relaxation cards, and
published shank geometry reduced to rigid rotation plus one bending coordinate.
The shank is not a fitted oscillator: mass, compliance and jack projection come
from its dimensions/material. This is a linearized reduction, not the full action.
Lengths, effective winding mass and EI remain fixed. Preset tensions are tuned
to A4=440 Hz by default to compensate rounded source values; --raw-tensions
preserves the published table. --concert-pitch tunes first partials by changing
physical tension, not oscillator frequencies. --stretch-tuning instead anchors
C4..B4 at the selected A4 and matches actual stiff-string 2:1 octave partials
outward across the supplied scale. It excludes --concert-pitch/--raw-tensions.
This is a declared octave policy, not a measured or perceptual tuning fit.
An explicit --scale always supplies the geometry/masses; absent either tuning option,
its tensions are preserved even with a preset. Preset hammer voicing still applies
unless --hammers supplies a complete per-key WoolFelt/Prony material file.
--hammers requires --render and a card for every key in the supplied scale, not
only the struck keys. Missing, duplicate or invalid cards refuse without fallback.
These cards change physical contact forces and relaxation, not an output EQ.
The preset shank and the scale's hammer mass/patch geometry remain unchanged.
See HAMMERS.md for the SI format; importing values does not certify measurements.
--string-sections overrides the inertial, bending and axial coefficients of every
selected speaking course using immutable fs-material source-claim layers.
Both input versions require an ordered, complete solid concentric section per key.
V2 also supplies each layer's solid bending viscosity and explicit rad/s support.
The section integral replaces the demonstration bending-loss spectrum on every
retained speaking/duplex mode and both transverse directions. The separately
authored common velocity drag remains; v1 retains the demonstration bending law.
Neither format is a surrogate for a helically wound bass string. Tension remains unchanged
unless explicit tuning is requested. The resulting scale is passed unchanged
to modal mechanics, hammer/bridge coupling and observation, with no timbre knob.
With --string-stretching, use section in place of numeric EA so extension and
longitudinal waves use the same layer integrals. Slope/strain limits remain
explicit. This mode excludes RT-0425 measured string damping, which is not
identified for a substituted material.
--hammer-footprints supplies point, uniform span or authored crown profile for
EVERY key. Profiles supply one to four ordered sites with longitudinal offsets,
face recession, local felt thickness and positive fractions of the original area.
Sites engage according to their gaps and retain independent felt/Prony histories
while sharing the hammer. Point/span thickness still comes from the scale.
An optional normal,key,primary,secondary row supplies a unit felt-face normal.
An oblique normal requires --string-polarization; the existing hammer trajectory
is preserved and contact work projects into both transverse string directions.
It requires --render and works with supplied materials, MIDI and existing pedals.
Published R_H requires the original uniform thickness. Varying thickness applies
the selected felt/Prony law locally and cannot silently reuse that source rate.
See HAMMER_FOOTPRINTS.md for SI rows and the fixed-normal contact model.
--dampers selects spatial viscous pads or finite-gap wool-felt pads; the file
header chooses the model. The steinway-d preset defaults to the existing
estimated finite-width viscous pads, removing the point image's artificial
undamped partials at its nodal station. 'estimated' selects these approximate
spans and drag explicitly; 'point' selects the original point approximation.
Without a preset, omitted --dampers retains point damping. A file must cover
every scale key with a pad or explicit free row and overrides either default.
Felt pads additionally require supplied geometry, travel speed and material
cards. MIDI/CSV key, sustain and sostenuto controls move their contact height;
felt contact, stored energy and signed actuator work enter the mechanics solve.
It requires --render. See DAMPERS.md for SI formats and model limits.
--fast-dampers carries the spatial viscous pads' flows on reduced coordinates:
each pad's stations move the board along one bridge direction, so a station
needs only its string's modes and one scalar. Same flows in exact arithmetic,
different rounding; about 3x less damper work. Off by default (the reference).
--string-polarization adds both transverse directions of each physical string.
Its complete per-key file supplies bridge sites, 3-D arms, string/hammer axes
and lateral damper ratios. The same board solve projects both bridge rows;
the primary row must agree with the board's existing bridge geometry.
It requires a geometric --render, using the selected P1, edge-cubic or crowned
motion field. Modal CSV lacks that motion. No lateral coupling or drag is guessed.
See STRING_POLARIZATION.md for the physical input format and scope.
--cavity supplies a sealed rectangular air cavity below a flat geometric
soundboard. Its source, dimensions, air state, retained pressure modes and
momentum damping are explicit. Uniform compression and standing waves react
on the played board in the same mechanical clock, alongside strings/contact.
The board may cover part of the top face; the remainder is a rigid wall.
It requires --render and full geometric motion. See CAVITY.md for the format;
no cavity geometry or sealed-enclosure assumption is inferred from the preset.
--string-stretching supplies linear, averaged-extension (stretch), or finite
longitudinal-wave (waves) selection for EVERY scale key. All nonlinear rows
supply axial rigidity EA in N (or section for --string-sections) and a
moderate-slope bound. A waves row also supplies axial band, mode cap, damping
ratio and material-strain bound; these
are not inferred from tension, EI or winding mass. Axial modes share the
existing string/contact clock and react through the moving transverse bridge;
the speaking end's axial force also drives the bridge through the declared
downbearing/bridge-cap projection. It requires --render; all-linear selection
preserves the original image. See STRING_STRETCHING.md for SI rows, explicit
retention limits and the distinction from a full longitudinal bridge model.
--dump-geometry/--dump-obj export the board model; export alone skips eigenanalysis.
--dump-board of a flat geometric board without --render also writes slope_x_per_m_sqrt_kg and
slope_y_per_m_sqrt_kg: the gradient (dw/dx, dw/dy) of each mass-normalized mode shape at the key's
bridge station, from the same retained nodal rotations and map as the bridge coefficient. A table
with them reads back as before: --board checks the two columns and drops them.
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
--board-reduction SPEC explicitly permits up to 512 source modes on a flat or crowned board.
SPEC is max_modes,keep_low_modes,Hz,... or
adaptive,max_modes,keep_low_modes,tolerance,Hz,... (0 < tolerance < 1).
It retains at most max_modes (1..128), including keep_low_modes exact low modes. Static and damped responses at 1..16 increasing target frequencies guide
the remaining basis at every admitted primary bridge and both directions of
each supplied --string-polarization frame. Targets must lie within
--board-band-hz; that source band is never widened implicitly. Adaptive selection
checks actual complex compliance between every pair of admitted bridge ports
at source resonances, midpoints, static and seed frequencies. Each error is
scaled by the geometric mean of the source driving-point magnitudes, including
the original relative error on the diagonal. Enrichment continues until the
tolerance is met; insufficient rank or work budget refuses preparation. This
finite source-model grid does not certify continuous-band or acoustic accuracy.
SPEC pod,max_modes,keep_low_modes,count,top_hz instead keeps the leading proper
orthogonal decomposition of every port's static and count (2..16) damped responses
log-spaced from the first unprotected source mode to top_hz, which may exceed
--board-band-hz: string partials drive the slice there too. See BOARD_REDUCTION.md.
Full projected wood damping, bridge motion and the
radiating field stay in the same basis.
For supplied downbearing, reduction uses the solved equilibrium tangent modes.
This requires --render and excludes --dump-board: modal CSV cannot preserve the
dense material operator. Flat-board mass controls remain flat-only. See BOARD_REDUCTION.md.
--modes defaults to auto: each string retains its physical partials through
45% of the OUTPUT sample rate, up to the 512-partial work ceiling. This avoids
the former 24-partial bass cutoff. A numeric value (for example --modes 24)
explicitly selects a smaller budget. Any segments still limited by the count
are reported; oversampling does not widen the output band. These budgets and the
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
This assumes an infinite baffle, with no lid/room scattering.
--rayleigh-loading adds physical air backreaction to a flat geometric pressure
render. Supply sides,min_hz,max_hz,samples: one or two independent half-spaces,
an explicit frequency band and an odd sample count from 33 to 257. The band
must cover the loaded board reference frequencies and stay below 45% of the
output rate. The full reciprocal Rayleigh impedance (up to 128 board modes) is
sampled linearly in frequency and fitted by the existing passive owner with
resonators that resolve the board's own travel time (spacing c/6D, D its largest
extent); training samples must be that close, so 90-2400 Hz on the concert
grand's 2.36 m board needs at least 193 samples. It is coupled on the mechanics clock. Supplied wood damping is
treated as intrinsic; the front microphone remains a one-side observation.
Two-sided loading excludes a sealed --cavity below the board. Crowned shells,
diagnostic-volume output and unsupported work/rate/fit budgets refuse.
Absent this option, pressure observation has no air backreaction. See
RAYLEIGH_LOADING.md for quadrature, fitting and finite-band limitations.
--radiation-ports explicitly bounds geometry-derived acoustic combinations for
--rayleigh-loading. All retained board coordinates and receiver inputs remain;
complete lifted complex/resistance errors must pass the same limits. An
insufficient rank refuses; it never increases rank or reduces the soundboard.
Pressure uses every mechanics substep with bandlimited acceleration and
anti-alias filtering before output-rate propagation. At 4x oversampling the
numerical delay is 84 audio samples (52 differentiation/decimation plus 32
fractional-flight interpolation), in addition to acoustic travel time.
Histories persist across blocks.
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
--pcm-full-scale-pa declares the pressure mapped to digital full scale (default 2 Pa).
--wav-format writes PCM16 (default), PCM24 or Float32 directly from computed
pressure, preserving quiet decay without an intermediate PCM16/CSV conversion.
The same explicit pressure scale applies to both channels. Float32 can retain
finite values beyond +/-1; attenuate them before fixed-point playback.
All output paths must be fresh, with existing parent directories. Equivalent
paths through dot components or symlinked parents refuse before preparation.
An over-range fixed-point render refuses before writing a clipped WAV; this option changes
encoding gain, not the mechanics, microphone position, or acoustic calibration.";

#[derive(Debug)]
struct Options {
    render: Option<String>, scale: Option<String>, board: Option<String>,
    board_geometry: Option<String>, performance: Option<String>, preset: Option<String>,
    hammers: Option<String>, string_sections: Option<String>, hammer_footprints: Option<String>, dampers: Option<String>,
    string_stretching: Option<String>, string_polarization: Option<String>, cavity: Option<String>,
    gas_record: Option<String>,
    gas_sources: thermochemical::Options,
    midi: Option<String>, midi_mapping: midi::Mapping,
    concert_pitch: Option<f64>, stretch_tuning: Option<f64>, raw_tensions: bool, fast_dampers: bool,
    rt0425_bridge_contacts: bool, estimated_bridge_contacts: bool, rt0425_hammer_stiffness: bool,
    rt0425_hammer_dissipation: bool, felt_pack_hammers: bool, felt_pack: Option<String>, legacy_felt_cards: bool,
    board_panel_pack: Option<String>, board_rib_pack: Option<String>, board_bridge_pack: Option<String>,
    rt0425_string_damping: bool,
    equilibrate_board_mass: bool, consistent_board_mass: bool, edge_cubic_board_mass: bool,
    mesh_divisions: usize, dump_geometry: Option<String>, dump_obj: Option<String>,
    board_band_hz: f64, observer_gain: f64,
    board_reduction: Option<board_geometry::ritz::RitzOptions>,
    rayleigh_loading: Option<rayleigh_loading::Options>,
    radiation_ports: Option<usize>,
    microphone: Option<[f64; 3]>, microphone_right: Option<[f64; 3]>, diagnostic_volume: bool,
    acoustic_refinement_levels: usize,
    dump_scale: Option<String>, dump_board: Option<String>,
    bridge_trace_csv: Option<String>, modal_pressure_csv: Option<String>, receiver_pressure_csv: Option<String>,
    pressure_basis_json: Option<String>,
    note: Option<u8>, velocity: Option<f64>, duration: f64,
    sample_rate: u32, substeps: usize, modes: usize, pcm_full_scale_pa: f64, help: bool,
    wav_format: wav_output::Format,
}
impl Default for Options {
    fn default() -> Self {
        Self { render: None, scale: None, board: None, board_geometry: None,
            performance: None, preset: None, hammers: None, string_sections: None, hammer_footprints: None, dampers: None, concert_pitch: None, stretch_tuning: None, raw_tensions: false, fast_dampers: false,
            rt0425_bridge_contacts: false, estimated_bridge_contacts: false, rt0425_hammer_stiffness: false,
            rt0425_hammer_dissipation: false, felt_pack_hammers: false, felt_pack: None, legacy_felt_cards: false,
            board_panel_pack: None, board_rib_pack: None, board_bridge_pack: None,
            rt0425_string_damping: false,
            equilibrate_board_mass: false, consistent_board_mass: false, edge_cubic_board_mass: false,
            string_stretching: None, string_polarization: None, cavity: None, gas_record: None,
            gas_sources: thermochemical::Options::default(),
            midi: None, midi_mapping: midi::Mapping::default(),
            mesh_divisions: 8, dump_geometry: None, dump_obj: None,
            board_band_hz: 400.0, observer_gain: 10_000.0, dump_scale: None,
            board_reduction: None, rayleigh_loading: None, radiation_ports: None,
            microphone: None, microphone_right: None, diagnostic_volume: false,
            acoustic_refinement_levels: 0,
            dump_board: None, bridge_trace_csv: None, modal_pressure_csv: None, receiver_pressure_csv: None,
            pressure_basis_json: None,
            note: None, velocity: None, duration: 6.0,
            sample_rate: 48_000, substeps: 4, modes: linear::MAX_STRING_MODES, pcm_full_scale_pa: 2.0, help: false,
            wav_format: wav_output::Format::default() }
    }
}
impl Options {
    fn parse(args: &[String]) -> Result<Self, String> {
        let mut options = Self::default();
        // The shared loader takes --gas-record and its source packs first.
        let (args, gas_sources) = thermochemical::Options::extract(args)?;
        options.gas_record = gas_sources.record.as_ref().map(|path| path.display().to_string());
        options.gas_sources = gas_sources;
        let mut seen = std::collections::BTreeSet::new();
        let mut args = args.iter();
        while let Some(flag) = args.next() {
            if flag == "--help" || flag == "-h" { options.help = true; continue; }
            if !seen.insert(flag.as_str()) { return Err(format!("duplicate option {flag}")); }
            if flag == "--diagnostic-volume" { options.diagnostic_volume = true; continue; }
            if flag == "--raw-tensions" { options.raw_tensions = true; continue; }
            if flag == "--fast-dampers" { options.fast_dampers = true; continue; }
            if flag == "--rt0425-bridge-contacts" { options.rt0425_bridge_contacts = true; continue; }
            if flag == "--estimated-bridge-contacts" { options.estimated_bridge_contacts = true; continue; }
            if flag == "--rt0425-hammer-stiffness" { options.rt0425_hammer_stiffness = true; continue; }
            if flag == "--rt0425-hammer-dissipation" { options.rt0425_hammer_dissipation = true; continue; }
            if flag == "--felt-pack-hammers" { options.felt_pack_hammers = true; continue; }
            if flag == "--legacy-felt-cards" { options.legacy_felt_cards = true; continue; }
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
                "--string-sections" => options.string_sections = Some(value.clone()),
                "--hammer-footprints" => options.hammer_footprints = Some(value.clone()),
                "--dampers" => options.dampers = Some(value.clone()),
                "--string-stretching" => options.string_stretching = Some(value.clone()),
                "--string-polarization" => options.string_polarization = Some(value.clone()),
                "--cavity" => options.cavity = Some(value.clone()),
                "--felt-pack" => options.felt_pack = Some(value.clone()),
                "--board-panel-pack" => options.board_panel_pack = Some(value.clone()),
                "--board-rib-pack" => options.board_rib_pack = Some(value.clone()),
                "--board-bridge-pack" => options.board_bridge_pack = Some(value.clone()),
                "--rayleigh-loading" => options.rayleigh_loading = Some(rayleigh_loading::Options::parse(value)?),
                "--radiation-ports" => options.radiation_ports = Some(value.parse().map_err(|_| invalid())?),
                "--concert-pitch" => options.concert_pitch = Some(value.parse().map_err(|_| invalid())?),
                "--stretch-tuning" => options.stretch_tuning = Some(value.parse().map_err(|_| invalid())?),
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
                "--modes" => options.modes = geometry::parse_partial_budget(value)?,
                "--pcm-full-scale-pa" => options.pcm_full_scale_pa = value.parse().map_err(|_| invalid())?,
                "--wav-format" => options.wav_format = wav_output::Format::parse(value)?,
                _ => return Err(format!("unknown option {flag}\n{USAGE}")),
            }
        }
        if seen.contains("--wav-format") && options.render.is_none() {
            return Err("--wav-format requires --render".into());
        }
        if options.note.is_some_and(|n| !(21..=108).contains(&n))
            || options.velocity.is_some_and(|v| !v.is_finite() || v <= 0.0 || v > 8.0)
            || !options.duration.is_finite() || !(0.001..=120.0).contains(&options.duration)
            || !(8_000..=192_000).contains(&options.sample_rate)
            || !(1..=16).contains(&options.substeps) || !(1..=linear::MAX_STRING_MODES).contains(&options.modes) {
            return Err("render control outside its finite admitted range".into());
        }
        if [options.concert_pitch, options.stretch_tuning].into_iter().flatten()
            .any(|f| !f.is_finite() || !(430.0..=450.0).contains(&f)) {
            return Err("concert pitch and stretch tuning require finite A4 in 430..450 Hz".into());
        }
        if (options.concert_pitch.is_some() && options.stretch_tuning.is_some())
            || (options.raw_tensions && (options.preset.is_none()
                || options.concert_pitch.is_some() || options.stretch_tuning.is_some())) {
            return Err("select at most one of --concert-pitch, --stretch-tuning and --raw-tensions; --raw-tensions requires a preset".into());
        }
        if options.board.is_some() && (options.board_geometry.is_some() || options.preset.is_some()) {
            return Err("modal board CSV excludes a preset and geometric board; --board-geometry may override a preset board".into());
        }
        if options.preset.as_deref().is_some_and(|p| p != "steinway-d") {
            return Err("unknown piano preset; available: steinway-d".into());
        }
        if ((options.rt0425_bridge_contacts || options.estimated_bridge_contacts) && !options.uses_preset_board())
            || (options.rt0425_bridge_contacts && options.estimated_bridge_contacts)
            || (options.rt0425_hammer_stiffness && (options.preset.is_none()
                || options.hammers.is_some() || options.render.is_none()))
            || (options.rt0425_hammer_dissipation && !options.rt0425_hammer_stiffness) {
            return Err("bridge-contact choices need the preset board and exclude each other; hammer stiffness needs a preset render without --hammers; hammer dissipation also requires source stiffness".into());
        }
        if options.legacy_felt_cards && (options.preset.is_none() || options.felt_pack_hammers
            || options.hammers.is_some() || options.rt0425_hammer_stiffness || options.rt0425_hammer_dissipation) {
            return Err("--legacy-felt-cards selects the preset's former card; it excludes the felt pack, --hammers and the RT-0425 hammer flags".into());
        }
        if options.felt_pack.is_some() && !options.plays_felt_pack() {
            return Err("--felt-pack supplies the felt-pack hammers, which the chosen hammer card does not play".into());
        }
        if (options.board_panel_pack.is_some() || options.board_rib_pack.is_some() || options.board_bridge_pack.is_some())
            && !options.uses_preset_board() {
            return Err("board packs build the preset Model D board; they need --preset without --board or --board-geometry".into());
        }
        if options.felt_pack_hammers && (options.preset.is_none() || options.hammers.is_some()
            || options.render.is_none() || options.rt0425_hammer_stiffness || options.rt0425_hammer_dissipation) {
            return Err("felt-pack hammers need a preset render without --hammers or the RT-0425 hammer flags".into());
        }
        if options.rt0425_string_damping && (options.preset.is_none()
            || options.scale.is_some() || options.string_sections.is_some() || options.render.is_none()) {
            return Err("RT-0425 string damping requires a preset render without a supplied scale".into());
        }
        if options.string_sections.as_deref().is_some_and(|path| path.is_empty() || path.starts_with("--")) {
            return Err("--string-sections requires a nonempty material-data path".into());
        }
        if options.hammers.is_some() && options.render.is_none() {
            return Err("--hammers requires --render; material input is not an export-only option".into());
        }
        if options.hammer_footprints.as_ref().is_some_and(|s| s.trim().is_empty() || options.render.is_none()) {
            return Err("--hammer-footprints requires --render and a complete nonempty specification path".into());
        }
        if options.dampers.as_ref().is_some_and(|s| s.trim().is_empty() || options.render.is_none()) {
            return Err("--dampers requires --render and estimated, point or a nonempty specification path".into());
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
        if options.radiation_ports.is_some_and(|ports|
            options.rayleigh_loading.is_none()
                || !(1..=fs_couple::render::plate::impact::radiation::MAX_PORTS).contains(&ports)) {
            return Err(format!("--radiation-ports requires --rayleigh-loading and an explicit ceiling of 1..{} ports",
                fs_couple::render::plate::impact::radiation::MAX_PORTS));
        }
        if let Some(load) = &options.rayleigh_loading {
            if !geometric || options.render.is_none() || options.diagnostic_volume {
                return Err("--rayleigh-loading requires a flat geometric pressure --render".into());
            }
            load.validate_clock(options.sample_rate, options.sample_rate * options.substeps as u32)?;
            if load.sides == 2 && options.cavity.is_some() {
                return Err("two-sided Rayleigh loading conflicts with a sealed cavity below the board; use one exterior side".into());
            }
        }
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
        if options.cavity.as_ref().is_some_and(|s| s.trim().is_empty()
            || s.starts_with("--") || options.render.is_none() || !geometric) {
            return Err("--cavity requires a complete nonempty specification and a geometric render; modal CSV has no cavity interface".into());
        }
        if options.acoustic_refinement_levels > 3 || (seen.contains("--acoustic-refinement-levels")
            && (!geometric || options.render.is_none() || options.diagnostic_volume)) {
            return Err("--acoustic-refinement-levels requires a flat P1 or edge-cubic geometric pressure render and levels 0..3".into());
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
            options.dampers.as_ref().filter(|s| !matches!(s.as_str(), "estimated" | "point")),
            options.string_stretching.as_ref(), options.string_polarization.as_ref(), options.cavity.as_ref()];
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
        if self.stretch_tuning.is_some() { return None; }
        self.concert_pitch.or_else(||
            (self.preset.is_some() && self.scale.is_none() && self.string_sections.is_none()
                && !self.raw_tensions).then_some(440.0))
    }
    fn tuning_report(&self) -> String {
        if let Some(f) = self.stretch_tuning {
            format!("tensions adjusted to A4={f} Hz, C4..B4 first-partial equal temperament and outward 2:1 stiff-string octave matching; L, mass, EI and contact geometry preserved")
        } else {
            self.tuning_hz().map_or_else(|| "input tensions preserved".to_owned(),
                |f| format!("tensions adjusted to A4={f} Hz first-partial equal temperament; L, mass and EI preserved"))
        }
    }
    /// The preset plays the felt pack's hammers unless another card is chosen:
    /// --hammers, the RT-0425 hammer flags, or --legacy-felt-cards.
    fn plays_felt_pack(&self) -> bool {
        self.felt_pack_hammers || (self.preset.is_some() && self.hammers.is_none()
            && !self.rt0425_hammer_stiffness && !self.rt0425_hammer_dissipation && !self.legacy_felt_cards)
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
    // Bind declared material before retuning or geometry preparation, so every
    // subsequent string/bridge/contact path sees exactly one frozen scale.
    if let Some(path) = options.string_sections.as_deref() {
        let data = std::fs::read_to_string(path).map_err(|e| format!("{path}: {e}"))?;
        scale = string_sections::apply(&data, &scale)?;
    }
    if let Some(reference) = options.stretch_tuning {
        stretch_tuning::tune(&mut scale, reference)?;
    } else if let Some(reference) = options.tuning_hz() {
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
        None if options.preset.as_deref() == Some("steinway-d") =>
            Some(engine::DamperSpecification::estimated(&scale)?),
        None | Some("point") => None,
        Some("estimated") => Some(engine::DamperSpecification::estimated(&scale)?),
        Some(path) => Some(engine::DamperSpecification::load(path, &scale)?),
    };
    let text = options.hammers.as_ref().map(|path| std::fs::read_to_string(path)
        .map_err(|e| format!("{path}: {e}"))).transpose()?;
    let mut piano = prepare_instrument_with_admitted_materials(scale, modes, options, text.as_deref(), stretching, polarization)?;
    if let Some(c) = physical_damping { piano.configure_bare_board_damping(c)?; }
    if let Some(spec) = &dampers { piano.configure_damper_spec(spec)?; }
    if options.fast_dampers && !piano.use_reduced_dampers() {
        return Err("--fast-dampers needs spatial dampers (--dampers estimated, a pad file, or a preset)".into());
    }
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
    // Supplied hammer cards (--hammers text) win over the default felt pack.
    let pack_hammers = (options.plays_felt_pack() && text.is_none()).then(|| {
        let card = selected_felt_card(options)?;
        scale.iter().map(|c| card.hammer(c.midi)).collect::<Result<Vec<_>, String>>()
    }).transpose()?;
    let scale: Vec<geometry::Course> = match &pack_hammers {
        Some(hammers) => scale.into_iter().zip(hammers)
            .map(|(mut c, h)| { c.hammer_mass_kg = h.mass_kg; c }).collect(),
        None => scale,
    };
    let imported = text.map(|text| hammer_materials::read(text,
        &scale.iter().map(|c| c.midi).collect::<Vec<_>>())).transpose()?;
    let footprints=options.hammer_footprints.as_deref().map(|path|
        linear::hammer_footprint::Specification::load(path,&scale)).transpose()?;
    let source_rates=match &pack_hammers {
        Some(hammers)=>Some(scale.iter().zip(hammers).map(|(c,h)|h.rate_per_string(c)).collect::<Vec<_>>()),
        None=>options.rt0425_hammer_dissipation.then(||
            scale.iter().map(|c|steinway_scale::hammer_relaxation_rt0425(c.midi))
                .collect::<Result<Vec<_>,_>>()).transpose()?,
    };
    let mut piano = if options.preset.is_some() {
        let materials = match imported {
            Some(materials) => materials,
            None => match &pack_hammers {
                Some(hammers) => scale.iter().zip(hammers).map(|(c, h)| h.material(c))
                    .collect::<Result<Vec<_>,_>>()?,
                None => scale.iter().map(if options.rt0425_hammer_dissipation {
                    steinway_scale::hammer_material_rt0425_damped
                } else if options.rt0425_hammer_stiffness {
                    steinway_scale::hammer_material_rt0425
                } else { steinway_scale::hammer_material }).collect::<Result<Vec<_>,_>>()?,
            },
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
    prepare_geometric_board_with_source_ports(text, keys, band_hz, equilibrate_mass,
        consistent_mass, edge_cubic_mass, acoustic_refinement_levels, retain_motion, reduction, None)
}
#[allow(clippy::too_many_arguments)]
fn prepare_geometric_board_with_source_ports(text: &str, keys: &[u8], band_hz: f64,
    equilibrate_mass: bool, consistent_mass: bool, edge_cubic_mass: bool,
    acoustic_refinement_levels: usize, retain_motion: bool,
    reduction: Option<&board_geometry::ritz::RitzOptions>,
    frames: Option<&[board_geometry::motion::SourceBridgeFrame]>) -> Result<board_geometry::PreparedBoard, String> {
    if let Some(reduction) = reduction {
        if crowned_board::is_crowned(text) {
            if equilibrate_mass || consistent_mass || edge_cubic_mass || acoustic_refinement_levels != 0 {
                return Err("flat-board mass controls require a flat geometric board".into());
            }
            return crowned_board::CrownedBoard::read(text)?
                .prepare_reduced_with_ports(keys, band_hz, retain_motion, reduction, frames);
        }
        return board_geometry::BoardGeometry::read(text)?
            .with_acoustic_refinement(acoustic_refinement_levels)?
            .prepare_reduced_with_ports(keys, band_hz, retain_motion, equilibrate_mass,
                consistent_mass, edge_cubic_mass, reduction, frames);
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
fn write_board_for_scale(modes: &[linear::BoardMode], scale: &[geometry::Course],
    slopes: Option<&[[[f64; 2]; 88]]>) -> Result<String, String> {
    let all = match slopes {
        Some(slopes) => board::write_with_slopes(modes, slopes)?,
        None => board::write(modes),
    };
    let mut out = String::new();
    for row in all.lines() {
        let admitted = row == board::HEADER || row == board::HEADER_WITH_SLOPES || row.split(',').nth(4)
            .and_then(|key| key.parse::<u8>().ok())
            .is_some_and(|key| scale.iter().any(|c| c.midi == key));
        if admitted { out.push_str(row); out.push('\n'); }
    }
    Ok(out)
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
    render_with_cavity(path, scale, modes, surface, options, stretching, polarization, physical_damping, None)
}
/// The selected air: the committed reference, or --gas-record resolved with
/// the --gas-component packs supplied (a thermochemical record needs them).
fn selected_air(options: &Options) -> Result<ambient::Ambient, String> {
    let loaded = options.gas_sources.load()?;
    let bindings = loaded.bindings();
    ambient::Ambient::selected_with(options.gas_record.as_deref(), loaded.registry(&bindings).as_ref())
}
#[allow(clippy::too_many_arguments)]
fn render_with_cavity(path: &str, scale: Vec<geometry::Course>, modes: &[linear::BoardMode],
    surface: Option<&[board_geometry::SurfaceSample]>, options: &Options,
    stretching: Option<&linear::string_stretching::Specification>,
    polarization: Option<&string_polarization::Prepared>, physical_damping: Option<&[f64]>,
    cavity: Option<&cavity::Projected>) -> Result<(), String> {
    if options.cavity.is_some() && cavity.is_none() {
        return Err("supplied cavity was not admitted and projected before rendering".into());
    }
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
    let mut piano = prepare_instrument_with_board_damping(scale, modes, options, stretching, polarization, physical_damping)?;
    if let Some(cavity) = cavity {
        piano.configure_cavity(&cavity.loaded(&piano.bank)?)?;
        println!("{}", cavity.report());
    }
    // The exterior air is data: the committed reference record or --gas-record.
    let ambient = selected_air(options)?;
    println!("{}", ambient.receipt());
    if let Some(load_options) = &options.rayleigh_loading {
        let load = rayleigh_loading::prepare(
            surface.ok_or("Rayleigh loading requires the played geometric surface")?,
            &piano.bank, load_options, ambient.medium(), piano.sample_rate(),
            options.radiation_ports)?;
        if let Some(basis) = &load.basis {
            piano.configure_projected_radiation(&load.model, basis)?;
        } else { piano.configure_radiation(&load.model)?; }
        println!("{}", load.report);
    }
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
        println!("String extension: {path}; {count} geometric speaking/duplex channels; {} retained longitudinal coordinates; same-tick contact and reciprocal transverse bridge work; the speaking end's axial force drives the bridge through the declared downbearing/bridge-cap projection. Zero channels means explicit all-linear selection. EA, wave retention and drag are supplied; slope/strain/iteration limits refuse rather than clamp. No real-time or calibration claim; see STRING_STRETCHING.md.",
            piano.bank.longitudinal_mode_count());
    }
    if let Some(path)=&options.hammer_footprints {
        println!("Hammer contact geometry: {path}; {} independent felt sites, unchanged total course area and hammer mass. No inferred/calibrated face width; see HAMMER_FOOTPRINTS.md.",
            piano.hammer_contact_count());
    }
    if let Some((strings,cells)) = piano.damper_resolution() {
        println!("Spatial viscous dampers: {} speaking-string pads, {} quadrature stations; source {}. No measured pad/action or real-time claim; see DAMPERS.md.",
            strings, cells, options.dampers.as_deref().unwrap_or("estimated preset defaults"));
    }
    if let Some((strings,sites)) = piano.felt_damper_resolution() {
        println!("Finite-gap felt dampers: {strings} speaking-string pads, {sites} independent contact sites; source {}. Supplied finite-speed travel, wool-felt crush and Prony memory, coupled on the mechanics clock; see DAMPERS.md.",
            options.dampers.as_deref().unwrap_or("supplied"));
    }
    let surface = if options.diagnostic_volume { None } else { surface };
    let left=options.microphone.unwrap_or([0.675,1.0,1.0]);
    let mut stream = match options.microphone_right {
        Some(right)=>audio::AudioStream::new_stereo(piano,score,
            surface.ok_or("stereo pressure requires a geometric soundboard surface")?,
            [left,right],ambient.medium())?,
        None=>audio::AudioStream::new(piano,score,surface,left,
            ambient.medium(),options.observer_gain)?,
    };
    // A supplied record must reach the sound: the microphone's medium,
    // Rayleigh loading, a gas,ambient cavity or the preset board's wood.
    if options.gas_record.is_some() && stream.microphone().is_none() && options.rayleigh_loading.is_none()
        && !cavity.is_some_and(cavity::Projected::reads_ambient) && !options.uses_preset_board() {
        return Err("--gas-record has no consumer: the diagnostic volume observer hears no gas, and no Rayleigh loading, gas,ambient cavity or preset board reads it".into());
    }
    // the bridge trace reads the instrument after every frame
    if bridge_row.is_some() { stream.disable_accelerator(); }
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
    stream.set_horizon(u64::from(count));
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
    if options.wav_format != wav_output::Format::Float32 {
        check_pcm_headroom(peak, options.pcm_full_scale_pa)?;
    }
    let wav = options.wav_format.encode(&pressure, rate, channels as u16, options.pcm_full_scale_pa)?;
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
        println!("Computed half-space pressure in Pa; {} full scale {} Pa, no peak normalization. Infinite baffle; no room/lid scattering or measured-SPL calibration. Air backreaction: {}.", options.wav_format.label(), options.pcm_full_scale_pa,
            if options.rayleigh_loading.is_some() { "explicit passive Rayleigh load" } else { "absent" });
    } else {
        println!("Diagnostic volume-velocity observer, gain {} Pa/(m^3/s); no peak normalization or calibrated SPL claim.", options.observer_gain);
    }
    let piano = stream.instrument();
    println!("{} string modes; {} board modes; {} above-band duplex segments omitted from dynamic retention (static attachment retained); {} string-direction segments omit in-band partials at the selected count ceiling {}.",
        piano.bank.modes.len(), piano.bank.board_count, piano.bank.omitted_duplex_modes,
        piano.bank.capped_string_segments, options.modes);
    println!("{seconds:.6} s, {channels} channel(s) rendered in {elapsed:.6} s; wall/audio ratio {:.4}; peak {peak:.6} Pa-equivalent; {} output, no clipped samples.", elapsed / seconds, options.wav_format.label());
    if options.wav_format == wav_output::Format::Float32 && peak > options.pcm_full_scale_pa {
        println!("Float WAV preserves pressure above nominal digital full scale; attenuate before fixed-point playback.");
    }
    println!("Input {:.9} J; stored {:.9} J; component losses {:.9} J; closure {:.3e} J; worst substep defect {:.3e} J.",
        piano.accounting.input_work_j, piano.energy_j(), piano.accounting.dissipated_j(),
        piano.accounting.input_work_j - piano.energy_j() - piano.accounting.dissipated_j(), piano.accounting.max_balance_error_j);
    println!("Felt loss {:.9} J, including {:.9} J time-dependent relaxation; shank damping {:.9} J.",
        piano.accounting.felt_loss_j, piano.accounting.felt_relaxation_loss_j, piano.accounting.shank_loss_j);
    if options.rayleigh_loading.is_some() {
        println!("Radiation storage {:.9} J; passive radiation loss {:.9} J, each included once in the combined balance.",
            piano.radiation_energy_j(), piano.accounting.radiation_loss_j);
    }
    if piano.has_cavity() {
        println!("Cavity storage {:.9} J; explicit air momentum loss {:.9} J, each included once in the combined balance.",
            piano.cavity_energy_j(), piano.accounting.cavity_loss_j);
    }
    if piano.damper_resolution().is_some() {
        println!("Spatial damper loss {:.9} J, already included in component losses; no output-envelope damping.",
            piano.accounting.damper_loss_j);
    }
    if piano.felt_damper_resolution().is_some() {
        println!("Felt damper storage {:.9} J; loss {:.9} J, including {:.9} J Prony relaxation; signed actuator work {:.9} J. Storage, loss and work are already included once in the combined totals.",
            piano.damper_energy_j(), piano.accounting.damper_loss_j,
            piano.accounting.damper_relaxation_loss_j, piano.accounting.damper_actuator_work_j);
    }
    Ok(())
}

/// A pack file as admitted bytes; a missing, truncated or corrupt file refuses.
fn pack_file(path: &str) -> Result<fs_matdb::NormalizedPack, String> {
    let bytes = std::fs::read(path).map_err(|e| format!("{path}: {e}"))?;
    fs_matdb::NormalizedPack::from_bytes(&bytes).map_err(|e| format!("{path}: {e}"))
}

/// The --felt-pack-hammers card: --felt-pack FILE, or the built-in Stulov set.
fn selected_felt_card(options: &Options) -> Result<hammer_felt::FeltCard, String> {
    match &options.felt_pack {
        Some(path) => hammer_felt::FeltCard::from_pack(&pack_file(path)?),
        None => hammer_felt::stulov_abel(),
    }
}

/// The preset board's materials: the command line's packs (or the built-in
/// Model D packs), in the wood moisture of the SAME air the piano radiates into.
fn preset_board_materials(options: &Options, ambient: &ambient::Ambient) -> Result<steinway_d::BoardMaterials, String> {
    let environment = steinway_d::BoardEnvironment::from_gas(&ambient.gas)?;
    let pack = |path: &Option<String>, built_in: &[u8]| match path {
        Some(path) => pack_file(path),
        None => steinway_d::load_pack(built_in),
    };
    steinway_d::model_d_materials_from(&pack(&options.board_panel_pack, steinway_d::SITKA_SPRUCE)?,
        &pack(&options.board_rib_pack, steinway_d::SUGAR_PINE)?,
        &pack(&options.board_bridge_pack, steinway_d::HARD_MAPLE)?, &environment)
}

fn run() -> Result<(), String> { run_with(&std::env::args().skip(1).collect::<Vec<_>>()) }
fn run_with(args: &[String]) -> Result<(), String> {
    let options = Options::parse(args)?;
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
    // A `gas,ambient` cavity holds the same selected air as the exterior.
    let air = selected_air(&options)?;
    let cavity = options.cavity.as_deref().map(cavity::Specification::load).transpose()?
        .map(|cavity| cavity.in_air(&air.gas));
    let scale_source = options.string_sections.as_deref().or(options.scale.as_deref()).unwrap_or(if options.preset.is_some() {
        "RT-0425 Appendix A wrapped-string MODEL: 84 published courses plus four estimated extensions"
    } else { "ESTIMATED demonstration" });
    let tuning_source = options.tuning_report();
    let preset = options.uses_preset_board().then(|| {
        let materials = preset_board_materials(&options, &air)?;
        steinway_d::build_with_materials(options.mesh_divisions, !options.estimated_bridge_contacts, &materials)
    }).transpose()?;
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
    if options.rayleigh_loading.is_some()
        && geometry_text.as_deref().is_some_and(crowned_board::is_crowned) {
        return Err("Rayleigh loading requires a flat FSB board; a crowned shell cannot be flattened for physical backreaction".into());
    }
    // Gradient of every mode shape at each key's bridge station, for a --dump-board of a flat geometric board.
    let mut board_slopes: Option<Vec<[[f64; 2]; 88]>> = None;
    let (modes, board_source, surface, polarization, physical_damping, cavity) = if let Some(text) = &geometry_text {
        let start = std::time::Instant::now();
        let source_ports = polarization.as_ref().map(|spec| spec.source_ports(&scale)).transpose()?;
        let want_slopes = options.dump_board.is_some() && options.render.is_none() && !crowned_board::is_crowned(text);
        let prepared = prepare_geometric_board_with_source_ports(text, &scale.iter().map(|c| c.midi).collect::<Vec<_>>(),
            options.board_band_hz, options.equilibrate_board_mass, options.consistent_board_mass,
            options.edge_cubic_board_mass, options.acoustic_refinement_levels,
            polarization.is_some() || cavity.is_some() || want_slopes,
            options.board_reduction.as_ref(), source_ports.as_deref())?;
        let projected = polarization.as_ref().map(|spec|
            spec.project(&scale, &prepared.modes, prepared.motion.as_ref())).transpose()?;
        let cavity = cavity.as_ref().map(|spec| spec.project(&prepared)).transpose()?;
        if want_slopes {
            let motion = prepared.motion.as_ref().ok_or("slope columns need the retained board motion")?;
            let keys: Vec<u8> = scale.iter().map(|c| c.midi).collect();
            board_slopes = Some(board_geometry::BoardGeometry::read(text)?.bridge_slopes(motion, &keys)?);
        }
        let model_name = if crowned_board::is_crowned(text) { "Crowned shell" } else { "Flat plate" };
        println!("{model_name} {:.6} m^2, {:.6} kg (panel+ribs/bridges), {} free DOFs, {} retained modes; source band (0,{}] Hz; preparation {:.6} s.",
            prepared.area_m2, prepared.mass_kg, prepared.free_dofs, prepared.modes.len(),
            options.board_band_hz, start.elapsed().as_secs_f64());
        if let Some(report) = &prepared.reduction {
            println!("Bridge-informed reduction: {} certified source modes -> {} retained coordinates; {} exact low modes; {} normalized static/harmonic snapshots; requested seeds {:?} Hz; maximum displacement projection residual {:.9e}. This residual does not certify transfer, acoustic or mesh convergence.",
                report.source_modes, prepared.modes.len(), report.protected_low_modes,
                report.snapshot_count, report.sample_hz, report.max_relative_snapshot_error);
            if let Some(check) = &report.adaptive_response {
                println!("Adaptive bridge response: maximum relative complex driving-point error {:.9e}; maximum normalized complex bridge-transfer error {:.9e} <= {:.9e}; {} frequencies, {} force ports, {} enrichment rounds. Finite-grid source-model comparison; no continuous-band, acoustic or mesh-convergence certificate.",
                    check.max_relative_driving_point_error, check.max_normalized_transfer_error, check.tolerance,
                    check.evaluation_hz.len(), check.port_count, check.enrichment_rounds);
            }
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
            Some(prepared.surface), projected, prepared.physical_damping, cavity)
    } else {
        (load_board(board_text.as_deref(), &scale)?, options.board.as_deref()
            .unwrap_or("AUTHORED illustrative modes; not measured Steinway geometry").to_owned(), None, None, None, None)
    };
    if modes.iter().any(|m| m.frequency_hz >= 0.45 * f64::from(options.sample_rate)) {
        return Err("soundboard mode at/above output retention ceiling; use an explicitly reduced board".into());
    }
    study_key(&scale, options.note)?;
    println!("String scale: {scale_source}; {tuning_source}.");
    println!("Soundboard: {board_source}.");
    if let Some(path) = &options.hammers {
        println!("Hammer materials supplied by {path}; no automatic coupon-fit or calibration claim. Scale mass/patch geometry unchanged.");
    } else if options.plays_felt_pack() {
        let card = selected_felt_card(&options)?;
        println!("Hammers from felt pack {} {}: Stulov's 88 measured Abel hammers, K u^p + R d(u^p)/dt per hammer with its mass.",
            card.pack_id(), card.content_hash());
    } else if options.preset.is_some() {
        println!("Per-key source-derived hammer loading envelopes{}; {}.",
            if options.rt0425_hammer_stiffness { " with RT-0425 stiffness per unison string" } else { " with legacy stiffness divided across the unison" },
            if options.rt0425_hammer_dissipation { "RT-0425 R_H power-law rate loss, no estimated crush or Prony terms" }
            else { "estimated unloading/crush and tangent-scaled Prony relaxation, not RT-0425 R_H" });
    }
    if options.preset.is_some() {
        println!("Published shank geometry -> rigid rotation + bending; reciprocal jack port and 1.5 mm let-off. Linearized action fragment; damping/backcheck estimated.");
    }
    let material_loss_courses = scale.iter().filter(|course|
        course.bonded_section.and_then(|section| section.bending_loss()).is_some()).count();
    if material_loss_courses > 0 {
        println!("String damping: material Kelvin-Voigt bending on {material_loss_courses}/{} courses, checked at every retained modal frequency; separately authored common velocity drag. Other courses retain the demonstration bending spectrum.", scale.len());
    } else {
        println!("String damping: {}.", if options.rt0425_string_damping {
            "RT-0425 per-key R_u and eta_u projected onto scalar stiff-string modes; four end keys extrapolated"
        } else { "estimated common rate and authored Maxwell bending spectrum" });
    }
    println!("Source authority belongs to the inputs, not the model name; imported files are not independently certified measurements.");
    if let Some(path) = &options.dump_scale {
        write_fresh_output(path, format!("# Source: {scale_source}; {tuning_source}\n{}", geometry::write_scale(&scale)).as_bytes())?;
    }
    if let Some(path) = &options.dump_board {
        write_fresh_output(path, format!("# Source: {board_source}\n{}",
            write_board_for_scale(&modes, &scale, board_slopes.as_deref())?).as_bytes())?;
    }
    if let Some(path) = &options.render {
        return render_with_cavity(path, scale, &modes, surface.as_deref(), &options,
            stretching.as_ref(), polarization.as_ref(), physical_damping.as_deref(), cavity.as_ref());
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

/// The command line, for a host that compiles these sources as a module and installs an accelerator
/// first (apps/instrument-lab/render-gpu-piano).
#[allow(dead_code)]
pub(crate) fn cli_main(args: &[String]) -> Result<(), String> { run_with(args) }

#[cfg(test)]
#[path = "rayleigh_loading_render_tests.rs"]
mod rayleigh_loading_render_tests;

#[cfg(test)]
#[path = "wav_output_render_tests.rs"]
mod wav_output_render_tests;

#[cfg(test)]
#[path = "ambient_render_tests.rs"]
mod ambient_render_tests;

#[cfg(test)]
#[path = "cavity_render_tests.rs"]
mod cavity_render_tests;

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
#[path = "longitudinal_render_tests.rs"]
mod longitudinal_render_tests;

#[cfg(test)]
#[path = "felt_damper_render_tests.rs"]
mod felt_damper_render_tests;

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
    fn wav_format_requires_a_render_and_an_explicit_supported_encoding() {
        assert_eq!(options(&[]).unwrap().wav_format, wav_output::Format::Pcm16);
        for (name, format) in [("pcm16", wav_output::Format::Pcm16),
            ("pcm24", wav_output::Format::Pcm24), ("float32", wav_output::Format::Float32)] {
            assert_eq!(options(&["--render", "piano.wav", "--wav-format", name]).unwrap().wav_format, format);
            assert!(options(&["--wav-format", name]).is_err());
        }
        for args in [vec!["--render", "piano.wav", "--wav-format"],
            vec!["--render", "piano.wav", "--wav-format", "float64"],
            vec!["--render", "piano.wav", "--wav-format", "pcm24", "--wav-format", "float32"]] {
            assert!(options(&args).is_err());
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
        assert!(options(&["--render", "piano.wav", "--dampers", "point"]).is_ok());
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
    fn preset_dampers_reach_nodal_partials_after_key_and_pedal_release() {
        // A free bridge isolates the twentieth speaking partial from board
        // loss. Its node at x/L=.35 defeats the original point damper.
        let c = geometry::Course { midi:69, unison:1, length_m:1.0,
            linear_density_kg_m:0.02, tension_n:100.0, flexural_rigidity_nm2:0.0,
            bonded_section:None,
            strike_fraction:0.12, hammer_mass_kg:0.008, felt_thickness_m:0.008,
            felt_area_m2:0.0001, duplex_length_m:0.0, detune_cents:0.0 };
        let board = [linear::BoardMode { frequency_hz:100.0, damping_ratio:0.0,
            bridge:[0.0;88], volume:0.1 }];
        let build = |selection: Option<&str>| {
            let o = Options { preset:Some("steinway-d".into()), modes:24,
                dampers:selection.map(str::to_owned), ..Options::default() };
            let mut p = prepare_instrument(vec![c], &board, &o).unwrap();
            p.bank.v[19] = 0.2;
            p.accounting.input_work_j = p.energy_j();
            p.silent_key_down(69).unwrap();
            p
        };
        let mut pad = build(None); let mut point = build(Some("point"));
        assert!(pad.damper_resolution().is_some());
        assert_eq!(point.damper_resolution(), None);
        for held in [true, false] {
            if !held {
                for p in [&mut pad, &mut point] {
                    p.set_sustain(1.0).unwrap(); p.note_off(69).unwrap();
                }
            }
            for _ in 0..64 { pad.step().unwrap(); point.step().unwrap(); }
            assert_eq!(pad.bank.q, point.bank.q);
            assert_eq!(pad.bank.v, point.bank.v);
            assert_eq!(pad.accounting.damper_loss_j, 0.0);
        }
        pad.set_sustain(0.0).unwrap(); point.set_sustain(0.0).unwrap();
        for _ in 0..128 { pad.step().unwrap(); point.step().unwrap(); }
        assert!(pad.accounting.damper_loss_j > 1e-6);
        assert!(point.accounting.damper_loss_j.abs() < 1e-20);
        assert!(pad.bank.energy() < point.bank.energy());
        for p in [&pad, &point] {
            assert!((p.accounting.input_work_j-p.energy_j()-p.accounting.dissipated_j()).abs() < 1e-10);
        }
    }
    #[test]
    fn explicit_free_dampers_override_the_preset_and_nonpreset_keeps_point_damping() {
        let c = geometry::demonstration_scale().unwrap()[48];
        let mut o = Options { modes:24, ..Options::default() };
        assert_eq!(prepare_instrument(vec![c], &board::demonstration(), &o)
            .unwrap().damper_resolution(), None);
        let dir = output_fixture(); let free = dir.join("free.fspd");
        write_fresh_output(free.to_str().unwrap(),
            b"frankensim-piano-dampers-v1\nfree,69\n").unwrap();
        o.preset = Some("steinway-d".into());
        o.dampers = Some(free.to_str().unwrap().to_owned());
        let mut p = prepare_instrument(vec![c], &board::demonstration(), &o).unwrap();
        assert_eq!(p.damper_resolution(), Some((0,0)));
        p.bank.v[0] = 0.01;
        for _ in 0..64 { p.step().unwrap(); }
        assert_eq!(p.accounting.damper_loss_j, 0.0);
        assert!(p.bank.energy() > 0.0);
    }
    #[test]
    fn rt0425_corrections_are_explicit_preset_inputs() {
        let o = options(&["--preset", "steinway-d", "--render", "piano.wav",
            "--rt0425-bridge-contacts", "--rt0425-hammer-stiffness"]).unwrap();
        assert!(o.rt0425_bridge_contacts && o.rt0425_hammer_stiffness);
        let plain = options(&["--preset", "steinway-d", "--render", "piano.wav"]).unwrap();
        assert!(!plain.rt0425_bridge_contacts && !plain.estimated_bridge_contacts);
        assert!(options(&["--preset", "steinway-d", "--render", "piano.wav", "--estimated-bridge-contacts"])
            .unwrap().estimated_bridge_contacts);
        for args in [
            vec!["--rt0425-bridge-contacts"],
            vec!["--estimated-bridge-contacts"],
            vec!["--preset", "steinway-d", "--render", "piano.wav", "--rt0425-bridge-contacts", "--estimated-bridge-contacts"],
            vec!["--preset", "steinway-d", "--board-geometry", "board.fsb", "--estimated-bridge-contacts"],
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
    fn board_and_felt_packs_are_command_line_data_in_the_selected_air() {
        let stamp = std::time::SystemTime::now().duration_since(std::time::UNIX_EPOCH).unwrap().as_nanos();
        let directory = std::env::temp_dir().join(format!("frankensim-piano-packs-{}-{stamp}", std::process::id()));
        std::fs::create_dir(&directory).unwrap();
        let write = |name: &str, bytes: &[u8]| {
            let path = directory.join(name);
            std::fs::write(&path, bytes).unwrap();
            path.to_str().unwrap().to_owned()
        };
        let engelmann = write("engelmann.fsmatpk", include_bytes!("packs/spruce-engelmann-fpl-gtr282.fsmatpk"));
        let truncated = write("truncated.fsmatpk", &steinway_d::SITKA_SPRUCE[..100]);
        let base = ["--preset", "steinway-d", "--render", "piano.wav"];
        let with = |extra: &[&str]| options(&base.iter().copied().chain(extra.iter().copied()).collect::<Vec<_>>());
        let reference = ambient::Ambient::reference();
        // No pack options: the built-in board, bit for bit, in the committed air.
        let built_in = preset_board_materials(&with(&[]).unwrap(), &reference).unwrap();
        assert_eq!(built_in, steinway_d::reference_materials().unwrap());
        // A panel pack on the command line is a data swap.
        let swapped = preset_board_materials(&with(&["--board-panel-pack", engelmann.as_str()]).unwrap(), &reference).unwrap();
        assert_eq!(swapped.panel.pack.pack_id, "spruce-engelmann-fpl-gtr282");
        assert_eq!(swapped.ribs, built_in.ribs);
        assert_ne!(swapped.panel, built_in.panel);
        for bad in [truncated.as_str(), "/nonexistent/panel.fsmatpk"] {
            let error = preset_board_materials(&with(&["--board-panel-pack", bad]).unwrap(), &reference).unwrap_err();
            assert!(error.contains(bad), "{error}");
        }
        // The wood reads the same air the piano radiates into: drier air, drier wood;
        // dry air is outside the isotherm and refuses rather than inventing a moisture.
        let drier = ambient::Ambient { gas: ambient::parse("293.15 101325 0.40").unwrap(), source: ambient::Source::CommandLine };
        let dried = preset_board_materials(&with(&[]).unwrap(), &drier).unwrap();
        assert!(dried.environment.moisture_content < built_in.environment.moisture_content);
        assert!(dried.panel.e_l > built_in.panel.e_l && dried.panel.density < built_in.panel.density);
        let dry = ambient::Ambient { gas: ambient::parse("293.15 101325 0").unwrap(), source: ambient::Source::CommandLine };
        assert!(preset_board_materials(&with(&[]).unwrap(), &dry).unwrap_err().contains("wood-emc-hailwood-horrobin-fpl-gtr282"));
        // Board packs need the preset board; --felt-pack needs the felt-pack card in play.
        assert!(options(&["--board-panel-pack", engelmann.as_str()]).is_err());
        assert!(with(&["--legacy-felt-cards", "--felt-pack", engelmann.as_str()]).is_err());
        let felt = write("felt.fsmatpk", include_bytes!("packs/hammer-felt-abel-stulov2005.fsmatpk"));
        // The felt pack is the preset's default card, so --felt-pack alone supplies it.
        assert!(selected_felt_card(&with(&["--felt-pack", felt.as_str()]).unwrap()).is_ok());
        let card = selected_felt_card(&with(&["--felt-pack-hammers", "--felt-pack", felt.as_str()]).unwrap()).unwrap();
        assert_eq!(card.pack_id(), "hammer-felt-abel-stulov2005");
        assert!(selected_felt_card(&with(&["--felt-pack-hammers", "--felt-pack", truncated.as_str()]).unwrap()).is_err());
        assert!(selected_felt_card(&with(&["--felt-pack-hammers", "--felt-pack", engelmann.as_str()]).unwrap()).unwrap().hammer(69).is_err());
    }
    #[test]
    fn the_preset_plays_the_felt_pack_by_default_and_the_legacy_card_on_request() {
        let base = ["--preset", "steinway-d", "--render", "piano.wav"];
        let with = |extra: &[&str]| options(&base.iter().copied().chain(extra.iter().copied()).collect::<Vec<_>>());
        let (default, legacy) = (with(&[]).unwrap(), with(&["--legacy-felt-cards"]).unwrap());
        assert!(default.plays_felt_pack() && !legacy.plays_felt_pack());
        assert!(with(&["--felt-pack-hammers"]).unwrap().plays_felt_pack());
        assert!(!with(&["--rt0425-hammer-stiffness"]).unwrap().plays_felt_pack());
        for extra in [
            vec!["--legacy-felt-cards", "--felt-pack-hammers"],
            vec!["--legacy-felt-cards", "--rt0425-hammer-stiffness"],
            vec!["--legacy-felt-cards", "--hammers", "felt.fsh"],
        ] { assert!(with(&extra).is_err(), "accepted {extra:?}"); }
        assert!(options(&["--legacy-felt-cards"]).is_err());
        // The golden bump in numbers: the same strike, felt pack against the former card.
        let modes = board::demonstration();
        for midi in [36u8, 60, 84] {
            let course = steinway_scale::courses().unwrap()[usize::from(midi - 21)];
            let mut pack = prepare_instrument_with_hammers(vec![course], &modes, &default, None).unwrap();
            let mut old = prepare_instrument_with_hammers(vec![course], &modes, &legacy, None).unwrap();
            pack.note_on(midi, 2.5).unwrap();
            old.note_on(midi, 2.5).unwrap();
            let (mut a, mut b) = (Vec::new(), Vec::new());
            for _ in 0..4800 {
                a.push(pack.step().unwrap());
                b.push(old.step().unwrap());
            }
            let rms = |v: &[f64]| (v.iter().map(|x| x * x).sum::<f64>() / v.len() as f64).sqrt();
            let difference: Vec<f64> = a.iter().zip(&b).map(|(x, y)| x - y).collect();
            eprintln!(
                "key {midi} at 2.5 m/s, 0.1 s: output rms felt pack {:.4e}, legacy {:.4e}; difference {:.4e} ({:.1} % of legacy); felt loss felt pack {:.4e} J (rate {:.4e} J), legacy {:.4e} J",
                rms(&a), rms(&b), rms(&difference), 100.0 * rms(&difference) / rms(&b),
                pack.accounting.felt_loss_j, pack.accounting.felt_relaxation_loss_j, old.accounting.felt_loss_j
            );
            assert!(rms(&difference) > 0.0, "key {midi}: the felt pack changed nothing");
            let balance = pack.accounting.input_work_j - pack.energy_j() - pack.accounting.dissipated_j();
            assert!(balance.abs() < 1e-7, "key {midi}: {balance:e}");
        }
    }
    #[test]
    fn felt_pack_hammers_are_an_explicit_preset_input_that_closes_its_work() {
        let o = options(&["--preset", "steinway-d", "--render", "piano.wav", "--felt-pack-hammers"]).unwrap();
        assert!(o.felt_pack_hammers);
        for args in [
            vec!["--felt-pack-hammers"],
            vec!["--preset", "steinway-d", "--felt-pack-hammers"],
            vec!["--preset", "steinway-d", "--render", "piano.wav", "--hammers", "felt.fsh", "--felt-pack-hammers"],
            vec!["--preset", "steinway-d", "--render", "piano.wav", "--rt0425-hammer-stiffness", "--felt-pack-hammers"],
            vec!["--preset", "steinway-d", "--render", "piano.wav", "--rt0425-hammer-stiffness",
                "--rt0425-hammer-dissipation", "--felt-pack-hammers"],
        ] { assert!(options(&args).is_err(), "accepted {args:?}"); }
        let source = options(&["--preset", "steinway-d", "--render", "piano.wav",
            "--rt0425-hammer-stiffness", "--rt0425-hammer-dissipation"]).unwrap();
        let course = steinway_scale::courses().unwrap()[48];
        let modes = board::demonstration();
        let mut pack = prepare_instrument_with_hammers(vec![course], &modes, &o, None).unwrap();
        let mut rt0425 = prepare_instrument_with_hammers(vec![course], &modes, &source, None).unwrap();
        pack.note_on(69, 2.0).unwrap();
        rt0425.note_on(69, 2.0).unwrap();
        let mut changed = false;
        for _ in 0..2400 {
            changed |= pack.step().unwrap() != rt0425.step().unwrap();
        }
        assert!(changed, "the felt pack changed nothing");
        assert!(pack.accounting.felt_relaxation_loss_j > 0.0);
        let balance = pack.accounting.input_work_j - pack.energy_j() - pack.accounting.dissipated_j();
        assert!(balance.abs() < 1e-7, "{balance:e}");
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
    fn default_and_auto_follow_the_output_band_with_a_numeric_preview_control() {
        let ordinary=options(&[]).unwrap();
        let automatic=options(&["--modes","auto"]).unwrap();
        assert_eq!(ordinary.modes,linear::MAX_STRING_MODES);
        assert_eq!(automatic.modes,ordinary.modes);
        let mut c=geometry::demonstration_scale().unwrap()[0];
        c.unison=1;c.duplex_length_m=0.0;
        let board=board::demonstration();
        let full=prepare_instrument(vec![c],&board,&ordinary).unwrap();
        let preview=prepare_instrument(vec![c],&board,
            &options(&["--modes","24"]).unwrap()).unwrap();
        assert_eq!(preview.bank.modes.len(),24);
        assert_eq!(preview.bank.capped_string_segments,1);
        assert!(full.bank.modes.len()>128);
        assert_eq!(full.bank.capped_string_segments,0);
        assert!(full.bank.modes.iter().all(|m|m.omega/std::f64::consts::TAU<=21_600.0));
    }
    #[test]
    fn acoustic_refinement_is_an_explicit_bounded_pressure_control() {
        let ordinary=options(&["--preset","steinway-d","--render","p.wav"]).unwrap();
        assert_eq!(ordinary.acoustic_refinement_levels,0);
        let refined=options(&["--preset","steinway-d","--render","p.wav",
            "--acoustic-refinement-levels","2"]).unwrap();
        assert_eq!(refined.acoustic_refinement_levels,2);
        let cubic=options(&["--preset","steinway-d","--render","p.wav",
            "--acoustic-refinement-levels","2","--edge-cubic-board-mass"]).unwrap();
        assert_eq!(cubic.acoustic_refinement_levels,2);
        assert!(cubic.edge_cubic_board_mass);
        for args in [vec!["--acoustic-refinement-levels","1"],
            vec!["--render","p.wav","--acoustic-refinement-levels","1"],
            vec!["--preset","steinway-d","--render","p.wav","--acoustic-refinement-levels","4"],
            vec!["--preset","steinway-d","--render","p.wav","--acoustic-refinement-levels","1","--diagnostic-volume"],
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
        o.dampers = Some(if spatial { "estimated" } else { "point" }.to_owned());
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
        let text = write_board_for_scale(&modes, &scale, None).unwrap();
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
    fn board_export_with_slopes_keeps_the_admitted_subset_and_the_slope_header() {
        let scale = vec![geometry::demonstration_scale().unwrap()[48]];
        let modes = board::demonstration();
        let slopes: Vec<[[f64; 2]; 88]> = (0..modes.len())
            .map(|m| std::array::from_fn(|k| [m as f64 + 0.5, k as f64])).collect();
        let text = write_board_for_scale(&modes, &scale, Some(&slopes)).unwrap();
        assert!(text.starts_with(board::HEADER_WITH_SLOPES));
        assert_eq!(text.lines().count(), modes.len() + 1);
        let row = text.lines().nth(1).unwrap().split(',').collect::<Vec<_>>();
        assert_eq!((row[4], row[6], row[7]), ("69", "5.00000000000000000e-1", "4.80000000000000000e1"));
        assert_eq!(board::read(&text, &[69]).unwrap().len(), modes.len());
        assert!(write_board_for_scale(&modes, &scale, Some(&slopes[..1])).is_err());
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
        // 48 from the FPL-GTR-282 board packs; the former hand-typed constants gave 43.
        assert_eq!(prepared.modes.len(), 48);
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
    fn stretch_selection_preserves_geometry_and_reaches_the_played_string_bank() {
        let o=options(&["--preset","steinway-d","--stretch-tuning","442"]).unwrap();
        assert_eq!(o.tuning_hz(),None);
        assert!(o.tuning_report().contains("2:1"));
        let source=steinway_scale::courses().unwrap();
        let tuned=selected_scale(None,&o).unwrap();
        let equal=selected_scale(None,&options(&["--preset","steinway-d","--concert-pitch","442"]).unwrap()).unwrap();
        let mut before=0.0_f64;let mut after=0.0_f64;
        for (a,b) in source.iter().zip(&tuned) {
            assert_eq!(*a,geometry::Course {tension_n:a.tension_n,..*b});
        }
        for index in 0..76 {
            let mismatch=|s:&[geometry::Course]| (s[index].partial_hz(2,s[index].tension_n)
                /s[index+12].partial_hz(1,s[index+12].tension_n)-1.0).abs();
            before=before.max(mismatch(&equal));after=after.max(mismatch(&tuned));
        }
        assert!(before>1e-4 && after<1e-10);
        // Prepared playback receives the new physical tension, not a frequency
        // override. Verify every retained speaking partial, including unisons.
        let played=vec![tuned[27],tuned[48],tuned[63]];
        let mut render_options=o;render_options.modes=4;
        let piano=prepare_instrument(played.clone(),&board::demonstration(),&render_options).unwrap();
        for port in piano.bank.strings.iter().filter(|p|!p.duplex) {
            let c=&played[port.course];
            let cents=(port.member as f64-0.5*(c.unison-1) as f64)*c.detune_cents;
            let tension=c.tension_at_cents(cents).unwrap();
            for (partial,index) in port.modes.clone().enumerate() {
                let actual=piano.bank.modes[index].omega/std::f64::consts::TAU;
                let expected=c.partial_hz(partial+1,tension);
                assert!((actual/expected-1.0).abs()<1e-11);
            }
        }
        assert!((tuned[48].partial_hz(1,tuned[48].tension_n)/442.0-1.0).abs()<1e-11);
        let imported=options(&["--preset","steinway-d","--scale","source.csv","--stretch-tuning","442"]).unwrap();
        assert_eq!(selected_scale(Some(&geometry::write_scale(&source)),&imported).unwrap(),tuned);
    }
    #[test]
    fn stretch_selection_rejects_conflicting_or_invalid_tuning_before_preparation() {
        for args in [vec!["--stretch-tuning"],vec!["--stretch-tuning","NaN"],
            vec!["--stretch-tuning","429.9"],vec!["--stretch-tuning","450.1"],
            vec!["--stretch-tuning","440","--stretch-tuning","442"],
            vec!["--stretch-tuning","440","--concert-pitch","440"],
            vec!["--preset","steinway-d","--stretch-tuning","440","--raw-tensions"]] {
            assert!(options(&args).is_err(),"accepted {args:?}");
        }
        let o=options(&["--stretch-tuning","440"]).unwrap();
        let missing=vec![geometry::demonstration_scale().unwrap()[0]];
        assert!(selected_scale(Some(&geometry::write_scale(&missing)),&o).unwrap_err()
            .contains("requires supplied octave neighbour"));
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
