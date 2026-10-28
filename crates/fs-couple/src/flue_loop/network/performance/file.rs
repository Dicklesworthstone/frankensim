//! File-driven jet, regional gas, tube graph and physical outlet pressure.
//!
//! The versioned source selects existing owners, including fractional travel,
//! optional gas-derived wide-tube loss, material wall motion, terminal reaction,
//! and the multirate Rayleigh receiver. Every physical selection is source data.
//! Local material walls, compact baffled loads and supplied mutual-radiation
//! operators retain their original limits; this is not a measured Boehm flute,
//! full shell model, resolved turbulent jet or broadband unbaffled qualification.
//!
//! Optional records after the duct and material walls, before boundaries or
//! radiation: `opening_losses COUNT` followed by COUNT records
//! `opening_loss NODE AREA_RATIO MAXIMUM_MACH SOURCE_TOKEN` (at most 64).
//! Each named node must carry a local RLC/relaxation terminal, either authored
//! or installed by `radiation compact-baffled`. Area comes from that terminal's
//! section and density/sound speed from its local gas, including mixtures.
//! The symmetric quasi-steady loss shares the current flow and original load
//! history; it does not replace radiation. Memoryless reflections, driven
//! moving apertures and mutual-radiation ports refuse this local closure.
//! AREA_RATIO is in (0,1], MAXIMUM_MACH in (0,0.3]; SOURCE_TOKEN is mandatory.

use super::{FlueMouthRadiation, FlueOutletConfig, FlueOutletPerformance, FlueRadiationInput};
use crate::bernoulli_aperture::fractional::FractionalTubePolicy;
use crate::bernoulli_aperture::network::{NetworkNode, TubeNetworkSpec, TubeSection};
use crate::bernoulli_aperture::viscothermal::ViscothermalSectionSpec;
use crate::bernoulli_aperture::radiation::environment::BaffledOutletBinding;
use crate::bernoulli_aperture::radiation::environment::multiport::input::{
    MimoRadiationInput, MimoRadiationInputError, MimoRadiationRecord,
};
use crate::bernoulli_aperture::radiation::environment::multiport::source::{
    ComputedBaffledInput, ComputedBaffledRecord, MimoSourceInput, MimoSourceRecord,
};
use crate::flue_loop::window_radiation::{WindowFlux, WindowRadiationLoad, cut_up_piston};
use crate::flue_loop::{FlueError, FlueGeometry, FlueVoice};
use crate::flue_loop::network::{FlueNetworkSpec, FlueSourcePaths};
use crate::pcm_wav::baffled::CircularOutletReceiver;
use crate::pcm_wav::baffled::openings::duct::{
    DuctOpening, DuctRadiationSpec, DuctReconstruction,
};
use crate::render::schedule::GestureCompileError;
use fs_aeroac::jetcard::{JetCard, mint_tonal_interim_card};
use fs_blake3::{ContentHash, hash_domain};
use fs_exec::CancelGate;
use fs_matdb::NormalizedMaterialCardPack;
use fs_material::gas::GasState;
use fs_material::gas_species::input::{GasInputError, GasRecord, ThermochemicalSources, read_gas_record};
use fs_scenario::gesture::{GestureSchedule, GestureTarget, GestureValue};
use std::str::{FromStr, Lines, SplitAsciiWhitespace};

mod duct;
mod wall;
mod material;
mod boundary;
mod window;
pub(crate) mod excitation;
use excitation::{FlueFlowInput, FlueJetInput};
pub use wall::{FlueWallRecord, MAX_FLUE_WALL_PACK_BYTES};
pub use window::FlueWindowRecord;
pub use boundary::{FlueBoundaryRecord, boundary_closing_force_path, flue_source_control_paths};

/// Versioned SI source grammar; the suffix is the existing gesture format.
pub const FLUE_PERFORMANCE_SCHEMA: &str = "frankensim-flue-performance-v1";
/// Complete source including an optional inline card and gesture suffix.
pub const MAX_FLUE_PERFORMANCE_BYTES: usize = 1024 * 1024;

/// Input failures retain their source line or original physics owner.
#[derive(Debug)]
pub enum FlueInputError {
    /// Malformed, unsupported or over-budget declaration.
    Input { /// One-based line, or zero for whole-file admission.
        line: usize, /// Refused condition.
        what: &'static str },
    /// Original material parsing, composition or phase refusal.
    Gas { /// One-based declaration line.
        line: usize, /// Shared material-owner error.
        source: GasInputError },
    /// The caller could not load or admit the exact declared material artifact.
    MaterialSource { /// One-based material declaration line.
        line: usize, /// Original source-loading message.
        source: String },
    /// Original material-state resolution, including missing claims or domains.
    Material { /// One-based material declaration line.
        line: usize, /// Original material owner's refusal.
        source: fs_material::state_point::MaterialStatePointError },
    /// Explicit mutual-radiation artifact loading or shared decoder refusal.
    RadiationSource { /// One-based radiation selection line.
        line: usize, /// Original source, gas, format or resource refusal.
        source: MimoRadiationInputError },
    /// Original jet-card schema or authority refusal.
    Card { /// One-based card line.
        line: usize, /// Shared source-card error.
        source: fs_aeroac::AeroacError },
    /// Geometry, jet, propagation or numerical admission failed.
    Physics(FlueError),
    /// Gesture, receiver or finite performance admission failed.
    Performance(GestureCompileError),
}
impl core::fmt::Display for FlueInputError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Input { line, what } => write!(f, "flue input line {line}: {what}"),
            Self::Gas { line, source } => write!(f, "flue gas line {line}: {source}"),
            Self::MaterialSource { line, source } => write!(f, "flue material source line {line}: {source}"),
            Self::Material { line, source } => write!(f, "flue material line {line}: {source}"),
            Self::RadiationSource { line, source } => write!(f, "flue radiation source line {line}: {source}"),
            Self::Card { line, source } => write!(f, "flue card line {line}: {source}"),
            Self::Physics(source) => write!(f, "flue physics: {source}"),
            Self::Performance(source) => write!(f, "flue performance: {source}"),
        }
    }
}
impl std::error::Error for FlueInputError {}
fn bad(line: usize, what: &'static str) -> FlueInputError { FlueInputError::Input { line, what } }

/// The held note's pressure: the highest the gesture's blowing tracks reach, the legacy target or
/// the declared source's pressure path. It sets the Reynolds number of the forced jet's integration.
pub(crate) fn held_blowing_pressure(schedule: &GestureSchedule, source: Option<&FlueSourcePaths>)
    -> f64
{
    schedule.tracks().iter()
        .filter(|track| match &track.target {
            GestureTarget::BlowingPressure => true,
            GestureTarget::Parameter(path) => source.is_some_and(|paths| *path == paths.pressure),
            _ => false,
        })
        .flat_map(|track| core::iter::once(&track.initial)
            .chain(track.events.iter().map(|event| &event.value)))
        .filter_map(|value| match value { GestureValue::PressurePa(p) => Some(*p), _ => None })
        .fold(0.0f64, f64::max)
}
fn poll(gate: &CancelGate) -> Result<(), FlueInputError> {
    if gate.is_requested() { Err(bad(0, "flue source construction cancelled")) } else { Ok(()) }
}

/// Exact source identity and clocks; none is a measured-instrument certificate.
#[derive(Clone, Copy, Debug)]
pub struct FluePerformanceInfo {
    /// Hash of all original source bytes, including inline card and gestures.
    pub input_hash: ContentHash,
    /// Jet and waveguide sample clock [Hz].
    pub mechanics_rate_hz: u32,
    /// Physical microphone output clock [Hz].
    pub sample_rate_hz: u32,
    /// Finite OUTPUT samples, including the caller's declared tail window.
    pub samples: u64,
    /// Positive pascals mapped to PCM full scale, never peak normalization.
    pub full_scale_pa: f64,
    /// Existing card owner's content identity.
    pub jet_card_hash: u64,
    /// Exact explicit selector, or an inline canonical card.
    pub jet_card_source: &'static str,
    /// Pressure, optional jet-angle and physical force assignments on one clock.
    pub compiled_controls: usize,
}

/// A complete source-to-microphone runtime and its original material records.
pub struct FlueFilePerformance {
    info: FluePerformanceInfo,
    renderer: FlueOutletPerformance,
    geometry: FlueGeometry,
    card: JetCard,
    gas_records: Vec<GasRecord>,
    wall_records: Vec<FlueWallRecord>,
    boundary_records: Vec<FlueBoundaryRecord>,
    radiation_bindings: Vec<BaffledOutletBinding>,
    mimo_radiation_record: Option<MimoRadiationRecord>,
    computed_radiation_record: Option<ComputedBaffledRecord>,
    window: Option<FlueWindowRecord>,
}
impl FlueFilePerformance {
    /// Parse the fixed-order source shown in `examples/flue.performance`.
    /// `jet_card tonal-interim` explicitly selects the retained stage-I card;
    /// `jet_card json {CANONICAL_JSON}` supplies a complete one-line card.
    /// An optional `windway none` or `windway slit LENGTH_M` row after `flue`
    /// selects the source-gas viscous windway. `windway slit-developing LENGTH_M` instead holds
    /// the jet by Shah's laminar developing-flow law (`bernoulli_aperture::developing`), refused
    /// unless the held note's pressure lies inside the law's validity. `windway unsteady LENGTH_M
    /// TRANSVERSE_CELLS CONVECTION_CELLS` instead resolves the evolving channel
    /// profile and jet convection, using explicit refinement counts. Its length is actual
    /// channel geometry, never an inferred lip channel for a transverse flute.
    /// Optional `mouth_loss contraction AREA_RATIO MAXIMUM_MACH SOURCE_TOKEN`
    /// after the windway row binds the existing separated-flow law to the actual
    /// mouth window and inlet gas. Omission or `mouth_loss none` retains the
    /// original linear mouth. The closure's source and Mach limit stay explicit.
    /// Optional `jet_speed card-centreline` after the `jet` row feeds the jet card the exit's centre
    /// speed (`FlueVoice::with_card_centreline_speed`, pexfp.135.5) where by default it reads the
    /// windway's bulk speed; the bulk speed still carries the flux. Omission keeps every default
    /// result bit-identical.
    /// Optional `jet_transfer nolle-price-2015 long-channel|short-channel CHAMFER_DEPTH_M` after the
    /// `jet` row replaces the island's receptivity chain by the published Nolle/Price jet transfer
    /// (`flue_loop::nolle_jet`), with the measured phase velocity of the channel class. Omission keeps
    /// the island's own chain and every default result bit-identical.
    /// Optional `jet_response forced-bickley SPACING WINDOW` in the same place instead replaces that
    /// chain by the forced Bickley jet's causal response (`flue_loop::forced_jet`). It is integrated
    /// at grid SPACING [b] over WINDOW [b/U0] at the Reynolds number of the gesture's highest
    /// blowing pressure; register physics is graded on 0.1 90. The two rows exclude each other.
    /// The `jet` row's DISPLACEMENT_GAIN is the default chain's authored gain. The cards' 90 is
    /// about 7-18x the forced jet's response at each note's speaking mode (mail 7491). It stands
    /// until the forced jet passes fg_003, fg_004, fj_001 and gas_direction_amount (ruling 7518,
    /// mail 7859), and is
    /// unused while either jet row is present.
    /// Gains, seed, dimensions, limits, clocks and the pressure phrase are data.
    /// The gesture suffix may also contain one typed `jet-angle` track in radians:
    /// it changes the existing geometric jet split at fixed cut-up and mouth area.
    /// `exterior same` shares the default gas; another gas record selects the
    /// same pressure-matched exterior for observation and any chosen reaction.
    /// `network` replaces
    /// `tube` with explicit branches, section gases and passive loads. The jet
    /// always uses the actual inlet section. A network explicitly names its
    /// radiating outlets; each outlet contributes its own accepted volume flow
    /// to the same physical receiver.
    /// Optional `mouth_radiation X_M Y_M WIDTH_CELLS CUT_UP_CELLS` after the
    /// outlet records also observes the solved outward source-mouth flow. Its
    /// rectangular size defaults to flue width by cut-up. Preceding optional
    /// `mouth_aperture rectangle WIDTH_M LENGTH_M` declares an independent
    /// acoustic opening without changing jet geometry. All apertures share
    /// the declared exterior and rigid baffle. No extra source load is inferred.
    /// Immediately following it, optional `mouth_reaction compact-baffled`
    /// replaces the old mouth end correction and linear loss with the same
    /// rectangle/exterior's checked compact self reaction. The last `flue`
    /// field (`mouth_loss_rel`) must then be zero. An explicit
    /// `mouth_reaction none` or omission preserves the original source load.
    /// The declared receiver band must satisfy `k*mouth_diagonal <= 0.5`,
    /// `dt*fmax <= 0.1` and the fixed 5% sampled response checks. Reaction and
    /// observation share the supplied mouth cells (at most 1024 in total).
    /// After the mouth records, optional `window wall-cut AXIAL_CELLS` replaces
    /// the legacy `0.61 a_eq` mouth inertance with the window's mass computed
    /// from the card's own channel, pocket and bore (see [`FlueWindowRecord`]),
    /// and keeps the channel in the graph only as a one-sample inlet line;
    /// `window none` or omission keeps it. It refuses a compact mouth reaction.
    /// `mouth_reaction window-cylinder` instead radiates that window's solved
    /// outer face from the head's outer cylinder in the declared exterior
    /// (see [`crate::flue_loop::window_radiation`]), replacing the authored
    /// linear loss: it needs `window wall-cut` and `mouth_loss_rel` zero.
    /// `mouth_reaction window-cylinder-piston` radiates one uniform piston over
    /// that face's angle and the window's cut-up instead, as a bracket.
    ///
    /// After the duct records, optional `walls COUNT FMIN_HZ FMAX_HZ` introduces
    /// COUNT rows `wall NODE isotropic PACK T_K OMEGA_RAD_S THICKNESS_M AXIAL LOSS`.
    /// AXIAL is `stress-free` or `strain-fixed`; LOSS is `elastic`,
    /// `kelvin-voigt`, or `maxwell COUNT MAX_DT_OVER_TAU MAX_ANGULAR_STEP
    /// MODULUS_PROPERTY TIME_PROPERTY ...`. Maxwell uses 1..=8 source-property
    /// pairs in Pa and seconds, with `young_modulus` as the equilibrium modulus.
    /// Each coefficient must cover the wall band; the actual mechanics clock
    /// must resolve the declared time and instantaneous-frequency limits.
    /// Alternatively `wall NODE orthotropic PACK T_K OMEGA_RAD_S THICKNESS_M
    /// AXIAL HOOP_AXIS AXIAL_AXIS LOSS` uses distinct one-based principal axes
    /// following the local cylindrical directions. Its LOSS is `elastic`,
    /// `proportional-kelvin-voigt TIME_PROPERTY`, or `normal-kelvin-voigt
    /// MAX_DT_OVER_TAU MAX_ANGULAR_STEP HH HA AA`. The last three properties
    /// supply the symmetric 2-D normal viscosity in Pa s. Free axial stress
    /// retains its exact Maxwell memory and residual Kelvin resistance; fixed
    /// axial strain selects the original hoop entry. These are the same
    /// source-backed reductions used by DG walls. No isotropic replacement is made.
    /// The same source selector and passive wall evolution serve reed ducts.
    /// Each source binds an existing empty degree-two cylindrical
    /// junction; geometry and gas assignments stay unchanged. The original wall
    /// owner covers half of each adjacent physical section, leaving unselected
    /// and endpoint half-strips rigid. `walls none` or no row retains rigid walls.
    ///
    /// Optional `boundaries COUNT` follows the wall records and binds every
    /// `duct_node driven` in ascending node order to a material cantilever slit.
    /// Each boundary has the following records (all numeric quantities are SI):
    /// ```text
    /// boundary NODE ASSEMBLY_NAME PART_NAME INITIAL_CLOSING_FORCE_N EXTERIOR_PRESSURE_PA
    /// boundary_plate PACK T_K OMEGA_RAD_S LENGTH_M WIDTH_M THICKNESS_M NX NY REST_GAP_M INITIAL_GAP_M INITIAL_VELOCITY_M_S
    /// boundary_mode FMIN_HZ FMAX_HZ MODE_INDEX MAX_SLIT_VARIATION MAX_SLOPE
    /// boundary_contact STIFFNESS_N_PER_M_ALPHA ALPHA INTERNAL_LOSS_S_PER_M
    /// boundary_slit LENGTH_M
    /// boundary_loss maxwell COUNT FMIN_HZ FMAX_HZ MAX_DT_OVER_TAU MAX_ANGULAR_STEP
    /// boundary_term MODULUS_PROPERTY TIME_PROPERTY
    /// ```
    /// `boundary_loss elastic` instead selects an explicit lossless solid.
    /// `boundary_loss kelvin-voigt` selects `kelvin_voigt_bending_viscosity`
    /// [Pa s] from the same card/query and projects it on the retained plate
    /// mode. This selection is exclusive with Maxwell memory, and even zero
    /// viscosity retains its source applicability checks. The
    /// `boundary_slit` row is optional; `none` selects the original separated
    /// Bernoulli flow. A supplied streamwise length adds the existing slit law
    /// using the actual chimney gas's viscosity and density. NX/NY are rectangle
    /// cell counts, with x=0 clamped and x=L the moving slit. Mode selection uses
    /// the original DKT pencil; no mass, stiffness, pressure area, frequency or
    /// damping coefficient is fitted from the instrument name. `young_modulus`
    /// is the equilibrium modulus when Maxwell terms are selected; every term's
    /// modulus [Pa] and relaxation time [s] is resolved from that same card/state.
    /// The finite exact structural block carries the original material memories.
    /// Every selected solid property must be constant at the declared temperature
    /// and applicable over the Maxwell use band, or the mode search band for an
    /// elastic or Kelvin-Voigt solid. A sampled frequency-dependent modulus is not silently frozen.
    ///
    /// [`boundary_closing_force_path`] derives the exact part identity and force
    /// dimensions from ASSEMBLY_NAME/PART_NAME. A canonical v2 `Parameter` track
    /// with that path and `force-n` values drives nonnegative closing force;
    /// ramps use the shared gesture owner. Omitting that track retains the declared
    /// initial force. Before `controls`, `control_source ASSEMBLY_NAME PART_NAME` binds
    /// the jet source through [`flue_source_control_paths`]. Use its pressure and
    /// optional angle paths in the same v2 schedule; v2 requires every track to
    /// be a resolved path and never mixes legacy target labels with force paths.
    /// Controls preserve all chimney waves, trapped pressure, contact and
    /// Maxwell history. The source is a uniform moving slit with one plate mode;
    /// it does not infer a circular pad, key linkage or moving three-dimensional
    /// opening, and its exterior pressure is prescribed rather than a second
    /// nonlinear mutual-radiation solve. Other passive outlets retain their loads.
    ///
    /// An optional `radiation prescribed` or `radiation compact-baffled` row then
    /// selects whether every observed outlet retains its supplied load or receives
    /// the existing compact passive load in the same exterior as the receiver.
    /// `radiation mimo ARTIFACT MAX_STATES MAX_MEMORY_BYTES PASSIVITY_POINTS`
    /// instead loads one explicit reciprocal rational operator for all declared
    /// outlets, in their original order. The shared artifact decoder retains its
    /// independently declared gas, geometry, identities and band; these must
    /// agree with the actual receiver and cover the mechanical Nyquist frequency.
    /// Preexisting reactive loads refuse replacement. Mutual feedback uses the
    /// existing terminal owner and does not fit or rescale the supplied operator.
    ///
    /// The original jet card keeps its own narrow lab scope. The complete
    /// authored composition remains EstimateOnly, including outside-card
    /// diagnostics; loading never promotes a flute accuracy or SPL claim.
    ///
    /// # Errors
    /// Malformed/over-budget input, missing card, original material, geometry,
    /// jet, fractional/loss/receiver domain, gesture or cancellation refusal.
    pub fn from_bytes(bytes: &[u8], max_block: usize, gate: &CancelGate)
        -> Result<Self, FlueInputError>
    {
        Self::from_bytes_with_gas_sources(bytes, max_block, gate, None)
    }

    /// The same path with explicitly loaded NASA-9/transport source packs.
    /// All gas forms use `fs_material::gas_species::input::read_gas_record`.
    pub fn from_bytes_with_gas_sources(bytes: &[u8], max_block: usize, gate: &CancelGate,
        sources: Option<&ThermochemicalSources<'_>>)
        -> Result<Self, FlueInputError>
    {
        Self::from_bytes_with_sources(bytes, max_block, gate, sources,
            |_| Err("wall material-card packs require an explicit source loader".into()))
    }

    /// Resolve supplied normalized material cards through the existing material
    /// state-point owner before constructing the jet, walls, bore and receiver.
    /// The caller loads declared packs; no filesystem I/O enters this library.
    /// Sources use explicit bare-SI T and omega queries and SingleClaimOnly
    /// selection. Missing, ambiguous, dimensionally incompatible or out-of-band
    /// properties refuse without a material-name fallback. Shared pack paths load
    /// once; each wall retains its independent solid temperature and query receipt.
    ///
    /// Original geometry, regional gases, loads and outlet addresses are preserved.
    /// Optional compact reaction uses the SAME exterior and band as the receiver.
    /// The complete file remains an estimated composition of the stated images.
    pub fn from_bytes_with_sources(bytes: &[u8], max_block: usize, gate: &CancelGate,
        sources: Option<&ThermochemicalSources<'_>>,
        load_material: impl FnMut(&str) -> Result<NormalizedMaterialCardPack, String>)
        -> Result<Self, FlueInputError>
    {
        Self::from_bytes_with_all_sources(bytes, max_block, gate, sources, load_material,
            |_| Err("MIMO radiation artifacts require an explicit source loader".into()))
    }

    /// The same native performance with imported or computed mutual radiation.
    /// Both loaders resolve paths relative to this performance's source directory;
    /// the library performs no filesystem I/O. The common MIMO decoder resolves
    /// its own gas record using the same optional thermochemical registry. A
    /// changed atmosphere requires a matching supplied operator for `mimo`;
    /// stale artifacts refuse instead of silently replacing their baked gas or
    /// coefficients. `computed-baffled` instead samples the actual source
    /// geometry and exterior through the shared physical bake owner, without
    /// calling the artifact loader.
    /// Walls, moving boundaries, mouth laws, controls and receiver clocks remain
    /// in their original owners. Omitting the MIMO selector never calls its loader.
    pub fn from_bytes_with_all_sources(bytes: &[u8], max_block: usize, gate: &CancelGate,
        sources: Option<&ThermochemicalSources<'_>>,
        mut load_material: impl FnMut(&str) -> Result<NormalizedMaterialCardPack, String>,
        mut load_radiation: impl FnMut(&str) -> Result<Vec<u8>, String>)
        -> Result<Self, FlueInputError>
    {
        poll(gate)?;
        if bytes.len() > MAX_FLUE_PERFORMANCE_BYTES || !(1..=65_536).contains(&max_block) {
            return Err(bad(0, "flue source exceeds its byte or callback budget"));
        }
        let text = std::str::from_utf8(bytes).map_err(|_| bad(0, "flue source must be UTF-8"))?;
        let (configuration, gestures) = text.split_once("\ngesture\n")
            .ok_or_else(|| bad(0, "missing LF-delimited gesture section"))?;
        let (header, body) = configuration.split_once('\n').ok_or_else(|| bad(1, "missing flue schema header"))?;
        if header != FLUE_PERFORMANCE_SCHEMA { return Err(bad(1, "unsupported flue source schema")); }
        let mut rows = Rows { lines: body.lines(), line: 1 };
        let mut row = rows.row("clock")?;
        let mechanics_rate_hz: u32 = row.number()?;
        let sample_rate_hz: u32 = row.number()?;
        row.end()?;
        if !(8_000..=192_000).contains(&sample_rate_hz) || mechanics_rate_hz == 0
            || mechanics_rate_hz % sample_rate_hz != 0
            || !(1..=16).contains(&(mechanics_rate_hz/sample_rate_hz)) {
            return Err(bad(rows.line, "output clock must be 8..192 kHz and mechanics an integer multiple in 1..=16"));
        }
        let mut row = rows.row("output")?;
        let samples: u64 = row.number()?;
        let full_scale_pa = row.scalar()?;
        let maximum_abs_pressure_pa = row.scalar()?;
        row.end()?;
        if full_scale_pa <= 0.0 { return Err(bad(rows.line, "full scale must be positive pascals")); }
        let mut gas_records = Vec::new();
        let mut row = rows.row("gas")?;
        let gas = row.gas(sources, &mut gas_records)?;
        row.end()?;
        let exterior_record = rows.record("exterior")?;
        let exterior = if exterior_record == "same" {
            gas
        } else {
            let mut row = Row { fields: exterior_record.split_ascii_whitespace(), line: rows.line };
            let state = row.gas(sources, &mut gas_records)?;
            row.end()?;
            state
        };
        let exterior_record = gas_records.last().expect("both gas routes retain their source record").clone();
        let mut duct = duct::Input::read(&mut rows, gas, sources, &mut gas_records, gate)?;
        let walls = wall::ParsedWalls::read(&mut rows)?;
        // Independent terminal separation loss, in series with whatever
        // physical radiation load the description selects below.
        let mut opening_losses = Vec::new();
        if rows.next_is("opening_losses") {
            let mut row = rows.row("opening_losses")?;
            let count: usize = row.number()?;
            row.end()?;
            if count > 64 { return Err(bad(rows.line, "opening_losses exceeds 64 endpoints")); }
            for _ in 0..count {
                let row = rows.row("opening_loss")?;
                opening_losses.push(crate::bernoulli_aperture::opening_loss::parse_record(row.fields)
                    .map_err(|what| bad(row.line, what))?);
            }
        }
        let boundaries = boundary::ParsedBoundaries::read(&mut rows)?;
        let (radiation_reaction, mimo_input) = if rows.next_is("radiation") {
            let record = rows.record("radiation")?;
            match record {
                "prescribed" => (false, None),
                "compact-baffled" => (true, None),
                _ => {
                    let Some((kind, arguments)) = record.split_once(char::is_whitespace) else {
                        return Err(bad(rows.line, "radiation requires prescribed, compact-baffled, mimo ARTIFACT LIMITS, or computed-baffled NUMERICAL_CONTROLS"));
                    };
                    let input = match kind {
                        "mimo" => MimoRadiationInput::parse(arguments).map(MimoSourceInput::Artifact),
                        "computed-baffled" => ComputedBaffledInput::parse(arguments).map(MimoSourceInput::ComputedBaffled),
                        _ => return Err(bad(rows.line, "unknown physical radiation source selection")),
                    }
                        .map_err(|source| FlueInputError::RadiationSource { line: rows.line, source })?;
                    (false, Some((rows.line, input)))
                }
            }
        } else { (false, None) };
        let inlet = duct.graph.inlet_section_index().map_err(|source|
            FlueInputError::Physics(FlueError::Physics(source)))?;
        let radius_m = duct.graph.sections[inlet].radius_m;
        let mut row = rows.row("flue")?;
        let geometry = FlueGeometry {
            flue_width_m: row.scalar()?, flue_height_m: row.scalar()?, cut_up_m: row.scalar()?,
            labium_offset_m: row.scalar()?, jet_angle_rad: row.scalar()?, mouth_loss_rel: row.scalar()?,
            bore_area_m2: core::f64::consts::PI * radius_m * radius_m,
        };
        row.end()?;
        let flow_input = FlueFlowInput::read(&mut rows)?;
        let card_source = rows.record("jet_card")?;
        let (card, jet_card_source) = if card_source == "tonal-interim" {
            (mint_tonal_interim_card(), "tonal-interim")
        } else if let Some(json) = card_source.strip_prefix("json ") {
            (JetCard::from_json(json), "inline-json")
        } else { return Err(bad(rows.line, "jet_card requires tonal-interim or json followed by a complete canonical card")); };
        let card = card.map_err(|source| FlueInputError::Card { line: rows.line, source })?;
        let mut row = rows.row("jet")?;
        let displacement_gain = row.scalar()?;
        let seed_amplitude_rel = row.scalar()?;
        let voice_seed = row.number()?;
        let voice_index = row.number()?;
        row.end()?;
        let jet_input = FlueJetInput::read(&mut rows)?;
        let mut row = rows.row("jet_limits")?;
        let minimum_jet_speed_m_s = row.scalar()?;
        let maximum_jet_mach = row.scalar()?;
        let max_jet_history_samples = row.number()?;
        row.end()?;
        let mut row = rows.row("wave_memory_bytes")?;
        let max_wave_memory_bytes: usize = row.number()?;
        row.end()?;
        if max_wave_memory_bytes > 256*1024*1024 { return Err(bad(rows.line, "flue wave memory exceeds 256 MiB")); }
        let mut row = rows.row("fractional")?;
        let policy = FractionalTubePolicy { maximum_frequency_hz: row.scalar()?,
            maximum_summed_phase_error_rad: row.scalar()? };
        row.end()?;
        let viscothermal = duct.read_losses(&mut rows)?;
        let mut row = rows.row("receiver")?;
        let receiver = CircularOutletReceiver { position_m: [row.scalar()?,row.scalar()?,row.scalar()?],
            radial_rings: row.number()?, angular_points: row.number()?, maximum_frequency_hz: row.scalar()? };
        row.end()?;
        if receiver.maximum_frequency_hz > policy.maximum_frequency_hz {
            return Err(bad(rows.line, "fractional phase band must cover the receiver's declared band"));
        }
        let (outlet_node, additional_outlets) = duct.read_outlets(&mut rows)?;
        let aperture_size_m = read_mouth_aperture(&mut rows)?;
        let mouth_radiation = if rows.next_is("mouth_radiation") {
            let mut row = rows.row("mouth_radiation")?;
            let input = FlueMouthRadiation { center_m: [row.scalar()?, row.scalar()?], aperture_size_m,
                cells: [row.number()?, row.number()?] };
            row.end()?;
            Some(input)
        } else { None };
        if aperture_size_m.is_some() && mouth_radiation.is_none() {
            return Err(bad(rows.line, "mouth_aperture requires mouth_radiation"));
        }
        let (mouth_reaction, window_reaction) = if rows.next_is("mouth_reaction") {
            let mut row = rows.row("mouth_reaction")?;
            let selected = match row.fields.next() {
                Some("none") => (false, None),
                Some("compact-baffled") => (true, None),
                Some("window-cylinder") => (false, Some(WindowFlux::WallWindow)),
                Some("window-cylinder-piston") => (false, Some(WindowFlux::CutUpPiston)),
                _ => return Err(bad(rows.line,
                    "mouth_reaction requires none, compact-baffled, window-cylinder or window-cylinder-piston")),
            };
            row.end()?;
            if selected.0 && mouth_radiation.is_none() {
                return Err(bad(rows.line, "mouth reaction requires an explicit mouth_radiation rectangle"));
            }
            selected
        } else { (false, None) };
        let window = match window::read(&mut rows)? {
            Some(_) if mouth_reaction => return Err(bad(rows.line,
                "window wall-cut replaces the legacy mouth inertance; a compact mouth reaction brings its own outer mass")),
            Some(cells) => Some(window::solve(&mut duct.graph, &duct.gases, mechanics_rate_hz,
                &geometry, cells, rows.line)?),
            None => None,
        };
        if window_reaction.is_some() && window.is_none() {
            return Err(bad(rows.line,
                "mouth_reaction window-cylinder radiates the window's outer face; declare window wall-cut"));
        }
        let control_source = boundary::read_control_source(&mut rows)?;
        let mut row = rows.row("controls")?;
        let max_controls: usize = row.number()?;
        let max_compile_work = row.number()?;
        row.end()?;
        if max_controls > 262_144 { return Err(bad(rows.line, "flue control budget exceeds 262144 ticks")); }
        if rows.lines.next().is_some() { return Err(bad(rows.line+1,"unexpected records before gesture section")); }
        for line in gestures.lines() {
            if let Some(count) = line.strip_prefix("events\t") {
                if count.parse::<usize>().ok().is_none_or(|n| n>16384) {
                    return Err(bad(rows.line+1,"gesture event count exceeds 16384-event input budget"));
                }
            }
        }
        let schedule = GestureSchedule::from_canonical_bytes(gestures.as_bytes())
            .map_err(|_| bad(rows.line+1, "invalid canonical pressure/jet-angle/force gesture suffix"))?;
        if schedule.to_canonical_bytes() != gestures.as_bytes() {
            return Err(bad(rows.line+1,"gesture suffix must be canonical without ignored trailing data"));
        }
        poll(gate)?;
        duct.graph.max_wave_memory_bytes = max_wave_memory_bytes;
        let source = FlueNetworkSpec { graph: duct.graph, gases: duct.gases,
            viscothermal, minimum_jet_speed_m_s, maximum_jet_mach, max_jet_history_samples };
        let mut material_sources = material::Sources::default();
        let (walls, wall_records) = match walls {
            Some(walls) => {
                let (walls, records) = walls.resolve(policy.maximum_frequency_hz, &mut material_sources, &mut load_material, gate)?;
                (Some(walls), records)
            }
            None => (None, Vec::new()),
        };
        poll(gate)?;
        let radiation = || {
            let mut openings = vec![DuctOpening { network: 0, node: outlet_node, center_m: [0.0; 2] }];
            openings.extend(additional_outlets.iter().map(|outlet| DuctOpening {
                network: 0, node: outlet.node, center_m: outlet.center_m,
            }));
            DuctRadiationSpec { openings, exterior: Some(exterior), receiver,
                reconstruction: DuctReconstruction::Differentiating {
                    mechanics_rate_hz, output_rate_hz: sample_rate_hz,
                },
            }
        };
        let (voice, radiation_bindings) = if radiation_reaction {
            let (voice, environment) = FlueVoice::from_baffled_environment(Some(&card),
                displacement_gain, seed_amplitude_rel, geometry, source, walls, Some(policy),
                radiation(), voice_seed, voice_index, gate).map_err(FlueInputError::Physics)?;
            for binding in environment.bindings() {
                binding.validate_receiver(voice.network().expect("constructed graph").graph(),
                    environment.radiation(), exterior)
                    .map_err(|error| FlueInputError::Physics(FlueError::Physics(
                        crate::acoustic_realize::AcousticRealizeError::Nonlinear(error.to_string()))))?;
            }
            (voice, environment.into_parts().1)
        } else {
            (FlueVoice::from_fractional_network(Some(&card), displacement_gain, seed_amplitude_rel,
                geometry, source, walls, policy, mechanics_rate_hz, voice_seed, voice_index, gate)
                .map_err(FlueInputError::Physics)?, Vec::new())
        };
        let held_pressure_pa = held_blowing_pressure(&schedule, control_source.as_ref());
        let voice = flow_input.attach(voice, held_pressure_pa)?;
        let voice = jet_input.attach(voice, held_pressure_pa, gate)?;
        let voice = match &window {
            Some(record) => voice.with_window_mass(record.ladder.mass_per_density,
                record.riser_mass_per_density).map_err(FlueInputError::Physics)?,
            None => voice,
        };
        let voice = match (window_reaction, &window) {
            (Some(flux), Some(record)) => {
                // The finest rung's outer face carries the flux distribution.
                let face = &record.ladder.rungs[record.ladder.rungs.len() - 1].outer_cells;
                let piston = [cut_up_piston(face, record.window.axial_length)];
                let cells = match flux {
                    WindowFlux::WallWindow => face.as_slice(),
                    WindowFlux::CutUpPiston => &piston[..],
                };
                let load = WindowRadiationLoad::new(
                    record.window.bore_radius + record.window.wall_thickness, cells, flux,
                    exterior, voice.sample_period_s(), gate).map_err(FlueInputError::Physics)?;
                voice.with_window_radiation(load).map_err(FlueInputError::Physics)?
            }
            _ => voice,
        };
        let mechanical_samples = samples.checked_mul(u64::from(mechanics_rate_hz / sample_rate_hz))
            .ok_or_else(|| bad(0, "moving boundary step window overflows"))?;
        let (voice, boundary_records) = boundaries.bind(voice, mechanical_samples,
            &mut material_sources, &mut load_material, gate)?;
        let voice = match control_source {
            Some(paths) => voice.with_source_control_paths(paths).map_err(FlueInputError::Physics)?,
            None => voice,
        };
        let (voice, mimo_radiation_record, computed_radiation_record) = match mimo_input {
            Some((line, input)) => {
                let context = radiation();
                let original = voice.network().expect("constructed physical flue graph").source();
                let (model, record) = input.build(&original.graph, &original.gases, &context,
                    &exterior_record, sources, &mut load_radiation, gate)
                    .map_err(|source| FlueInputError::RadiationSource { line, source })?;
                let voice = voice.with_mimo_environment(model, context, input.limits(), gate)
                    .map_err(FlueInputError::Physics)?;
                match record {
                    MimoSourceRecord::Artifact(record) => (voice, Some(record), None),
                    MimoSourceRecord::ComputedBaffled(record) => (voice, None, Some(record)),
                }
            }
            None => (voice, None, None),
        };
        let voice = if opening_losses.is_empty() { voice } else {
            voice.with_opening_losses(opening_losses, gate).map_err(FlueInputError::Physics)?
        };
        let mut renderer = FlueOutletPerformance::with_multirate_radiation(voice, schedule,
            FlueOutletConfig { sample_rate_hz, samples, max_block, max_controls, max_compile_work,
                outlet_node, receiver, maximum_abs_pressure_pa },
            FlueRadiationInput { exterior, additional_outlets }, mechanics_rate_hz)
            .map_err(FlueInputError::Performance)?;
        if let Some(mouth) = mouth_radiation {
            renderer = renderer.with_mouth_radiation(mouth).map_err(FlueInputError::Performance)?;
        }
        if mouth_reaction {
            renderer = renderer.with_mouth_reaction(gate).map_err(FlueInputError::Performance)?;
        }
        poll(gate)?;
        let info = FluePerformanceInfo { input_hash: hash_domain("org.frankensim.fs-couple.flue-performance.v1",bytes),
            mechanics_rate_hz, sample_rate_hz, samples, full_scale_pa,
            jet_card_hash: card.content_hash(), jet_card_source, compiled_controls: renderer.controls.len() };
        Ok(Self { info, renderer, geometry, card, gas_records, wall_records, boundary_records,
            radiation_bindings, mimo_radiation_record, computed_radiation_record, window })
    }
    /// Original identity and finite clocks, separate from the evolving solver.
    #[must_use] pub const fn info(&self) -> FluePerformanceInfo { self.info }
    /// Actual physical runtime; inspect its original and lowered gas/tube inputs.
    #[must_use] pub const fn renderer(&self) -> &FlueOutletPerformance { &self.renderer }
    /// Supplied mouth geometry, with bore area derived from the inlet section.
    #[must_use] pub const fn geometry(&self) -> FlueGeometry { self.geometry }
    /// The complete selected card, including its narrower lab-validity scope.
    #[must_use] pub const fn card(&self) -> &JetCard { &self.card }
    /// Gas records and thermochemical receipts, in source order.
    #[must_use] pub fn gas_records(&self) -> &[GasRecord] { &self.gas_records }
    /// Original normalized solid sources, in material-wall assignment order.
    #[must_use] pub fn wall_records(&self) -> &[FlueWallRecord] { &self.wall_records }
    /// Original material cards and typed force paths of moving chimney boundaries.
    #[must_use] pub fn boundary_records(&self) -> &[FlueBoundaryRecord] { &self.boundary_records }
    /// Exact installed reaction loads, in the observed outlet order.
    #[must_use] pub fn radiation_bindings(&self) -> &[BaffledOutletBinding] { &self.radiation_bindings }
    /// Exact artifact bytes' identity and declared source gas/geometry, separate
    /// from the shared terminal model's evolving state and passivity diagnostics.
    #[must_use] pub const fn mimo_radiation_record(&self) -> Option<&MimoRadiationRecord> {
        self.mimo_radiation_record.as_ref()
    }
    /// Actual geometry/gas and numerical results of computed mutual radiation.
    #[must_use] pub const fn computed_radiation_record(&self) -> Option<&ComputedBaffledRecord> {
        self.computed_radiation_record.as_ref()
    }
    /// The window derived from the card and its computed mass, when selected.
    #[must_use] pub const fn window(&self) -> Option<&FlueWindowRecord> { self.window.as_ref() }
    /// Move the original physical pressure producer into streaming/PCM output.
    #[must_use] pub fn into_renderer(self) -> FlueOutletPerformance { self.renderer }
}

/// Network acoustic aperture record; no change to the flue source.
pub(in crate::flue_loop) fn read_mouth_aperture(rows: &mut Rows<'_>)
    -> Result<Option<[f64; 2]>, FlueInputError>
{
    if !rows.next_is("mouth_aperture") { return Ok(None); }
    let mut row = rows.row("mouth_aperture")?;
    if row.fields.next() != Some("rectangle") {
        return Err(bad(row.line, "mouth_aperture requires rectangle WIDTH_M LENGTH_M"));
    }
    let size = [row.scalar()?, row.scalar()?];
    if size.iter().any(|x| *x <= 0.0) {
        return Err(bad(row.line, "mouth aperture dimensions must be positive finite"));
    }
    row.end()?;
    Ok(Some(size))
}

pub(in crate::flue_loop) struct Rows<'a> {
    pub(in crate::flue_loop) lines: Lines<'a>,
    pub(in crate::flue_loop) line: usize,
}
impl<'a> Rows<'a> {
    pub(in crate::flue_loop) fn next_is(&self, expected: &str) -> bool {
        self.lines.clone().next().is_some_and(|line|
            line.trim_ascii().split_ascii_whitespace().next() == Some(expected))
    }
    pub(in crate::flue_loop) fn record(&mut self, expected: &'static str) -> Result<&'a str, FlueInputError> {
        self.line += 1;
        let text = self.lines.next().ok_or_else(|| bad(self.line,"missing source record"))?.trim_ascii();
        let key = text.split_ascii_whitespace().next().unwrap_or("");
        if key != expected { return Err(bad(self.line, expected)); }
        Ok(text[key.len()..].trim_ascii())
    }
    pub(in crate::flue_loop) fn row(&mut self, expected: &'static str) -> Result<Row<'a>, FlueInputError> {
        let text = self.record(expected)?;
        Ok(Row { fields: text.split_ascii_whitespace(), line: self.line })
    }
}
pub(in crate::flue_loop) struct Row<'a> {
    pub(in crate::flue_loop) fields: SplitAsciiWhitespace<'a>,
    pub(in crate::flue_loop) line: usize,
}
impl Row<'_> {
    pub(in crate::flue_loop) fn number<T: FromStr>(&mut self) -> Result<T, FlueInputError> {
        self.fields.next().and_then(|s| s.parse().ok()).ok_or_else(|| bad(self.line,"missing or invalid numeric field"))
    }
    pub(in crate::flue_loop) fn scalar(&mut self) -> Result<f64, FlueInputError> {
        let x: f64 = self.number()?;
        if !x.is_finite() { return Err(bad(self.line,"physical scalars must be finite")); }
        Ok(x)
    }
    pub(in crate::flue_loop) fn end(mut self) -> Result<(), FlueInputError> {
        if self.fields.next().is_some() { Err(bad(self.line,"unexpected trailing fields")) } else { Ok(()) }
    }
    pub(in crate::flue_loop) fn gas(&mut self, sources: Option<&ThermochemicalSources<'_>>, records: &mut Vec<GasRecord>)
        -> Result<GasState, FlueInputError>
    {
        let record = read_gas_record(&mut self.fields, None, sources)
            .map_err(|source| FlueInputError::Gas { line: self.line, source })?;
        let gas = *record.state();
        records.push(record);
        Ok(gas)
    }
}
