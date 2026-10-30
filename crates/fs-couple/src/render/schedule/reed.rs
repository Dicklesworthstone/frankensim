//! Authored reed/duct performances over the existing physical and gesture owners.
//!
//! The file supplies primitive SI parameters, a static segmented duct, an explicit
//! dry/humid-air or named-mixture default, optional segment gases, and one canonical
//! blowing-pressure track. An optional gas_field block supplies primitive axial
//! temperatures and mole fractions through explicitly loaded source packs.
//! An alternative gas_steady block generates that same field from the existing
//! conservative thermal/species owner, using the authored bore geometry.
//! An explicit-species source can declare `molecular`, `molecular_resolution`
//! and `molecular_wall` records before `compile_limits`, with either its
//! homogeneous ambient or the same continuous authored/injected/steady field.
//! A molecular field requires the explicit continuous spatial-resolution axis.
//! The shared molecular reservoirs then reach both pressure feedback and the
//! physical terminal and side mouths, independently of the listener's exterior
//! absorption law. Molecular side mouths apply that law at their own ranges.
//! It does not infer material identity
//! or turn a pitch into a medium property. Pressure commands replace the fixture
//! attack envelope; all retained reed and bore state survives release and resume.
//! Optional `slit_length LENGTH_M` after the reed record selects the shared
//! quasisteady viscous slit from the exact inlet gas, including axial fields.
//! Legacy `listener_m` selects bore pressure plus the compact-jet proxy.
//! `baffled_outlet` selects Rayleigh observation of the admitted openings.
//! Neither is measured calibration; the underlying reed mechanics are shared.
//!
//! Both observers then apply one observer absorption law over their range,
//! chosen as in `acoustic_realize`: a legacy (T, p, RH) medium gets ISO
//! 9613-1, and a supplied mixture or thermochemical medium gets
//! Stokes-Kirchhoff with its omission receipted. An optional
//! `observer_law none|iso9613|stokes-kirchhoff|relaxation
//! iso9613-terrestrial-air` record after the observation declares the law
//! instead. `none` is a receipted lossless choice for analytic oracles.
//! The law acts on radiated pressure only. The proxy's bore pressure is an
//! internal observation that travels no distance, so it passes unfiltered.

use crate::acoustic_realize::AcousticRealizeError;
use crate::air_path::{
    GasAbsorptionModel, ObserverFilter, ObserverFilterReceipt, ObserverOmission, legacy_relative_humidity,
    observer_band_top_hz, omitted_observer_absorption,
};
use crate::driving_point::regional::outlet::{OutletError, OutletReceiverInfo};
use fs_blake3::{ContentHash, hash_domain};
use fs_duct::{Duct, HoleState, Segment, Termination};
use fs_material::gas::GasState;
use fs_material::gas_species::input::{
    GasInputError, GasRecord, ThermochemicalSources, read_gas_record,
};
use fs_scenario::{BeatingReed, gesture::GestureSchedule};
use std::str::{FromStr, Lines, SplitAsciiWhitespace};
mod field_input;
mod molecular_input;
mod outlet;
use super::{
    GestureCompileError, PressureGestureBinding, ScheduledControl, ScheduledRenderer,
    compile_pressure_gestures,
};
use crate::pcm_wav::observation::{ObservationError, ObservationPolicy, PressureRenderer};
use crate::render::{ControlDelta, ReedBoreVoice, ReedLay, RenderContext, RenderError, RenderVoice};
use crate::thin_plate::PlateBank;
pub use field_input::steady::SteadyGasPerformance;

/// Versioned header; the schedule retains the shared canonical wire format.
pub const REED_PERFORMANCE_SCHEMA: &str = "frankensim-reed-performance-v1";
/// Bound input bytes before reading or invoking count-driven decoders.
pub const MAX_REED_PERFORMANCE_BYTES: usize = 1024 * 1024;
/// Authored axial and tone-hole records. A measured wind lattice needs one
/// record per bore station pair plus one per chimney and per split bore, so the
/// cap must admit a complete published bore table without decimating it.
const MAX_SEGMENTS: usize = 256;
const MAX_EVENTS: usize = 16_384;
const MAX_CONTROLS: usize = 262_144;
const MAX_WORK: u64 = 16_777_216;

/// Syntax, model and scheduling refusals retain their owning diagnostics.
#[derive(Debug)]
pub enum ReedPerformanceError {
    /// An observation is diagnostic-only or has inconsistent physical media.
    Observation(ObservationError),
    /// Invalid source record, with a one-based header line (zero for the file).
    Input {
        /// Header line, or zero for whole-file/schedule admission.
        line: usize,
        /// Failed condition.
        what: &'static str,
    },
    /// Material admission retains its original cause and source line.
    Gas {
        /// One-based ambient, segment, field or exterior record line.
        line: usize,
        /// Original species, composition or state refusal.
        source: fs_material::MaterialError,
    },
    /// The source-bound stationary heat/species owner refused publication.
    SteadyGas {
        /// One-based gas_steady header line.
        line: usize,
        /// Original thermal, caloric, phase, balance, work or mean-flow refusal.
        source: fs_material::gas_field::steady::SteadyGasError,
    },
    /// The original terminal transfer or physical receiver refused.
    Outlet(OutletError),
    /// The existing reed, TMM or gas realization refused.
    Model(AcousticRealizeError),
    /// The existing pressure compiler refused.
    Schedule(GestureCompileError),
    /// The existing common renderer refused.
    Render(RenderError),
}
impl core::fmt::Display for ReedPerformanceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Observation(e) => write!(f, "reed observation: {e}"),
            Self::Input { line, what } => write!(f, "reed input line {line}: {what}"),
            Self::Gas { line, source } => write!(f, "reed gas line {line}: {source}"),
            Self::SteadyGas { line, source } => write!(f, "reed steady gas line {line}: {source}"),
            Self::Model(e) => write!(f, "reed model: {e}"),
            Self::Outlet(e) => write!(f, "reed outlet: {e}"),
            Self::Schedule(e) => write!(f, "reed input: {e}"),
            Self::Render(e) => write!(f, "reed render: {e}"),
        }
    }
}
impl std::error::Error for ReedPerformanceError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Observation(source) => Some(source),
            Self::Gas { source, .. } => Some(source),
            Self::SteadyGas { source, .. } => Some(source),
            Self::Model(source) => Some(source),
            Self::Outlet(source) => Some(source),
            _ => None,
        }
    }
}
fn bad(line: usize, what: &'static str) -> ReedPerformanceError {
    ReedPerformanceError::Input { line, what }
}

/// The observer law a performance file declares with `observer_law`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum DeclaredObserverLaw {
    Lossless,
    Iso9613,
    StokesKirchhoff,
    IsoAirRelaxation,
}

/// The observer law applied to a performance's emitted pressure.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReedObserverReceipt {
    /// The law, its causal realization and its certificates. The baffled
    /// receiver's propagation delays use `filter.delay_speed_m_s`, so its
    /// delay plus the filter's low-frequency delay is `r / c_0` once. The
    /// legacy proxy applies no propagation delay.
    pub filter: ObserverFilterReceipt,
    /// The law is applied at the receiver's distance from the outlet centre.
    /// This is the largest difference from it over the outlet's Rayleigh
    /// elements [m]. `None` for the legacy proxy.
    pub range_spread_m: Option<f64>,
    /// The legacy proxy applies no propagation delay, so its filter's
    /// low-frequency delay (dispersion plus `filter.latency_s`) is
    /// uncompensated. A consumer placing the proxy's radiated part at its
    /// range must delay it by this much [s]: `filter.receiver_delay_s`, exactly
    /// what the baffled receiver applies, so both total `r / c_0`. Negative
    /// when the filter's delay exceeds the propagation delay at this range.
    /// `None` for the baffled receiver, which applies it itself.
    pub consumer_delay_s: Option<f64>,
    /// The proxy's sample is bore pressure plus radiation (plates and jet at
    /// `listener_m`). The law and its latency act on the radiation only; the
    /// bore pressure has no range and passes unfiltered and undelayed, so no
    /// single delay of the summed stream places both parts. `true` for the
    /// proxy, `false` for the baffled receiver, whose output is all radiated.
    pub internal_part_unpropagated: bool,
}

/// The medium between instrument and listener: the baffled receiver's
/// exterior gas (by default the outlet's own gas), or the ambient for the
/// legacy proxy.
fn observer_medium(
    baffled: Option<([f64; 3], Option<GasState>)>,
    ambient: &GasState,
    outlet: &GasState,
) -> GasState {
    match baffled {
        None => *ambient,
        Some((_, Some(exterior))) => exterior,
        Some((_, None)) => *outlet,
    }
}

/// Nearest and farthest element ranges of a baffled circular outlet of
/// `radius_m`, seen from `position_m`.
fn element_ranges(position_m: [f64; 3], radius_m: f64) -> (f64, f64) {
    let rho = position_m[0].hypot(position_m[1]);
    let z = position_m[2];
    ((rho - radius_m).max(0.0).hypot(z), (rho + radius_m).hypot(z))
}

/// Largest `|element range - centre range|` over a baffled circular outlet of
/// `radius_m`, seen from `position_m`.
fn range_spread(position_m: [f64; 3], radius_m: f64) -> f64 {
    let centre = position_m[0].hypot(position_m[1]).hypot(position_m[2]);
    let (near, far) = element_ranges(position_m, radius_m);
    (far - centre).max(centre - near)
}

/// A baffled receiver delays each element by its range over the receipt's
/// delay speed. At the nearest element that delay must keep the receiver's
/// one-sample minimum after the filter's low-frequency delay is absorbed.
fn admit_receiver_delay(
    observer: &ObserverFilter,
    baffled: Option<([f64; 3], Option<GasState>)>,
    radius_m: f64,
    c0: f64,
) -> Result<(), ReedPerformanceError> {
    let Some((position, _)) = baffled else { return Ok(()) };
    let receipt = observer.receipt();
    let rate = f64::from(receipt.sample_rate_hz);
    let (near, _) = element_ranges(position, radius_m);
    let available_samples = near / c0 * rate;
    let needed_samples = near * receipt.low_frequency_delay_s / receipt.range_m * rate + 1.0;
    check_receiver_delay(receipt.range_m, available_samples, needed_samples, receipt.budget_energy)
        .map_err(ReedPerformanceError::Model)
}

/// The receiver-delay admission itself, on explicit numbers. A 5 cm receiver
/// at 48 kHz under a dispersive law refuses here; the budget sets where the
/// law's dispersion can be certified for each law, gas, clock and range.
fn check_receiver_delay(
    range_m: f64,
    available_samples: f64,
    needed_samples: f64,
    budget_energy: f64,
) -> Result<(), AcousticRealizeError> {
    if available_samples < needed_samples {
        return Err(AcousticRealizeError::ObserverLatency {
            range_m,
            available_samples,
            needed_samples,
            budget_energy,
        });
    }
    Ok(())
}

/// ISO 9613-1's admission of the observer medium and its relative humidity,
/// from the medium's own record: a legacy record's U.S. Standard Atmosphere
/// air, or a thermochemical record's composition. `None` for a named mixture
/// record, whose default law here is Stokes-Kirchhoff, or a resolved field.
fn iso_admission(
    record: Option<&GasRecord>,
    gas: &GasState,
) -> Result<Option<(fs_material::iso9613::TerrestrialAir, f64)>, AcousticRealizeError> {
    use fs_material::iso9613::{TerrestrialAir, iso9613_relative_humidity};
    match record {
        Some(GasRecord::Legacy(state)) => Ok(legacy_relative_humidity(state).and_then(|rh| {
            let (rebuilt, air) =
                TerrestrialAir::moist_air(state.temperature, state.pressure, rh).ok()?;
            (rebuilt == *gas).then_some((air, rh))
        })),
        Some(GasRecord::Thermochemical(resolved)) => {
            let medium = fs_material::relaxation::RelaxingMedium::from_thermochemical(resolved)
                .map_err(|error| AcousticRealizeError::Ambient(error.to_string()))?;
            let air = TerrestrialAir::admit(&medium).map_err(AcousticRealizeError::NonTerrestrialAir)?;
            let rh = iso9613_relative_humidity(gas.temperature, gas.pressure, gas.water_mole_fraction)
                .map_err(|error| AcousticRealizeError::Ambient(error.to_string()))?;
            Ok(Some((air, rh)))
        }
        Some(GasRecord::Mixture(..)) | None => Ok(None),
    }
}

/// The observer filter for the medium `gas`, as `acoustic_realize` chooses:
/// declared, or by default ISO 9613-1 for a legacy air record inside ISO's
/// window, otherwise Stokes-Kirchhoff.
fn observer_filter(
    declared: Option<DeclaredObserverLaw>,
    gas: &GasState,
    records: &[GasRecord],
    range_m: f64,
    rate: u32,
) -> Result<ObserverFilter, ReedPerformanceError> {
    let record = records.iter().find(|record| record.state() == gas);
    let filter = match declared {
        Some(DeclaredObserverLaw::Lossless) => ObserverFilter::lossless(gas, range_m, rate),
        Some(DeclaredObserverLaw::StokesKirchhoff) => {
            ObserverFilter::stokes_kirchhoff(gas, range_m, rate)
        }
        Some(DeclaredObserverLaw::Iso9613) => match iso_admission(record, gas) {
            Ok(Some((air, rh))) => ObserverFilter::iso9613(&air, gas, rh, range_m, rate),
            Ok(None) => Err(AcousticRealizeError::InvalidDescription {
                what: "observer_law iso9613 needs a legacy or thermochemical observer gas record",
            }),
            Err(error) => Err(error),
        },
        Some(DeclaredObserverLaw::IsoAirRelaxation) => match record {
            Some(GasRecord::Thermochemical(resolved)) => {
                fs_material::relaxation::RelaxingMedium::from_thermochemical(resolved)
                    .and_then(|medium| {
                        fs_material::relaxation::RelaxationAbsorption::iso9613_terrestrial_air()?
                            .evaluate(&medium)
                    })
                    .map_err(|error| AcousticRealizeError::Ambient(error.to_string()))
                    .and_then(|state| ObserverFilter::molecular_relaxation(&state, range_m, rate))
            }
            _ => Err(AcousticRealizeError::InvalidDescription {
                what: "observer_law relaxation needs a thermochemical observer gas record",
            }),
        },
        None => match iso_admission(record, gas) {
            Ok(Some((air, rh)))
                if matches!(record, Some(GasRecord::Legacy(_)))
                    && fs_material::iso9613::iso9613_absorption(&air, gas, rh, 1.0).is_ok() =>
            {
                ObserverFilter::iso9613(&air, gas, rh, range_m, rate)
            }
            _ => ObserverFilter::stokes_kirchhoff(gas, range_m, rate),
        },
    };
    filter.map_err(ReedPerformanceError::Model)
}

/// The absorption the listener gas applies over the observer's range up to
/// its band top, read from that gas's own record. A gas without one (a
/// resolved field cell) is unrecorded.
fn observer_omission(
    medium: &GasState,
    records: &[GasRecord],
    range_m: f64,
    band_top_hz: f64,
) -> Result<ObserverOmission, AcousticRealizeError> {
    match records.iter().find(|record| record.state() == medium) {
        Some(record) => omitted_observer_absorption(record, range_m, band_top_hz),
        None => Ok(ObserverOmission::unrecorded(medium, range_m, band_top_hz)),
    }
}

/// Source facts; neither a sourced-material receipt nor a fidelity certificate.
#[derive(Clone, Copy, Debug)]
pub struct ReedPerformanceInfo {
    /// Domain-separated hash of the complete exact input bytes.
    pub input_hash: ContentHash,
    /// Actual solver clock [Hz], independent of the observer's output clock.
    pub sample_rate_hz: u32,
    /// Finite half-open physical window [0, samples).
    pub samples: u64,
    /// Declared PCM scale [Pa]; never a per-part ensemble gain.
    pub full_scale_pa: f64,
    /// Static axial and side-branch records in their authored order.
    pub segments: usize,
    /// Side branches included in the same TMM reflectance.
    pub tone_holes: usize,
    /// Whether the source has nonzero reed mass and invokes massive mechanics.
    pub massive_reed: bool,
    /// Authored pressure commands, excluding the initial held value.
    pub gesture_events: usize,
    /// Compiled sample-addressed controls, including the initial value.
    pub compiled_controls: usize,
}

/// Explicit spatial image of a prescribed material field.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GasFieldResolution {
    /// Existing piecewise-constant midpoint cells per authored segment.
    Midpoint(usize),
    /// Shared fourth-order Magnus steps per smooth material/geometry span.
    Continuous(fs_duct::RiccatiSteps),
}

/// Finite pressure producer retaining the physical solver and schedule together.
///
/// It implements the common pressure interface directly, so decimation and
/// heterogeneous ensembles retain its horizon without extracting an unbounded
/// raw renderer. Admission never advances a sample. A late physics failure
/// inherits the existing renderer's poison-on-failure contract.
pub struct ReedPerformance {
    info: ReedPerformanceInfo,
    renderer: ScheduledRenderer,
    ambient_gas: GasState,
    segment_gases: Vec<GasState>,
    terminal: Option<outlet::TerminalObservation>,
    observation_policy: ObservationPolicy,
    gas_records: Vec<GasRecord>,
    gas_field: Option<(fs_material::gas_field::AxialGasField, GasFieldResolution)>,
    steady_gas: Option<SteadyGasPerformance>,
    observer: ObserverFilter,
    observer_range_spread_m: Option<f64>,
    observer_omission: Result<ObserverOmission, AcousticRealizeError>,
    observer_law_declared: bool,
}
impl ReedPerformance {
    /// Decode and admit all records before building or advancing the runtime.
    ///
    /// The schedule owns the entire pressure history, so the fixture's extra
    /// attack envelope is disabled. The six reed fields are opening, width,
    /// closing pressure, mass, stiffness and damping ratio. Zero mass/stiffness
    /// retain the owner's quasistatic/derived-stiffness limits. Four optional
    /// trailing fields give a measured lay of a massive reed: stiffness
    /// [N/m^alpha], exponent, onset opening [m] and Hunt-Crossley loss [s/m]; a
    /// six-field row keeps the default slit lay byte for byte. Static tone-hole
    /// openings are explicit fractions; they are never silently clamped here.
    ///
    /// # Errors
    /// Malformed/oversized input, unsupported or unobserved gestures, exhausted
    /// compilation budgets, or the existing physical owner's domain refusal.
    pub fn from_bytes(bytes: &[u8], max_block: usize) -> Result<Self, ReedPerformanceError> {
        Self::from_bytes_with_policy(bytes, max_block, ObservationPolicy::Diagnostic)
    }

    /// Select acceptance explicitly. Legacy listener_m is diagnostic-only;
    /// a baffled receiver retains the SAME flanged feedback and exterior gas.
    pub fn from_bytes_with_policy(
        bytes: &[u8],
        max_block: usize,
        policy: ObservationPolicy,
    ) -> Result<Self, ReedPerformanceError> {
        Self::from_bytes_with_gas_sources(bytes, max_block, policy, None)
    }

    /// Resolve ambient, segment, exterior and optional gas_field declarations
    /// against explicit caller-loaded sources. Original records, including the
    /// field knots, remain available for sidecars. A file-defined field uses
    /// ambient static pressure and conflicts with explicit segment gas records.
    pub fn from_bytes_with_gas_sources(
        bytes: &[u8],
        max_block: usize,
        policy: ObservationPolicy,
        sources: Option<&ThermochemicalSources<'_>>,
    ) -> Result<Self, ReedPerformanceError> {
        Self::from_bytes_with_environment(bytes, max_block, policy, sources, None)
    }

    /// Bind a source-retaining axial temperature/composition field to the same
    /// authored reed, geometry, gestures and physical terminal observation.
    /// This accepts both prescribed fields and converged steady breath fields.
    /// The exact inlet state drives the reed; the same regional duct drives
    /// feedback and the optional baffled microphone, including an explicitly
    /// authored exterior gas. No proxy gain or second acoustic line is added.
    ///
    /// The field must cover the complete authored length at the ambient static
    /// pressure. Explicit per-segment gas records conflict and refuse. The
    /// field remains frozen during this finite performance; mean-flow
    /// acoustics and time-dependent thermal transport are not implied.
    #[allow(clippy::too_many_arguments)]
    pub fn from_bytes_with_gas_field(
        bytes: &[u8],
        max_block: usize,
        policy: ObservationPolicy,
        sources: Option<&ThermochemicalSources<'_>>,
        field: fs_material::gas_field::AxialGasField,
        cells_per_segment: usize,
    ) -> Result<Self, ReedPerformanceError> {
        Self::from_bytes_with_field_resolution(
            bytes,
            max_block,
            policy,
            sources,
            field,
            GasFieldResolution::Midpoint(cells_per_segment),
        )
    }

    /// Select continuous Magnus or midpoint field propagation while retaining
    /// the SAME authored mechanics, pressure schedule, terminal load, and
    /// physical microphone. No observer or acoustic source is rebuilt through
    /// a separate material path. The field remains frozen during the render.
    #[allow(clippy::too_many_arguments)]
    pub fn from_bytes_with_field_resolution(
        bytes: &[u8],
        max_block: usize,
        policy: ObservationPolicy,
        sources: Option<&ThermochemicalSources<'_>>,
        field: fs_material::gas_field::AxialGasField,
        resolution: GasFieldResolution,
    ) -> Result<Self, ReedPerformanceError> {
        Self::from_bytes_with_environment(
            bytes,
            max_block,
            policy,
            sources,
            Some((field, resolution)),
        )
    }

    fn from_bytes_with_environment(
        bytes: &[u8],
        max_block: usize,
        policy: ObservationPolicy,
        sources: Option<&ThermochemicalSources<'_>>,
        mut gas_field: Option<(fs_material::gas_field::AxialGasField, GasFieldResolution)>,
    ) -> Result<Self, ReedPerformanceError> {
        if bytes.len() > MAX_REED_PERFORMANCE_BYTES || !(1..=65_536).contains(&max_block) {
            return Err(bad(
                0,
                "source exceeds 1 MiB or callback capacity is outside 1..=65536",
            ));
        }
        let text = std::str::from_utf8(bytes).map_err(|_| bad(0, "source must be UTF-8"))?;
        let (header, schedule_text) = text.split_once("\nschedule\n").ok_or_else(|| {
            bad(
                0,
                "expected an LF-delimited schedule record and canonical gesture bytes",
            )
        })?;
        let mut r = Reader {
            lines: header.lines(),
            line: 0,
        };
        r.row(REED_PERFORMANCE_SCHEMA)?.finish()?;
        let mut row = r.row("audio")?;
        let rate: u32 = row.parse()?;
        let samples: u64 = row.parse()?;
        let full_scale_pa = row.positive()?;
        row.finish()?;
        if !(1..=192_000).contains(&rate)
            || samples == 0
            || samples > 28_800_000
            || samples > u64::from(rate) * 600
        {
            return Err(bad(
                r.line,
                "audio requires a rate in 1..=192000 and at most 600 seconds/28800000 positive samples",
            ));
        }
        let mut row = r.row("ambient")?;
        let mut gas_records = Vec::new();
        let air = row.gas_with_sources(None, sources, &mut gas_records)?;
        row.finish()?;
        if gas_field
            .as_ref()
            .is_some_and(|(field, _)| field.static_pressure_pa() != air.pressure)
        {
            return Err(bad(
                r.line,
                "axial gas field must have the authored ambient static pressure",
            ));
        }
        let mut row = r.row("reed")?;
        let reed = BeatingReed {
            rest_opening_m: row.positive()?,
            width_m: row.positive()?,
            closing_pressure_pa: row.positive()?,
            mass_kg: row.scalar()?,
            stiffness_n_m: row.scalar()?,
            damping_ratio: row.scalar()?,
            blowing_pressure_pa: 0.0,
            attack_s: 0.0,
        };
        // Optional measured lay: stiffness [N/m^alpha], exponent, onset opening [m]
        // and Hunt-Crossley loss [s/m]; absent, the default slit lay applies.
        let lay = if row.more() {
            Some(ReedLay {
                stiffness: row.positive()?,
                alpha: row.scalar()?,
                onset_m: row.scalar()?,
                loss_s_per_m: row.scalar()?,
            })
        } else {
            None
        };
        row.finish()?;
        if !crate::reed_bore::reed_parameters_valid(reed) {
            return Err(bad(
                r.line,
                "reed mass, stiffness and damping must be finite and nonnegative",
            ));
        }
        if let Some(lay) = lay {
            if reed.mass_kg <= 0.0 {
                return Err(bad(r.line, "a measured lay needs a massive reed"));
            }
            if lay.alpha < 1.0 || lay.loss_s_per_m < 0.0 {
                return Err(bad(
                    r.line,
                    "lay exponent must be at least 1 and its loss nonnegative",
                ));
            }
            if lay.onset_m < 0.0 || lay.onset_m >= reed.rest_opening_m {
                return Err(bad(
                    r.line,
                    "lay onset must lie between zero and the reed's rest opening",
                ));
            }
        }
        let slit_length_m = if r.lines.clone().next()
            .and_then(|line| line.split_ascii_whitespace().next()) == Some("slit_length") {
            let mut row = r.row("slit_length")?;
            let length = row.positive()?;
            row.finish()?;
            Some(length)
        } else { None };
        let tongue_contact = if r.lines.clone().next()
            .and_then(|line| line.split_ascii_whitespace().next()) == Some("tongue_contact") {
            let mut row = r.row("tongue_contact")?;
            let area_m2 = row.positive()?;
            let effective_area_m2 = row.positive()?;
            let force_n = row.positive()?;
            let release_start_s = row.scalar()?;
            let release_s = row.positive()?;
            row.finish()?;
            if reed.mass_kg <= 0.0 {
                return Err(bad(r.line, "tongue_contact requires a massive reed"));
            }
            if release_start_s < 0.0 {
                return Err(bad(
                    r.line,
                    "tongue_contact requires a nonnegative release_start_s",
                ));
            }
            let modal_factor = cantilever_tip_modal_factor(area_m2, effective_area_m2)
                .ok_or_else(|| bad(r.line, "tongue_contact area exceeds the vibrating reed blade area"))?;
            Some(TongueContactInput {
                modal_factor,
                force_n,
                release_start_s,
                release_s,
            })
        } else {
            None
        };
        let observation =
            outlet::ObservationInput::read_with_gas_sources(&mut r, sources, &mut gas_records)?;
        observation.validate_policy(policy)?;
        let declared_law = if r.lines.clone().next()
            .and_then(|line| line.split_ascii_whitespace().next()) == Some("observer_law") {
            let mut row = r.row("observer_law")?;
            let law = match row.word()? {
                "none" => DeclaredObserverLaw::Lossless,
                "iso9613" => DeclaredObserverLaw::Iso9613,
                "stokes-kirchhoff" => DeclaredObserverLaw::StokesKirchhoff,
                "relaxation" => {
                    if row.word()? != "iso9613-terrestrial-air" {
                        return Err(bad(r.line, "observer relaxation data must be iso9613-terrestrial-air"));
                    }
                    DeclaredObserverLaw::IsoAirRelaxation
                }
                _ => return Err(bad(
                    r.line,
                    "observer_law must be none, iso9613, stokes-kirchhoff or relaxation iso9613-terrestrial-air",
                )),
            };
            row.finish()?;
            Some(law)
        } else {
            None
        };
        let observer_range_m = observation.range_m();
        let baffled = observation.baffled();
        let mut row = r.row("termination")?;
        let termination = match row.word()? {
            "closed" => Termination::Closed,
            "ideal-open" => Termination::IdealOpen,
            "unflanged" => Termination::UnflangedOpen,
            "flanged" => Termination::FlangedOpen,
            _ => return Err(bad(r.line, "unsupported static duct termination")),
        };
        row.finish()?;
        observation.validate_termination(termination, r.line)?;
        let mut row = r.row("segments")?;
        let count = row.count(MAX_SEGMENTS)?;
        row.finish()?;
        if count == 0 {
            return Err(bad(r.line, "at least one axial segment is required"));
        }
        let mut segments = Vec::with_capacity(count);
        let mut segment_gases = Vec::with_capacity(count);
        let mut explicit_segment_gas = false;
        let mut holes = 0;
        for _ in 0..count {
            let mut row = r.next()?;
            let segment = match row.word()? {
                "cylinder" => Segment::Cylinder {
                    radius: row.positive()?,
                    length: row.positive()?,
                },
                "cone" => Segment::Cone {
                    inlet_radius: row.positive()?,
                    outlet_radius: row.positive()?,
                    length: row.positive()?,
                },
                "hole" => {
                    let hole_radius = row.positive()?;
                    let chimney_height = row.positive()?;
                    let bore_radius = row.positive()?;
                    let sigma = row.scalar()?;
                    if hole_radius >= bore_radius || !(0.0..=1.0).contains(&sigma) {
                        return Err(bad(
                            r.line,
                            "hole radius must be below bore radius and opening must be in [0,1]",
                        ));
                    }
                    holes += 1;
                    let state = if sigma == 0.0 {
                        HoleState::Closed
                    } else if sigma == 1.0 {
                        HoleState::Open
                    } else {
                        HoleState::Vent(sigma)
                    };
                    Segment::ToneHole {
                        hole_radius,
                        chimney_height,
                        bore_radius,
                        state,
                    }
                }
                _ => return Err(bad(r.line, "expected cylinder, cone or hole geometry")),
            };
            let gas = match row.words.next() {
                None => air,
                Some("gas") => {
                    if gas_field.is_some() {
                        return Err(bad(
                            row.line,
                            "explicit segment gas conflicts with the complete axial gas field",
                        ));
                    }
                    explicit_segment_gas = true;
                    row.gas_with_sources(Some(air.pressure), sources, &mut gas_records)?
                }
                Some(_) => {
                    return Err(bad(
                        row.line,
                        "expected an optional gas record after geometry",
                    ));
                }
            };
            if gas.pressure != air.pressure {
                return Err(bad(
                    row.line,
                    "regional duct gas must have the ambient static pressure; mean flow is not modeled",
                ));
            }
            row.finish()?;
            segments.push(segment);
            segment_gases.push(gas);
        }
        if matches!(segments.first(), Some(Segment::ToneHole { .. }))
            || matches!(segments.last(), Some(Segment::ToneHole { .. }))
        {
            return Err(bad(
                r.line,
                "the inlet and termination must have explicit axial segments",
            ));
        }
        let mut steady_gas = None;
        let field_kind = r
            .lines
            .clone()
            .next()
            .and_then(|line| line.split_ascii_whitespace().next());
        if matches!(field_kind, Some("gas_field" | "gas_steady")) {
            if gas_field.is_some() {
                return Err(bad(
                    r.line + 1,
                    "file gas field conflicts with the injected axial gas field",
                ));
            }
            if explicit_segment_gas {
                return Err(bad(
                    r.line + 1,
                    "file gas field conflicts with explicit segment gas records",
                ));
            }
            if field_kind == Some("gas_steady") {
                let (generated, resolution) = field_input::steady::read(
                    &mut r,
                    sources,
                    air.pressure,
                    &segments,
                    &mut gas_records,
                )?;
                gas_field = Some((generated.solution().field().clone(), resolution));
                steady_gas = Some(generated);
            } else {
                gas_field = Some(field_input::read(
                    &mut r,
                    sources,
                    air.pressure,
                    &mut gas_records,
                )?);
            }
        }
        let molecular = if r.lines.clone().next()
            .and_then(|line| line.split_ascii_whitespace().next()) == Some("molecular") {
            if explicit_segment_gas {
                return Err(bad(r.line + 1,
                    "molecular bore propagation requires homogeneous ambient species or a continuous field, no segment gas"));
            }
            if matches!(&gas_field, Some((_, GasFieldResolution::Midpoint(_)))) {
                return Err(bad(r.line + 1,
                    "molecular gas fields require explicit continuous spatial resolution, not midpoint cells"));
            }
            Some(molecular_input::read(
                &mut r, &gas_records[0], gas_field.as_ref().map(|(field, _)| field),
            )?)
        } else { None };
        let mut row = r.row("compile_limits")?;
        let max_work: u64 = row.parse()?;
        let max_controls = row.count(MAX_CONTROLS)?;
        row.finish()?;
        if max_work > MAX_WORK {
            return Err(bad(r.line, "compilation exceeds 16777216 source visits"));
        }
        if r.lines.next().is_some() {
            return Err(bad(r.line + 1, "unexpected trailing header record"));
        }
        let schedule = decode_schedule(schedule_text)?;
        let [track] = schedule.tracks() else {
            return Err(bad(0, "exactly one pressure track is required"));
        };
        if fs_scenario::gesture::path::legacy_member(&track.target, "blowing_pressure") != Some(0)
            || schedule.control_rate_hz == 0
            || schedule.control_rate_hz > rate
        {
            return Err(bad(
                0,
                "only blowing pressure is supported, at a control rate no greater than mechanics",
            ));
        }
        let last_tick = (samples - 1) * u64::from(schedule.control_rate_hz) / u64::from(rate);
        // One track emits at most one control per tick. Bound the compiler
        // before it can allocate, including long constant-valued schedules.
        if last_tick >= MAX_CONTROLS as u64 {
            return Err(bad(
                0,
                "pressure compilation exceeds the 262144 control-tick budget",
            ));
        }
        let last_time = last_tick as f64 / f64::from(schedule.control_rate_hz);
        if track.events.iter().any(|event| event.time_s > last_time) {
            return Err(bad(
                0,
                "an authored pressure command starts after the last observed control tick",
            ));
        }
        let mut controls = compile_pressure_gestures(
            &schedule,
            &[PressureGestureBinding {
                track: track.id.clone(),
                voice: 0,
            }],
            rate,
            samples,
            max_work,
        )
        .map_err(ReedPerformanceError::Schedule)?;
        if let Some(tongue) = tongue_contact {
            if tongue.release_start_s > last_time {
                return Err(bad(
                    0,
                    "tongue_contact release starts after the last observed control tick",
                ));
            }
            tongue.append_controls(
                &mut controls,
                &schedule,
                &track.id,
                reed,
                rate,
                last_tick,
            )?;
        }
        if controls.len() > max_controls {
            return Err(bad(
                0,
                "compiled pressure controls exceed the declared storage budget",
            ));
        }
        let info = ReedPerformanceInfo {
            input_hash: hash_domain("org.frankensim.fs-couple.reed-performance.v1", bytes),
            sample_rate_hz: rate,
            samples,
            full_scale_pa,
            segments: count,
            tone_holes: holes,
            massive_reed: reed.mass_kg > 0.0,
            gesture_events: track.events.len(),
            compiled_controls: controls.len(),
        };
        let duct = Duct { segments };
        let outlet_radius_m = duct.segments.last().map_or(0.0, Segment::outlet_radius);
        let band_top_hz = observer_band_top_hz(observation.declared_band_hz(), rate);
        let (voice, terminal, observer, observer_omission) = if let Some(molecular) = molecular {
            let outlet_gas = match &gas_field {
                Some((field, _)) => *field.state_at(field.length_m())
                    .map_err(|source| ReedPerformanceError::Gas { line: 0, source })?.state(),
                None => air,
            };
            let medium = observer_medium(baffled, &air, &outlet_gas);
            let exterior_gas = baffled.and_then(|(_, exterior)| exterior).unwrap_or(outlet_gas);
            let observer =
                observer_filter(declared_law, &medium, &gas_records, observer_range_m, rate)?;
            admit_receiver_delay(&observer, baffled, outlet_radius_m, medium.sound_speed)?;
            let omission = observer_omission(&medium, &gas_records, observer_range_m, band_top_hz);
            let point = match &gas_field {
                Some((field, GasFieldResolution::Continuous(steps))) => {
                    let point = molecular.bind_field(&duct, field, *steps, exterior_gas)?;
                    // Read-only metadata uses the same authored-segment inlet
                    // samples as the continuous equilibrium field route. The
                    // actual molecular operator keeps its local quadrature.
                    let mut position = 0.0;
                    segment_gases = duct.segments.iter().map(|segment| {
                        let sample = field.state_at(position)
                            .map_err(|source| ReedPerformanceError::Gas { line: 0, source })?;
                        position += segment.axial_length();
                        Ok(*sample.state())
                    }).collect::<Result<_, ReedPerformanceError>>()?;
                    point
                }
                None => molecular.bind(&duct, exterior_gas)?,
                Some((_, GasFieldResolution::Midpoint(_))) => {
                    unreachable!("molecular field resolution admitted before compilation")
                }
            };
            let (voice, terminal) = observation.build_molecular(
                point, reed, termination, rate, samples as usize,
                observer.receipt().delay_speed_m_s.unwrap_or(medium.sound_speed),
                declared_law, &gas_records,
            )?;
            (voice, terminal, observer, omission)
        } else if let Some((field, resolution)) = &gas_field {
            use crate::driving_point::regional::{RegionalDuct, field::GasFieldDuctError};
            let regional = match resolution {
                GasFieldResolution::Midpoint(cells) => {
                    RegionalDuct::from_gas_field(&duct, field, *cells)
                }
                GasFieldResolution::Continuous(steps) => {
                    RegionalDuct::from_continuous_gas_field(&duct, field, *steps)
                }
            }
            .map_err(|error| match error {
                GasFieldDuctError::Material(source) => {
                    ReedPerformanceError::Gas { line: 0, source }
                }
                GasFieldDuctError::Duct(source) => {
                    ReedPerformanceError::Model(AcousticRealizeError::Duct(source))
                }
            })?;
            segment_gases = regional.segment_gases().to_vec();
            let medium = observer_medium(baffled, &air, regional.outlet_gas());
            let observer =
                observer_filter(declared_law, &medium, &gas_records, observer_range_m, rate)?;
            admit_receiver_delay(&observer, baffled, outlet_radius_m, medium.sound_speed)?;
            let omission = observer_omission(&medium, &gas_records, observer_range_m, band_top_hz);
            let (voice, terminal) = observation.build_regional(
                &duct,
                regional,
                reed,
                termination,
                rate,
                samples as usize,
                observer.receipt().delay_speed_m_s.unwrap_or(medium.sound_speed),
            )?;
            (voice, terminal, observer, omission)
        } else {
            let outlet_gas = *segment_gases.last().expect("at least one admitted segment");
            let medium = observer_medium(baffled, &air, &outlet_gas);
            let observer =
                observer_filter(declared_law, &medium, &gas_records, observer_range_m, rate)?;
            admit_receiver_delay(&observer, baffled, outlet_radius_m, medium.sound_speed)?;
            let omission = observer_omission(&medium, &gas_records, observer_range_m, band_top_hz);
            let (voice, terminal) = observation.build(
                &duct,
                &segment_gases,
                reed,
                termination,
                rate,
                samples as usize,
                observer.receipt().delay_speed_m_s.unwrap_or(medium.sound_speed),
            )?;
            (voice, terminal, observer, omission)
        };
        let observer_range_spread_m =
            baffled.map(|(position, _)| range_spread(position, outlet_radius_m));
        let mut voice = voice;
        if let Some(lay) = lay {
            voice.set_lay(lay).map_err(ReedPerformanceError::Model)?;
        }
        if let Some(length_m) = slit_length_m {
            voice.set_viscous_slit_length(length_m).map_err(ReedPerformanceError::Model)?;
        }
        let context = RenderContext::new(vec![RenderVoice::ReedBore(voice)], max_block);
        let renderer = ScheduledRenderer::new(context, controls, max_controls)
            .map_err(ReedPerformanceError::Render)?;
        renderer
            .validate_sample_rate(rate)
            .map_err(ReedPerformanceError::Render)?;
        renderer
            .validate_sample_count(samples)
            .map_err(ReedPerformanceError::Render)?;
        Ok(Self {
            info,
            renderer,
            ambient_gas: air,
            segment_gases,
            terminal,
            observation_policy: policy,
            gas_records,
            gas_field,
            steady_gas,
            observer,
            observer_range_spread_m,
            observer_omission,
            observer_law_declared: declared_law.is_some(),
        })
    }

    /// Original source facts, unchanged by output-clock conversion.
    #[must_use]
    pub const fn info(&self) -> ReedPerformanceInfo {
        self.info
    }
    /// Admission scope; diagnostic output must not be used as acceptance audio.
    pub const fn observation_policy(&self) -> ObservationPolicy {
        self.observation_policy
    }
    /// The observer absorption law applied to the emitted pressure.
    #[must_use]
    pub const fn observer(&self) -> ReedObserverReceipt {
        let filter = *self.observer.receipt();
        ReedObserverReceipt {
            filter,
            range_spread_m: self.observer_range_spread_m,
            consumer_delay_s: match self.observer_range_spread_m {
                None => Some(filter.receiver_delay_s),
                Some(_) => None,
            },
            internal_part_unpropagated: self.observer_range_spread_m.is_none(),
        }
    }
    /// The absorption the listener gas applies over the observer's range up to
    /// its band top, and what the applied law leaves of it (rulings 6510,
    /// 6532). ISO 9613-1 and the relaxation law leave nothing out;
    /// Stokes-Kirchhoff, the default for a mixture, a thermochemical gas or a
    /// field, leaves out the gas's relaxation; `none` leaves out everything.
    /// The library only reports it: a sidecar writer refuses through
    /// [`ObserverOmission::admit`], which this applies for an undeclared
    /// Stokes-Kirchhoff law. A declared law always passes.
    ///
    /// # Errors
    /// The figure's own refusal (a state outside its law's window), or
    /// [`AcousticRealizeError::ObserverOmission`].
    pub fn observer_admission(&self) -> Result<&ObserverOmission, AcousticRealizeError> {
        let omission = self.observer_omission.as_ref().map_err(Clone::clone)?;
        if matches!(self.observer.receipt().absorption.law, GasAbsorptionModel::StokesKirchhoff) {
            omission.admit(true, self.observer_law_declared)?;
        }
        if let Some(terminal) = &self.terminal {
            terminal.admit_molecular_hole_observers(self.observer_law_declared)?;
        }
        Ok(omission)
    }
    /// Whether the file declares its observer law (`observer_law ...`).
    #[must_use]
    pub const fn observer_law_declared(&self) -> bool {
        self.observer_law_declared
    }
    /// Read-only access to sample-addressed controls and retained state clocks.
    #[must_use]
    pub const fn renderer(&self) -> &ScheduledRenderer {
        &self.renderer
    }
    /// Actual selected channel, retaining source-gas viscosity and density.
    /// The field inlet can differ from the first propagation cell's midpoint.
    pub fn viscous_slit(&self) -> Option<&crate::bernoulli_aperture::viscous::ViscousSlit> {
        self.renderer.context().reed_voice(0).and_then(ReedBoreVoice::viscous_slit)
    }
    /// Original molecular gas, wall and geometry owner used for both bore
    /// feedback and physical outlet pressure. None preserves the existing
    /// equilibrium regional or axial-field path.
    #[must_use]
    pub fn molecular_driving_point(
        &self,
    ) -> Option<&crate::driving_point::molecular::MolecularDrivingPoint> {
        self.renderer.context().reed_voice(0)
            .and_then(ReedBoreVoice::molecular_driving_point)
    }
    /// The immutable default medium. Segments without `gas` inherit this state;
    /// explicit segment records may override it. No air-only absorption is added.
    #[must_use]
    pub const fn ambient_gas(&self) -> &GasState {
        &self.ambient_gas
    }
    /// Actual propagation gases, in authored segment order for regional input
    /// or resolved midpoint-cell order for a midpoint axial field. Continuous
    /// fields expose authored-segment inlet samples; their operator evaluates
    /// propagation at its exact local quadrature points. The field's exact
    /// inlet state drives reed density/impedance; it can differ from the first
    /// cell's midpoint. An explicit exterior separately supplies radiation
    /// loading and observation. All are stationary.
    #[must_use]
    pub fn segment_gases(&self) -> &[GasState] {
        &self.segment_gases
    }
    /// Source-retaining axial field and its explicit resolution count.
    /// See gas_field_resolution for the selected midpoint/continuous meaning.
    /// None selects the unchanged uniform/piecewise-regional input route.
    #[must_use]
    pub fn gas_field(&self) -> Option<(&fs_material::gas_field::AxialGasField, usize)> {
        self.gas_field.as_ref().map(|(field, resolution)| {
            (
                field,
                match resolution {
                    GasFieldResolution::Midpoint(cells) => *cells,
                    GasFieldResolution::Continuous(steps) => steps.0,
                },
            )
        })
    }
    /// The selected field integration image, independent of output block size.
    pub fn gas_field_resolution(&self) -> Option<GasFieldResolution> {
        self.gas_field.as_ref().map(|(_, resolution)| *resolution)
    }
    /// The admitted stationary heat/species problem and its actual convergence
    /// evidence, when gas_steady generated the field used by this renderer.
    /// Prescribed fields and legacy regional input return None.
    pub fn steady_gas(&self) -> Option<&SteadyGasPerformance> {
        self.steady_gas.as_ref()
    }
    /// All authored gas records, including field knots and source receipts, in
    /// read order. The field model separately retains every declared source ID.
    #[must_use]
    pub fn gas_records(&self) -> &[GasRecord] {
        &self.gas_records
    }
    /// Remaining samples on the mechanical clock.
    #[must_use]
    pub fn remaining_samples(&self) -> u64 {
        self.info.samples - self.renderer.samples_rendered()
    }
}
impl PressureRenderer for ReedPerformance {
    fn samples_rendered(&self) -> u64 {
        self.renderer.samples_rendered()
    }
    fn max_block_len(&self) -> usize {
        self.renderer.context().max_block_len()
    }
    fn validate_sample_rate(&self, rate: u32) -> Result<(), RenderError> {
        self.renderer.validate_sample_rate(rate)
    }
    fn validate_sample_count(&self, samples: u64) -> Result<(), RenderError> {
        if self.last_outlet_error().is_some() {
            return Err(RenderError::Poisoned);
        }
        self.renderer.validate_sample_count(samples)?;
        if samples > self.remaining_samples() {
            return Err(RenderError::Sizing {
                what: "requested pressure window exceeds the reed performance horizon",
            });
        }
        Ok(())
    }
    fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        // The older shared renderer's error envelope predates terminal observers.
        // The typed owner error remains available through pressure_block and
        // last_outlet_error; this compatibility path retains its full diagnostic.
        self.pressure_block(output).map_err(|error| match error {
            ReedPerformanceError::Render(error) => error,
            error => RenderError::Voice(AcousticRealizeError::Nonlinear(error.to_string())),
        })
    }
}

fn decode_schedule(text: &str) -> Result<GestureSchedule, ReedPerformanceError> {
    // The shared decoder reserves from declared counts, including trailing ones.
    let lines = text.lines().count();
    for line in text.lines() {
        if let Some(value) = line.strip_prefix("tracks\t") {
            if value.parse::<usize>().ok() != Some(1) {
                return Err(bad(0, "exactly one canonical pressure track is required"));
            }
        }
        if let Some(value) = line.strip_prefix("events\t") {
            let count = value
                .parse::<usize>()
                .map_err(|_| bad(0, "invalid gesture event count"))?;
            if count > MAX_EVENTS || count > lines {
                return Err(bad(0, "gesture count exceeds source/16384-event budget"));
            }
        }
    }
    let schedule = GestureSchedule::from_canonical_bytes(text.as_bytes())
        .map_err(|e| ReedPerformanceError::Schedule(GestureCompileError::Gesture(e)))?;
    if schedule.to_canonical_bytes() != text.as_bytes() {
        return Err(bad(
            0,
            "gesture bytes must be canonical without ignored fields or trailing records",
        ));
    }
    Ok(schedule)
}
struct Reader<'a> {
    lines: Lines<'a>,
    line: usize,
}
struct Row<'a> {
    words: SplitAsciiWhitespace<'a>,
    line: usize,
}
impl<'a> Reader<'a> {
    fn next(&mut self) -> Result<Row<'a>, ReedPerformanceError> {
        self.line += 1;
        let words = self
            .lines
            .next()
            .ok_or_else(|| bad(self.line, "missing header record"))?
            .split_ascii_whitespace();
        Ok(Row {
            words,
            line: self.line,
        })
    }
    fn row(&mut self, key: &str) -> Result<Row<'a>, ReedPerformanceError> {
        let mut row = self.next()?;
        if row.word()? != key {
            return Err(bad(self.line, "unexpected header record"));
        }
        Ok(row)
    }
}
impl<'a> Row<'a> {
    fn gas(&mut self, pressure: Option<f64>) -> Result<GasState, ReedPerformanceError> {
        self.gas_with_sources(pressure, None, &mut Vec::new())
    }
    fn gas_with_sources(
        &mut self,
        pressure: Option<f64>,
        sources: Option<&ThermochemicalSources<'_>>,
        records: &mut Vec<GasRecord>,
    ) -> Result<GasState, ReedPerformanceError> {
        let record =
            read_gas_record(&mut self.words, pressure, sources).map_err(|error| match error {
                GasInputError::Input(what) => bad(self.line, what),
                GasInputError::Material(source) => ReedPerformanceError::Gas {
                    line: self.line,
                    source,
                },
            })?;
        let gas = *record.state();
        records.push(record);
        Ok(gas)
    }
    fn more(&self) -> bool {
        self.words.clone().next().is_some()
    }
    fn word(&mut self) -> Result<&'a str, ReedPerformanceError> {
        self.words
            .next()
            .ok_or_else(|| bad(self.line, "missing field"))
    }
    fn parse<T: FromStr>(&mut self) -> Result<T, ReedPerformanceError> {
        self.word()?
            .parse()
            .map_err(|_| bad(self.line, "malformed field"))
    }
    fn scalar(&mut self) -> Result<f64, ReedPerformanceError> {
        let x: f64 = self.parse()?;
        if !x.is_finite() {
            return Err(bad(self.line, "physical scalars must be finite"));
        }
        Ok(x)
    }
    fn positive(&mut self) -> Result<f64, ReedPerformanceError> {
        let x = self.scalar()?;
        if x <= 0.0 {
            return Err(bad(self.line, "physical dimension/scale must be positive"));
        }
        Ok(x)
    }
    fn count(&mut self, cap: usize) -> Result<usize, ReedPerformanceError> {
        let n = self.parse()?;
        if n > cap {
            return Err(bad(self.line, "count exceeds the input budget"));
        }
        Ok(n)
    }
    fn finish(mut self) -> Result<(), ReedPerformanceError> {
        if self.words.next().is_some() {
            return Err(bad(self.line, "unexpected extra field"));
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug)]
struct TongueContactInput {
    modal_factor: f64,
    force_n: f64,
    release_start_s: f64,
    release_s: f64,
}

impl TongueContactInput {
    fn append_controls(
        &self,
        controls: &mut Vec<ScheduledControl>,
        schedule: &GestureSchedule,
        track_id: &str,
        reed: BeatingReed,
        audio_rate_hz: u32,
        last_tick: u64,
    ) -> Result<(), ReedPerformanceError> {
        let control_rate_hz = schedule.control_rate_hz;
        let face = crate::reed_bore::reed_pressure_face(reed);
        let release_tick = ((self.release_start_s * f64::from(control_rate_hz)).round() as u64)
            .min(last_tick);
        let release_pressure_pa = schedule
            .sample(track_id, release_tick)
            .map_err(|e| ReedPerformanceError::Schedule(GestureCompileError::Gesture(e)))?;
        let pre_release_close_n = self.modal_factor * self.force_n;
        let held_downward_n = pre_release_close_n + face * release_pressure_pa;
        let mut previous = None;
        for tick in 0..=last_tick {
            let t = tick as f64 / f64::from(control_rate_hz);
            let force_n = if t <= self.release_start_s {
                let p_m = schedule
                    .sample(track_id, tick)
                    .map_err(|e| ReedPerformanceError::Schedule(GestureCompileError::Gesture(e)))?;
                (held_downward_n - face * p_m).max(0.0)
            } else {
                let elapsed = t - self.release_start_s;
                if elapsed < self.release_s {
                    (pre_release_close_n * (1.0 - elapsed / self.release_s)).max(0.0)
                } else {
                    0.0
                }
            };
            let bits = force_n.to_bits();
            if previous != Some(bits) {
                let sample = (u128::from(tick) * u128::from(audio_rate_hz))
                    .div_ceil(u128::from(control_rate_hz)) as u64;
                controls
                    .try_reserve(1)
                    .map_err(|_| bad(0, "cannot allocate tongue_contact controls"))?;
                controls.push(ScheduledControl {
                    sample,
                    delta: ControlDelta::SetReedClosingForce { voice: 0, force_n },
                });
                previous = Some(bits);
            }
        }
        Ok(())
    }
}

/// Fundamental clamped-free Euler-Bernoulli cantilever mode projection factor
/// $\Phi_{\text{tongue}}(\eta) = \frac{1}{\eta} \int_{1-\eta}^1 \psi_1(\xi)\,d\xi$
/// for a uniform tongue contact strip of area $A_{\text{tongue}}$ at the reed tip,
/// where $S_r = c_1 A_{\text{reed}}$ is the effective pressure area of the reed
/// (`Chatziioannou & van Walstijn 2009` Eq. 6; `Ducasse 2003` p. 62).
fn cantilever_tip_modal_factor(area_m2: f64, effective_area_m2: f64) -> Option<f64> {
    const BETA1: f64 = 1.875_104_068_711_961_2;
    const SIGMA1: f64 = 0.734_095_513_758_888_6;
    const C1: f64 = 2.0 * SIGMA1 / BETA1;
    let eta = C1 * area_m2 / effective_area_m2;
    if !eta.is_finite() || eta <= 0.0 || eta > 1.0 {
        return None;
    }
    let psi_integral = |xi: f64| -> f64 {
        let bxi = BETA1 * xi;
        ((bxi.sinh() - bxi.sin()) - SIGMA1 * (bxi.cosh() + bxi.cos() - 2.0)) / (2.0 * BETA1)
    };
    let factor = (psi_integral(1.0) - psi_integral(1.0 - eta)) / eta;
    (factor.is_finite() && factor > 0.0).then_some(factor)
}

#[cfg(test)]
impl ReedPerformance {
    /// A fresh filter for this performance's observer law, rebuilt from its
    /// receipt alone, to apply the same law to an unabsorbed reference.
    pub(super) fn reference_observer(&self, medium: &GasState) -> ObserverFilter {
        let filter = self.observer().filter;
        ObserverFilter::for_law(filter.absorption.law, medium, None, filter.range_m, filter.sample_rate_hz)
            .expect("the receipt's law rebuilds")
    }
}

#[cfg(test)]
mod gas_tests;

#[cfg(test)]
mod regional_tests;

#[cfg(test)]
mod field_tests;

#[cfg(test)]
mod segment_cap_tests;

#[cfg(test)]
mod lay_tests;

#[cfg(test)]
mod viscous_tests;

