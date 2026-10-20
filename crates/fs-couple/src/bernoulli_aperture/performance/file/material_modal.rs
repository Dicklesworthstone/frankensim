//! Material-card modal walls, two real duct networks, and exterior pressure.
//!
//! This source adapter has no instrument presets. Each blade retains its
//! original regional material chart, every requested signed mode, pressure
//! traction and full selected viscosity matrix. Shared owners perform all
//! geometry, flow, nonlinear contact, propagation, radiation and scheduling.
//! Flat DKT walls and the effective uniform slit/contact remain approximations;
//! source provenance and matched baffled observation do not qualify a recording.
//! An optional `propagator midpoint|exact-zoh` record after `audio` selects the
//! structural clock image before material modes are admitted. Omission retains
//! midpoint. ExactZoh keeps the original material M/K/full C, with pressure and
//! contact held per mechanical sample; it does not certify nonlinear time error.
mod plate;

use super::{PlateValveInputError, Reader, bad, checkpoint, decode_schedule, gas,
    read_receiver, MAX_PLATE_VALVE_PERFORMANCE_BYTES, MAX_CONTROLS};
use super::duct::DuctInput;
use crate::bernoulli_aperture::fractional::FractionalTubePolicy;
use crate::bernoulli_aperture::modal::ModalAperture;
use crate::bernoulli_aperture::paired::{BladePropagator, BladeState, PairedApertureLimits};
use crate::bernoulli_aperture::performance::AperturePerformanceConfig;
use crate::bernoulli_aperture::radiation::environment::BaffledOutletBinding;
use crate::bernoulli_aperture::slit_channel::{SlitChannel, SlitChannelSpec, SlitExit};
use crate::bernoulli_aperture::two_sided::{TwoSidedLimits, combined_inlet_state};
use crate::bernoulli_aperture::two_sided::supply::combined_inlet_state_with_supply;
use crate::bernoulli_aperture::two_sided::observer::TwoSidedPressure;
use crate::bernoulli_aperture::two_sided::observer::schedule::ScheduledTwoSidedPressure;
use crate::pcm_wav::baffled::openings::duct::{DuctOpening, DuctRadiationSpec, DuctReconstruction};
use crate::render::schedule::GestureCompileError;
use fs_blake3::{ContentHash, hash_domain};
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;
use fs_matdb::NormalizedMaterialCardPack;
use fs_material::gas::GasState;
use fs_material::gas_species::input::{GasRecord, ThermochemicalSources};
use plate::{MaterialCache, LossSelection};

pub use plate::{ModalValveMaterialRecord, MAX_MODAL_VALVE_PACK_BYTES};

/// Exact grammar selector. Filenames and instrument names never select physics.
pub const MATERIAL_MODAL_VALVE_SCHEMA: &str = "frankensim-material-modal-valve-v1";

/// Physical source facts, independent of callback partition and PCM scaling.
#[derive(Clone, Copy, Debug)]
pub struct MaterialModalValveInfo {
    /// Exact performance bytes; normalized material artifacts are separate below.
    pub input_hash: ContentHash,
    /// Actual mechanical clock [Hz].
    pub mechanics_rate_hz: u32,
    /// Original midpoint image or exact held-force material structural image.
    pub propagator: BladePropagator,
    /// Output microphone clock [Hz].
    pub output_rate_hz: u32,
    /// Finite output frame count, including the authored release/ringdown window.
    pub samples: u64,
    /// Fixed pascals corresponding to PCM full scale.
    pub full_scale_pa: f64,
    /// Number of independent source plates, one or two.
    pub blades: usize,
    /// Complete retained modal coordinate count across both source plates.
    pub modes: usize,
    /// Exact original mesh nodes, summed across source plates.
    pub nodes: usize,
    /// Original triangles, summed across source plates.
    pub triangles: usize,
}

/// File-built material valve retaining the original scheduled physical system.
pub struct MaterialModalValvePerformance {
    info: MaterialModalValveInfo,
    renderer: ScheduledTwoSidedPressure,
    materials: Vec<ModalValveMaterialRecord>,
    gases: Vec<GasRecord>,
    ambient: GasState,
    radiation: Vec<BaffledOutletBinding>,
}

impl MaterialModalValvePerformance {
    /// Resolve every named material artifact through the caller, with no I/O in
    /// the library. Cache repeated paths as one immutable original artifact.
    /// All acoustic media use the existing gas-record decoder and optional
    /// thermochemical source registry. Material T/omega are separately declared;
    /// ambient temperature never rewrites a solid's source state.
    ///
    /// Builds the valve against BOTH actual inlet impedances, including the
    /// optional retained no-slip supply passage. The single common exterior
    /// declaration installs each compact radiation load exactly once and drives
    /// the pressure observer. Source events execute at mechanical sample times.
    ///
    /// # Errors
    /// Missing/ambiguous/out-of-domain source properties, malformed or oversized
    /// records, incomplete geometry, inadmissible actual modes or full damping,
    /// gas/clock/load mismatch, nonlinear/observation limits, or cancellation.
    /// No model is returned with dropped modes, fallback coefficients or an
    /// internal-pressure substitute for the requested exterior microphone.
    pub fn from_bytes_with_sources(
        bytes: &[u8],
        max_block: usize,
        gate: &CancelGate,
        gas_sources: Option<&ThermochemicalSources<'_>>,
        mut load_material: impl FnMut(&str) -> Result<NormalizedMaterialCardPack, String>,
    ) -> Result<Self, PlateValveInputError> {
        checkpoint(gate)?;
        if bytes.len() > MAX_PLATE_VALVE_PERFORMANCE_BYTES || !(1..=65536).contains(&max_block) {
            return Err(bad(0, "material-modal input exceeds 4 MiB or callback capacity is outside 1..=65536"));
        }
        let text = core::str::from_utf8(bytes).map_err(|_| bad(0, "input must be UTF-8"))?;
        let (header, gestures) = text.split_once("\nschedule\n")
            .ok_or_else(|| bad(0, "missing LF-delimited schedule and canonical gesture suffix"))?;
        let mut rows = Reader { lines: header.lines(), line: 0 };
        rows.row(MATERIAL_MODAL_VALVE_SCHEMA)?.finish()?;
        let mut row = rows.row("audio")?;
        let mechanics: u32 = row.parse()?;
        let output: u32 = row.parse()?;
        let samples: u64 = row.parse()?;
        let full_scale_pa = row.scalar()?;
        row.finish()?;
        if output == 0 || output > 192000 || mechanics == 0 || mechanics % output != 0
            || !(1..=16).contains(&(mechanics / output)) || samples == 0
            || samples > 28_800_000 || samples > 600 * u64::from(output) || full_scale_pa <= 0.0 {
            return Err(bad(rows.line, "audio needs positive clocks, an integer mechanics/output ratio 1..=16, positive full scale and at most 600 seconds/28800000 output frames"));
        }
        let mechanical_samples = samples.checked_mul(u64::from(mechanics / output))
            .ok_or_else(|| bad(rows.line, "mechanical sample window overflow"))?;
        let dt = 1.0 / f64::from(mechanics);
        let propagator = read_propagator(&mut rows)?;
        let mut gases = Vec::new();
        let mut row = rows.row("ambient")?;
        let ambient = gas::read_with_sources(&mut row, None, gas_sources, &mut gases)?;
        row.finish()?;
        let mut row = rows.row("exterior")?;
        // A full GAS_RECORD (T P RH or a mixture), as documented; passing the
        // ambient pressure here would switch the legacy record to `T RH`.
        let exterior = gas::read_with_sources(&mut row, None, gas_sources, &mut gases)?;
        if exterior.pressure.to_bits() != ambient.pressure.to_bits() {
            return Err(bad(row.line, "exterior and ambient must share one static pressure; mean interface flow is not modeled"));
        }
        row.finish()?;
        rows.row("upstream")?.finish()?;
        let up = DuctInput::read_with_gas_sources(&mut rows, &ambient, gate, gas_sources, &mut gases)?
            .into_two_sided(&ambient, dt, gate)?;
        rows.row("downstream")?.finish()?;
        let down = DuctInput::read_with_gas_sources(&mut rows, &ambient, gate, gas_sources, &mut gases)?
            .into_two_sided(&ambient, dt, gate)?;
        let sides = [up, down];
        let (inlet, mut impedance) = combined_inlet_state(&sides).map_err(PlateValveInputError::Physics)?;

        let mut row = rows.row("propagation")?;
        let policies = match row.word()? {
            "integer" => None,
            "fractional" => Some([FractionalTubePolicy {
                maximum_frequency_hz: row.scalar()?,
                maximum_summed_phase_error_rad: row.scalar()?,
            }; 2]),
            _ => return Err(bad(row.line, "propagation must be integer or fractional BAND_HZ PHASE_ERROR_RAD")),
        };
        row.finish()?;
        let mut row = rows.row("supply")?;
        let channel = match row.word()? {
            "none" => None,
            "slit" => Some(SlitChannel::from_gas(&inlet, SlitChannelSpec {
                length_m: row.scalar()?, width_m: row.scalar()?, height_m: row.scalar()?,
                transverse_cells: row.count(4096)?, sample_period_s: dt, exit: SlitExit::StaticPressure,
            }).map_err(PlateValveInputError::Physics)?),
            _ => return Err(bad(row.line, "supply must be none or slit LENGTH WIDTH HEIGHT CELLS")),
        };
        row.finish()?;
        if let Some(channel) = &channel {
            impedance = combined_inlet_state_with_supply(&sides, channel)
                .map_err(PlateValveInputError::Physics)?.1;
        }
        let mut row = rows.row("slit_length")?;
        let slit = optional_positive(&mut row)?;
        row.finish()?;
        let mut row = rows.row("receiver")?;
        let receiver = read_receiver(&mut row)?;
        row.finish()?;
        for side in &sides {
            if side.viscothermal.iter().any(|s| s.maximum_frequency_hz < receiver.maximum_frequency_hz) {
                return Err(bad(rows.line, "every selected viscothermal band must cover the receiver band"));
            }
        }
        let count = rows.count("openings", 128)?;
        if count == 0 { return Err(bad(rows.line, "at least one physical radiating opening is required")); }
        let mut openings = Vec::with_capacity(count);
        for _ in 0..count {
            let mut row = rows.row("opening")?;
            let network = match row.word()? {
                "upstream" => 0, "downstream" => 1,
                _ => return Err(bad(row.line, "opening must name upstream or downstream")),
            };
            openings.push(DuctOpening { network, node: row.parse()?,
                center_m: [row.scalar()?, row.scalar()?] });
            row.finish()?;
        }
        let radiation_spec = DuctRadiationSpec { openings, exterior: Some(exterior), receiver,
            reconstruction: DuctReconstruction::Differentiating {
                mechanics_rate_hz: mechanics, output_rate_hz: output,
            } };
        let mut row = rows.row("valve")?;
        let rest_opening = row.scalar()?;
        row.finish()?;
        let mut row = rows.row("limits")?;
        let limits = PairedApertureLimits {
            maximum_abs_displacement_m: row.scalar()?,
            maximum_abs_velocity_m_s: row.scalar()?,
            maximum_energy_j: row.scalar()?,
            max_evaluations: row.count(65536)?,
            velocity_tolerance_m_s: row.scalar()?,
            relative_tolerance: row.scalar()?,
        };
        row.finish()?;
        let mut row = rows.row("coupling")?;
        let coupling = TwoSidedLimits {
            maximum_energy_j: row.scalar()?,
            energy_absolute_tolerance_j: row.scalar()?,
            energy_relative_tolerance: row.scalar()?,
            flow_absolute_tolerance_m3_s: row.scalar()?,
            flow_relative_tolerance: row.scalar()?,
        };
        row.finish()?;
        let mut row = rows.row("contact")?;
        let stiffness = row.scalar()?;
        let alpha = row.scalar()?;
        let contact_loss = row.scalar()?;
        row.finish()?;
        let contact = Obstacle::new(vec![-1.0], 1, 1, vec![0.0], vec![1.0],
            stiffness, alpha, "material-modal file: declared effective uniform slit contact".into())
            .and_then(|c| c.with_internal_loss(contact_loss))
            .map_err(|e| PlateValveInputError::Physics(
                crate::acoustic_realize::AcousticRealizeError::Nonlinear(e.to_string())))?;
        let loss = LossSelection::read(&mut rows)?;
        let blades = rows.count("blades", 2)?;
        if blades == 0 { return Err(bad(rows.line, "one or two independently supplied blades are required")); }
        let mut cache = MaterialCache::default();
        let mut plates = Vec::with_capacity(blades);
        let mut states = Vec::<BladeState>::new();
        let mut materials = Vec::new();
        let (mut nodes, mut triangles) = (0, 0);
        for blade in 0..blades {
            let (source, options, initial, records) = plate::read(
                &mut rows, blade, states.len(), &loss, &mut cache, &mut load_material, gate)?;
            nodes += source.chart().mesh.nodes.len();
            triangles += source.chart().mesh.tris.len();
            states.extend(initial);
            materials.extend(records);
            plates.push((source, options));
        }
        let modes = states.len();
        let mut row = rows.row("compile_limits")?;
        let max_compile_work: u64 = row.parse()?;
        let max_controls = row.count(MAX_CONTROLS)?;
        row.finish()?;
        if max_compile_work > 16_777_216 {
            return Err(bad(rows.line, "pressure compilation exceeds 16777216 visits"));
        }
        if rows.lines.next().is_some() {
            return Err(bad(rows.line + 1, "unexpected trailing material-modal record"));
        }
        let schedule = decode_schedule(gestures)?;
        let mut valve = match loss {
            LossSelection::Elastic => ModalAperture::from_material_plates_with_propagator(
                plates, rest_opening, inlet, impedance, dt, mechanical_samples,
                limits, states, contact, propagator, gate),
            LossSelection::Tensor { properties, limits: projection } => {
                let mut viscous = Vec::with_capacity(plates.len());
                for (source, options) in plates {
                    let law = crate::plate_viscosity::PlateViscosity::from_chart(
                        &source, properties.each_ref().map(String::as_str), gate)
                        .map_err(PlateValveInputError::Physics)?;
                    viscous.push((source, options, law, projection));
                }
                ModalAperture::from_viscous_material_plates_with_propagator(
                    viscous, rest_opening, inlet, impedance, dt, mechanical_samples,
                    limits, states, contact, propagator, gate)
            }
        }.map_err(PlateValveInputError::Physics)?;
        if let Some(length) = slit {
            valve.set_material_viscous_slit(length).map_err(PlateValveInputError::Physics)?;
        }
        let (source, radiation) = match channel {
            Some(channel) => TwoSidedPressure::from_baffled_environment_with_supply(
                valve.into(), sides, [None, None], channel, coupling, policies,
                radiation_spec, max_block, gate),
            None => TwoSidedPressure::from_baffled_environment(
                valve.into(), sides, [None, None], coupling, policies,
                radiation_spec, max_block, gate),
        }.map_err(|e| PlateValveInputError::Gesture(GestureCompileError::Render(e)))?;
        let renderer = ScheduledTwoSidedPressure::new(source, schedule, AperturePerformanceConfig {
            sample_rate_hz: output, samples, max_block, max_compile_work, max_controls,
        }).map_err(PlateValveInputError::Gesture)?;
        checkpoint(gate)?;
        Ok(Self { info: MaterialModalValveInfo {
            input_hash: hash_domain("org.frankensim.fs-couple.material-modal-valve.v1", bytes),
            mechanics_rate_hz: mechanics, propagator, output_rate_hz: output, samples, full_scale_pa,
            blades, modes, nodes, triangles,
        }, renderer, materials, gases, ambient, radiation })
    }

    /// Immutable original source dimensions, coordinates, clocks and PCM scale.
    pub const fn info(&self) -> MaterialModalValveInfo { self.info }
    /// Original mode/duct/radiation owners through the complete scheduled source.
    pub const fn renderer(&self) -> &ScheduledTwoSidedPressure { &self.renderer }
    /// Exact normalized material artifacts and regional declaration addresses.
    pub fn material_records(&self) -> &[ModalValveMaterialRecord] { &self.materials }
    /// Original ambient, regional and exterior gas declarations.
    pub fn gas_records(&self) -> &[GasRecord] { &self.gases }
    /// Default interior state; actual inlet/remote states stay in the two ducts.
    pub const fn ambient_gas(&self) -> &GasState { &self.ambient }
    /// Actual compact reaction models, in the authored physical-opening order.
    pub fn radiation_bindings(&self) -> &[BaffledOutletBinding] { &self.radiation }
    /// Move the complete loaded mechanics, source clock and microphone history.
    pub fn into_renderer(self) -> ScheduledTwoSidedPressure { self.renderer }
}

fn read_propagator(rows: &mut Reader<'_>) -> Result<BladePropagator, PlateValveInputError> {
    if rows.lines.clone().next().and_then(|line| line.split_ascii_whitespace().next())
        != Some("propagator") {
        return Ok(BladePropagator::Midpoint);
    }
    let mut row = rows.row("propagator")?;
    let selected = match row.word()? {
        "midpoint" => BladePropagator::Midpoint,
        "exact-zoh" => BladePropagator::ExactZoh,
        _ => return Err(bad(row.line, "propagator must be midpoint or exact-zoh")),
    };
    row.finish()?;
    Ok(selected)
}

fn optional_positive(row: &mut super::Row<'_>) -> Result<Option<f64>, PlateValveInputError> {
    let word = row.word()?;
    if word == "none" { return Ok(None); }
    let value: f64 = word.parse().map_err(|_| bad(row.line, "expected none or positive length in metres"))?;
    if !value.is_finite() || value <= 0.0 { return Err(bad(row.line, "slit length must be positive finite metres")); }
    Ok(Some(value))
}

#[cfg(test)]
mod tests;
