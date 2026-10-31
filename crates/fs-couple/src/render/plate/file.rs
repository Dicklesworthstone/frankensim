//! Flat plate geometry/material/force files lowered through the existing owners.
//!
//! Each triangle selects an explicit rotated orthotropic section. Shared nodes
//! mean perfect bonding, not layer offsets, joints or delamination. Supports and
//! the unit-total-force nodal footprint are explicit input. Eigenfrequencies,
//! masses, force participation and signed radiation areas come from the chart
//! reduction, never from an authored oscillator list. The admitted search window
//! need not include every physical mode; no continuum convergence is implied.
//! The exterior medium is one shared gas record: its density and sound speed
//! reach the observer, and a record that does not resolve refuses (never air).
//! The optional `radiation-load none|baffled-piston` record after `observer`
//! selects the existing compact reaction model from the SAME gas and clock.
//! Omission keeps the original unloaded image; no loss model is inferred.
//!
//! A section is either bare numbers ending in an explicit `fixture LABEL`
//! (`section THICK RHO E1 E2 NU12 G12 ANGLE fixture LABEL`), or pack-backed:
//! `section material isotropic THICK PACK STATE [AXIS=VALUE ...]` or
//! `section material orthotropic12|orthotropic21 THICK ANGLE PACK STATE
//! [AXIS=VALUE ...]`, resolved through `crate::material_source` from a
//! caller-supplied material-card pack compiled at the named STATE. Bare
//! numbers without a label refuse (bead pexfp.32). `section_provenance` keeps
//! each section's pack, hash, state and point, or its fixture label.

use std::collections::BTreeSet;
use std::str::{FromStr, Lines, SplitAsciiWhitespace};
use fs_blake3::{ContentHash, hash_domain};
use fs_material::MaterialError;
use fs_material::gas::GasState;
use fs_matdb::NormalizedMaterialCardPack;
use crate::material_source::{MaterialRecord, fixture_label, plate_section};
use crate::thin_plate::PlateMaterialModel;
use fs_material::gas_species::input::{GasInputError, GasRecord, ThermochemicalSources, read_gas_record};
use fs_plate::{AssemblyOptions, EdgeSupport, PlateChart, PlateMesh, PlateSection};
use crate::acoustic_realize::AcousticRealizeError;
use crate::render::{ControlDelta, RenderContext, RenderError, RenderVoice};
use crate::render::schedule::{ScheduledControl, ScheduledRenderer};
use crate::thin_plate::{PlateChartRadiation, certified_chart_radiators};
use super::{CompactPlateVoice, PlateVoiceConfig};

/// Versioned grammar for the source geometry, mechanics and force history.
pub const PLATE_PERFORMANCE_SCHEMA: &str = "frankensim-plate-performance-v2";
/// Maximum source bytes read or parsed by this input path.
pub const MAX_PLATE_PERFORMANCE_BYTES: usize = 4 * 1024 * 1024;
const MAX_NODES: usize = 512;
const MAX_TRIANGLES: usize = 2048;
const MAX_SECTIONS: usize = 64;
const MAX_MODES: usize = 64;
const MAX_EVENTS: usize = 16_384;

/// Input/geometry/reduction/runtime refusal; no partial runtime is returned.
#[derive(Debug)]
pub enum PlatePerformanceError {
    /// Syntax, domain or resource admission failed at a source line.
    Input {
        /// One-based source line (zero for whole-file admission).
        line: usize,
        /// Refused condition.
        what: &'static str,
    },
    /// The declared exterior gas was refused by its material owner.
    Gas {
        /// One-based source line of the `ambient` record.
        line: usize,
        /// Original material refusal.
        source: MaterialError,
    },
    /// The existing plate geometry or section owner refused.
    Geometry(fs_plate::PlateError),
    /// The existing plate reduction refused.
    Reduction(AcousticRealizeError),
    /// The scheduled voice refused.
    Render(RenderError),
    /// A section's material record, fixture label or caller-supplied pack refused.
    MaterialSource {
        /// One-based source line of the `section` record.
        line: usize,
        /// Original record, loader or material-owner refusal.
        source: String,
    },
}
impl core::fmt::Display for PlatePerformanceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Input { line, what } => write!(f, "plate input line {line}: {what}"),
            Self::Gas { line, source } => write!(f, "plate input line {line}: ambient gas refused: {source}"),
            Self::Geometry(e) => write!(f, "plate geometry/section: {e}"),
            Self::Reduction(e) => write!(f, "plate reduction: {e}"),
            Self::Render(e) => write!(f, "plate render: {e}"),
            Self::MaterialSource { line, source } => write!(f, "plate input line {line}: material section: {source}"),
        }
    }
}
impl std::error::Error for PlatePerformanceError {}

/// Facts about the admitted source, not a physical-validation certificate.
#[derive(Clone, Copy, Debug)]
pub struct PlatePerformanceInfo {
    /// Exact source bytes hashed in a plate-performance-specific domain.
    pub input_hash: ContentHash,
    /// Explicit audio clock [Hz].
    pub sample_rate_hz: u32,
    /// Half-open render interval [0, samples).
    pub samples: u64,
    /// Declared pascals mapped to PCM full scale.
    pub full_scale_pa: f64,
    /// Resolved exterior gas the observer uses.
    pub gas: GasState,
    /// Observation distance [m].
    pub listener_m: f64,
    /// Time of flight to the listener, r/c [s].
    pub propagation_delay_s: f64,
    /// Whether the input explicitly requested the compact piston reaction.
    pub baffled_piston_loading: bool,
    /// Geometry counts of the actual assembled chart.
    pub nodes: usize,
    /// Triangle count of the actual assembled chart.
    pub triangles: usize,
    /// Authored constitutive section count.
    pub sections: usize,
    /// Requested upper limit, not an assertion of modal completeness.
    pub requested_modes: usize,
    /// Actual positive in-window modes retained by the reduction.
    pub retained_modes: usize,
    /// Authored sample-addressed force assignments.
    pub force_events: usize,
}

/// Complete, admitted geometry-derived plate performance.
pub struct PlatePerformance {
    info: PlatePerformanceInfo,
    gas_record: GasRecord,
    renderer: ScheduledRenderer,
    section_provenance: Vec<String>,
}
impl PlatePerformance {
    /// Read a bounded flat triangular plate and its physical force history.
    /// Every source record is consumed before the eigensolve. The chart/section,
    /// reducer, plate voice and scheduler retain their own admission boundaries.
    /// Caller-supplied meshes must be conforming and nonoverlapping: this parser
    /// does not certify arbitrary triangle-soup topology or source materials.
    /// An explicit `radiation-load baffled-piston` attaches the existing reaction
    /// filter from the same exterior gas. Reduction is offline and currently
    /// noncancellable; its fixed owner solver policy is not a wall-time bound.
    ///
    /// # Errors
    /// Bad syntax/limits, gas, geometry/section, eigenproblem or rendering admission.
    pub fn from_bytes(bytes: &[u8], max_block: usize) -> Result<Self, PlatePerformanceError> {
        Self::from_bytes_with_gas_sources(bytes, max_block, None)
    }

    /// [`PlatePerformance::from_bytes`] with an explicit thermochemical source
    /// registry for `ambient thermochemical` records. Without one, such a record
    /// refuses; supplied sources never select a gas by themselves.
    ///
    /// # Errors
    /// As [`PlatePerformance::from_bytes`], plus missing or refused gas sources.
    pub fn from_bytes_with_gas_sources(
        bytes: &[u8],
        max_block: usize,
        sources: Option<&ThermochemicalSources<'_>>,
    ) -> Result<Self, PlatePerformanceError> {
        Self::from_bytes_with_all_sources(bytes, max_block, sources,
            |_| Err("a pack-backed material section needs an explicit material-card loader".into()))
    }

    /// As [`PlatePerformance::from_bytes_with_gas_sources`], with `load_material`
    /// supplying each material-card pack a `section material` record names (by
    /// its path relative to the input). Each pack is loaded once.
    ///
    /// # Errors
    /// As [`PlatePerformance::from_bytes`], plus a loader refusal, a card compiled
    /// at a state other than the record's, or a missing section property, by name.
    pub fn from_bytes_with_all_sources(
        bytes: &[u8],
        max_block: usize,
        sources: Option<&ThermochemicalSources<'_>>,
        mut load_material: impl FnMut(&str) -> Result<NormalizedMaterialCardPack, String>,
    ) -> Result<Self, PlatePerformanceError> {
        let input = Parsed::read_with_materials(bytes, max_block, sources, &mut load_material)?;
        let bodies = certified_chart_radiators(&input.chart, &input.options)
            .map_err(PlatePerformanceError::Reduction)?;
        let retained_modes = bodies.len();
        let voice = if input.baffled_piston_loading {
            CompactPlateVoice::from_radiators_with_piston_loads(bodies, input.initial_force_n, input.config)
        } else {
            CompactPlateVoice::from_radiators(bodies, input.initial_force_n, input.config)
        }
            .map_err(PlatePerformanceError::Render)?;
        let info = PlatePerformanceInfo {
            input_hash: hash_domain("org.frankensim.fs-couple.plate-performance.v1", bytes),
            sample_rate_hz: input.config.sample_rate_hz, samples: input.samples,
            full_scale_pa: input.full_scale_pa, gas: input.config.gas,
            listener_m: input.config.listener_m,
            propagation_delay_s: input.config.listener_m / input.config.gas.sound_speed,
            baffled_piston_loading: input.baffled_piston_loading,
            nodes: input.chart.mesh.nodes.len(),
            triangles: input.chart.mesh.tris.len(), sections: input.sections,
            requested_modes: input.options.n_modes, retained_modes,
            force_events: input.events.len(),
        };
        let context = RenderContext::new(vec![RenderVoice::CompactPlate(voice)], max_block);
        let renderer = ScheduledRenderer::new(context, input.events, input.max_events)
            .map_err(PlatePerformanceError::Render)?;
        Ok(Self { info, gas_record: input.gas_record, renderer, section_provenance: input.section_provenance })
    }

    /// Each section's source, in section order: `pack PATH HASH state STATE at
    /// AXIS=VALUE,...` for a pack-backed section, `fixture LABEL` for labelled
    /// bare numbers.
    #[must_use]
    pub fn section_provenance(&self) -> &[String] { &self.section_provenance }

    /// Retained source and reduction counts for the output consumer.
    #[must_use]
    pub const fn info(&self) -> PlatePerformanceInfo { self.info }

    /// The exterior gas record, with any thermochemical source receipts.
    #[must_use]
    pub const fn gas_record(&self) -> &GasRecord { &self.gas_record }

    /// Transfer the complete physical state to existing rendering/WAV consumers.
    #[must_use]
    pub fn into_renderer(self) -> ScheduledRenderer { self.renderer }
}

struct Parsed {
    gas_record: GasRecord,
    chart: PlateChart,
    options: PlateChartRadiation,
    config: PlateVoiceConfig,
    baffled_piston_loading: bool,
    samples: u64,
    full_scale_pa: f64,
    initial_force_n: f64,
    sections: usize,
    section_provenance: Vec<String>,
    events: Vec<ScheduledControl>,
    max_events: usize,
}
impl Parsed {
    fn read(
        bytes: &[u8],
        max_block: usize,
        sources: Option<&ThermochemicalSources<'_>>,
    ) -> Result<Self, PlatePerformanceError> {
        Self::read_with_materials(bytes, max_block, sources,
            &mut |_| Err("a pack-backed material section needs an explicit material-card loader".into()))
    }
    #[allow(clippy::too_many_lines)] // fixed-order grammar; all records precede reduction
    fn read_with_materials(
        bytes: &[u8],
        max_block: usize,
        sources: Option<&ThermochemicalSources<'_>>,
        load_material: &mut dyn FnMut(&str) -> Result<NormalizedMaterialCardPack, String>,
    ) -> Result<Self, PlatePerformanceError> {
        if bytes.len() > MAX_PLATE_PERFORMANCE_BYTES || !(1..=65_536).contains(&max_block) {
            return Err(bad(0, "source exceeds 4 MiB or callback capacity is outside 1..=65536"));
        }
        let text = std::str::from_utf8(bytes).map_err(|_| bad(0, "source must be UTF-8"))?;
        let mut r = Reader { lines: text.lines(), line: 0 };
        r.row(PLATE_PERFORMANCE_SCHEMA)?.finish()?;
        let mut audio = r.row("audio")?;
        let sample_rate_hz: u32 = audio.parse()?;
        let samples: u64 = audio.parse()?;
        let full_scale_pa = audio.scalar()?;
        audio.finish()?;
        if sample_rate_hz == 0 || samples == 0 || samples > u64::from(sample_rate_hz) * 600
            || samples > crate::pcm_wav::stream::MAX_PCM16_WAV_SAMPLES || full_scale_pa <= 0.0
        { return Err(bad(r.line, "audio needs a positive rate, full scale and RIFF-sized duration at most 600 seconds")); }
        let mut ambient = r.row("ambient")?;
        let gas_record = read_gas_record(&mut ambient.fields, None, sources).map_err(|error| match error {
            GasInputError::Input(what) => bad(ambient.line, what),
            GasInputError::Material(source) => PlatePerformanceError::Gas { line: ambient.line, source },
        })?;
        ambient.finish()?;
        let gas = *gas_record.state();
        let listener_m: f64 = r.one("observer")?;
        let baffled_piston_loading = if r.lines.clone().next()
            .and_then(|line| line.split_ascii_whitespace().next()) == Some("radiation-load")
        {
            let mut row = r.row("radiation-load")?;
            let enabled = match row.word()? {
                "none" => false,
                "baffled-piston" => true,
                _ => return Err(bad(row.line, "radiation-load must be none or baffled-piston")),
            };
            row.finish()?;
            enabled
        } else { false };
        let mut limits = r.row("limits")?;
        let max_modes = limits.count(MAX_MODES)?;
        let nyquist_guard_fraction = limits.scalar()?;
        let maximum_abs_force_n = limits.scalar()?;
        let maximum_abs_pressure_pa = limits.scalar()?;
        let max_events = limits.count(MAX_EVENTS)?;
        limits.finish()?;
        let config = PlateVoiceConfig { sample_rate_hz, max_modes, nyquist_guard_fraction,
            gas, listener_m, maximum_abs_force_n, maximum_abs_pressure_pa };
        config.validate().map_err(PlatePerformanceError::Render)?;
        let mut mechanics = r.row("mechanics")?;
        let support = match mechanics.word()? {
            "simply-supported" => EdgeSupport::SimplySupported,
            "clamped" => EdgeSupport::Clamped,
            _ => return Err(bad(mechanics.line, "support must be simply-supported or clamped")),
        };
        let pretension = mechanics.scalar()?;
        let damping_ratio = mechanics.scalar()?;
        let low_hz = mechanics.scalar()?;
        let high_hz = mechanics.scalar()?;
        let n_modes = mechanics.count(max_modes)?;
        mechanics.finish()?;
        if pretension < 0.0 || damping_ratio < 0.0 || low_hz < 0.0 || high_hz <= low_hz
            || high_hz > 0.5 * f64::from(sample_rate_hz) * nyquist_guard_fraction || n_modes == 0
        { return Err(bad(r.line, "nonnegative prestress/damping, ordered sub-Nyquist window and positive mode count required")); }
        let eigenvalue_window = ((core::f64::consts::TAU * low_hz).powi(2),
            (core::f64::consts::TAU * high_hz).powi(2));
        let sections_count = r.count("sections", MAX_SECTIONS)?;
        if sections_count == 0 { return Err(bad(r.line, "at least one section is required")); }
        let mut sections = Vec::with_capacity(sections_count);
        let mut section_provenance = Vec::with_capacity(sections_count);
        let mut packs: Vec<(String, NormalizedMaterialCardPack)> = Vec::new();
        for _ in 0..sections_count {
            let mut row = r.row("section")?;
            let line = row.line;
            let ingress = |e: crate::material_source::MaterialIngressError|
                PlatePerformanceError::MaterialSource { line, source: e.0 };
            if row.fields.clone().next() == Some("material") {
                row.word()?;
                let model = match row.word()? {
                    "isotropic" => PlateMaterialModel::Isotropic,
                    "isotropic-young-shear" => PlateMaterialModel::IsotropicYoungShear,
                    "isotropic-bulk-shear" => PlateMaterialModel::IsotropicBulkShear,
                    "orthotropic12" => PlateMaterialModel::Orthotropic12,
                    "orthotropic21" => PlateMaterialModel::Orthotropic21,
                    _ => return Err(bad(line, "a material section needs isotropic, isotropic-young-shear, isotropic-bulk-shear, orthotropic12 or orthotropic21")),
                };
                let thickness = row.scalar()?;
                let angle = if model.is_isotropic() { 0.0 } else { row.scalar()? };
                let record = MaterialRecord::parse(&row.rest()).map_err(ingress)?;
                if !packs.iter().any(|(path, _)| *path == record.pack) {
                    if packs.len() >= MAX_SECTIONS { return Err(bad(line, "too many material packs")); }
                    let pack = load_material(&record.pack)
                        .map_err(|source| PlatePerformanceError::MaterialSource { line, source })?;
                    packs.push((record.pack.clone(), pack));
                }
                let pack = &packs.iter().find(|(path, _)| *path == record.pack).expect("loaded above").1;
                let card = record.card(pack).map_err(ingress)?;
                let point = record.query_point(card).map_err(ingress)?;
                let (section, _) = plate_section(card, &point, model, thickness, angle).map_err(ingress)?;
                sections.push(section);
                section_provenance.push(record.provenance(pack));
            } else {
                let thickness = row.scalar()?;
                let density = row.scalar()?;
                let e1 = row.scalar()?;
                let e2 = row.scalar()?;
                let nu12 = row.scalar()?;
                let g12 = row.scalar()?;
                let angle = row.scalar()?;
                let label = fixture_label(&row.rest()).map_err(ingress)?;
                sections.push(PlateSection::orthotropic_plane_stress_at_angle(
                    e1, e2, nu12, g12, thickness, density, angle,
                ).map_err(PlatePerformanceError::Geometry)?);
                section_provenance.push(format!("fixture {label}"));
            }
        }
        let nodes_count = r.count("nodes", MAX_NODES)?;
        if nodes_count < 3 { return Err(bad(r.line, "at least three nodes required")); }
        let mut nodes = Vec::with_capacity(nodes_count);
        for _ in 0..nodes_count {
            let mut row = r.row("node")?;
            nodes.push((row.scalar()?, row.scalar()?));
            row.finish()?;
        }
        let triangles_count = r.count("triangles", MAX_TRIANGLES)?;
        if triangles_count == 0 { return Err(bad(r.line, "at least one triangle required")); }
        let mut triangles = Vec::with_capacity(triangles_count);
        let mut assigned = Vec::with_capacity(triangles_count);
        for _ in 0..triangles_count {
            let mut row = r.row("triangle")?;
            let tri = [row.parse()?, row.parse()?, row.parse()?];
            let section: usize = row.parse()?;
            let value = sections.get(section).copied()
                .ok_or_else(|| bad(row.line, "triangle references unknown section"))?;
            row.finish()?;
            triangles.push(tri);
            assigned.push(value);
        }
        let mesh = PlateMesh::from_unstructured(nodes, triangles).map_err(PlatePerformanceError::Geometry)?;
        let support_count = r.count("supports", nodes_count)?;
        let mut supports = Vec::with_capacity(support_count);
        let mut seen = BTreeSet::new();
        for _ in 0..support_count {
            let node: usize = r.one("support")?;
            if node >= nodes_count || !seen.insert(node) { return Err(bad(r.line, "support node must exist and be unique")); }
            supports.push(node);
        }
        let chart = PlateChart::with_boundary_and_regions(mesh, sections[0], supports, Vec::new())
            .and_then(|chart| chart.with_element_sections(assigned))
            .map_err(PlatePerformanceError::Geometry)?;
        let weights_count = r.count("footprint", nodes_count)?;
        if weights_count == 0 { return Err(bad(r.line, "force footprint must be explicit and nonempty")); }
        let mut unit_force_weights = vec![0.0; nodes_count];
        seen.clear();
        let mut weight_sum = 0.0_f64;
        for _ in 0..weights_count {
            let mut row = r.row("weight")?;
            let node: usize = row.parse()?;
            let weight = row.scalar()?;
            if node >= nodes_count || !seen.insert(node) { return Err(bad(row.line, "force node must exist and be unique")); }
            row.finish()?;
            unit_force_weights[node] = weight;
            weight_sum += weight;
        }
        if !weight_sum.is_finite() || (weight_sum - 1.0).abs() > 1e-10 {
            return Err(bad(r.line, "nodal force weights must sum to one; they are never normalized silently"));
        }
        let initial_force_n: f64 = r.one("initial_force_n")?;
        check_force(initial_force_n, maximum_abs_force_n, r.line)?;
        let event_count = r.count("events", max_events)?;
        let mut events = Vec::with_capacity(event_count);
        for _ in 0..event_count {
            let mut row = r.row("force")?;
            let sample: u64 = row.parse()?;
            let force_n = row.scalar()?;
            if sample >= samples { return Err(bad(row.line, "force event lies outside the render interval")); }
            check_force(force_n, maximum_abs_force_n, row.line)?;
            row.finish()?;
            events.push(ScheduledControl { sample, delta: ControlDelta::SetPlateForce { voice: 0, force_n } });
        }
        if r.lines.next().is_some() { return Err(bad(r.line + 1, "unexpected trailing record")); }
        let options = PlateChartRadiation { assembly: AssemblyOptions { pretension, support },
            eigenvalue_window, n_modes, damping_ratio, unit_force_weights };
        Ok(Self { gas_record, chart, options, config, baffled_piston_loading,
            samples, full_scale_pa, initial_force_n,
            sections: sections_count, section_provenance, events, max_events })
    }
}
fn check_force(value: f64, maximum: f64, line: usize) -> Result<(), PlatePerformanceError> {
    if !value.is_finite() || value.abs() > maximum { return Err(bad(line, "force exceeds its finite newton limit")); }
    Ok(())
}
fn bad(line: usize, what: &'static str) -> PlatePerformanceError { PlatePerformanceError::Input { line, what } }
struct Reader<'a> { lines: Lines<'a>, line: usize }
impl<'a> Reader<'a> {
    fn row(&mut self, key: &str) -> Result<Row<'a>, PlatePerformanceError> {
        self.line += 1;
        let line = self.lines.next().ok_or_else(|| bad(self.line, "missing record"))?;
        let mut fields = line.split_ascii_whitespace();
        if fields.next() != Some(key) { return Err(bad(self.line, "unexpected record kind/order")); }
        Ok(Row { fields, line: self.line })
    }
    fn one<T: FromStr>(&mut self, key: &str) -> Result<T, PlatePerformanceError> {
        let mut row = self.row(key)?;
        let value = row.parse()?;
        row.finish()?;
        Ok(value)
    }
    fn count(&mut self, key: &str, maximum: usize) -> Result<usize, PlatePerformanceError> {
        let mut row = self.row(key)?;
        let count = row.count(maximum)?;
        row.finish()?;
        Ok(count)
    }
}
struct Row<'a> { fields: SplitAsciiWhitespace<'a>, line: usize }
impl<'a> Row<'a> {
    fn word(&mut self) -> Result<&'a str, PlatePerformanceError> {
        self.fields.next().ok_or_else(|| bad(self.line, "missing field"))
    }
    fn parse<T: FromStr>(&mut self) -> Result<T, PlatePerformanceError> {
        self.word()?.parse().map_err(|_| bad(self.line, "invalid number"))
    }
    fn scalar(&mut self) -> Result<f64, PlatePerformanceError> {
        let value: f64 = self.parse()?;
        if !value.is_finite() { return Err(bad(self.line, "physical numbers must be finite")); }
        Ok(value)
    }
    fn count(&mut self, maximum: usize) -> Result<usize, PlatePerformanceError> {
        let value = self.parse()?;
        if value > maximum { return Err(bad(self.line, "count exceeds input budget")); }
        Ok(value)
    }
    fn finish(mut self) -> Result<(), PlatePerformanceError> {
        if self.fields.next().is_some() { return Err(bad(self.line, "unexpected extra field")); }
        Ok(())
    }
    /// Every remaining field, consuming the row.
    fn rest(&mut self) -> Vec<&'a str> { self.fields.by_ref().collect() }
}

#[cfg(test)]
mod tests;
