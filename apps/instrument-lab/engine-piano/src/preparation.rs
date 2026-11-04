//! Cold physical scale preparation for the existing live piano, before its
//! reference Bank builds string modes, static bridge loading and hammer data.
//! Reuses the offline CSV owner and multi-partial tension solver verbatim.
//! No rendered-signal pitch shift, measured tuning curve or new string law.
use crate::{Description, refmods::{geometry, stretch_tuning}};
pub use stretch_tuning::partials::{OctaveReport, OctaveWeights, Profile};

pub mod materials;

/// Bound the text ingress before parsing. The canonical 88-course SI scale is
/// much smaller; this also leaves room for authored source comments.
pub const MAX_SCALE_BYTES: usize = 256 * 1024;

/// Resolve the same bounded atmosphere declaration used by file-driven physics
/// and the live drum: `T P RH` or `mixture T P mole|mass N SPECIES FRACTION…`.
/// This cold ingress delegates every physical coefficient and domain check to
/// fs-material. Canonical thermochemical sources enter through
/// [`resolve_atmosphere_sources`]; trailing text always refuses.
pub fn resolve_atmosphere(text: &str) -> Result<fs_material::gas::GasState, String> {
    resolve_atmosphere_record(text).map(|record| *record.state())
}

pub(crate) fn resolve_atmosphere_record(text: &str)
    -> Result<fs_material::gas_species::input::GasRecord, String> {
    if text.len() > 8192 { return Err("gas declaration exceeds 8 KiB".into()); }
    let mut fields = text.split_ascii_whitespace();
    let gas = fs_material::gas_species::input::read_gas_record(&mut fields, None, None)
        .map_err(|e| e.to_string())?;
    if fields.next().is_some() { return Err("unexpected fields after the gas declaration".into()); }
    Ok(gas)
}

/// Resolve an arbitrary gas from the same canonical source packet as the live
/// drum and kit. Consecutive species/NASA-9/transport pack triples declare each
/// species explicitly; no chemical name selects thermodynamics or transport.
/// The shared reader owns byte limits, composition, source validity and phase
/// admission. The returned record retains all original models and receipts.
///
/// `water_species` is empty only when no water identity is declared. `phase`
/// must be `frozen` or `liquid-water`, including its original saturation checks.
/// Sources may only accompany a complete `thermochemical ...` declaration;
/// they cannot be silently ignored in favor of the built-in air/mixture path.
pub fn resolve_atmosphere_sources(text: &str, bytes: &[u8], lengths: &[u32],
    water_species: &str, phase: &str)
    -> Result<fs_material::gas_species::input::GasRecord, String> {
    use fs_material::gas_species::input::{ThermochemicalPhase, packs::read_source_packet};
    let phase = match phase {
        "frozen" => ThermochemicalPhase::FrozenGas,
        "liquid-water" => ThermochemicalPhase::LiquidWaterLimit,
        _ => return Err("gas source phase must be frozen or liquid-water".into()),
    };
    read_source_packet(text, bytes, lengths,
        (!water_species.is_empty()).then_some(water_species), phase)
        .map_err(|e| e.to_string())
}

/// One explicit preparation policy. The shared solver changes only tension;
/// C4..B4 stay anchored and missing octave-neighbour geometry is never invented.
#[derive(Clone, Copy, Debug)]
pub struct Tuning {
    concert_a_hz: f64,
    profile: Profile,
}
impl Tuning {
    pub fn new(concert_a_hz: f64, bass: [f64; 4], treble: [f64; 4],
        max_partial_hz: Option<f64>) -> Result<Self, String> {
        if !concert_a_hz.is_finite() || !(430.0..=450.0).contains(&concert_a_hz) {
            return Err("concert A must be finite and in 430..450 Hz".into());
        }
        if max_partial_hz.is_some_and(|hz| !hz.is_finite() || hz <= 0.0) {
            return Err("tuning band must be finite and positive".into());
        }
        Ok(Self { concert_a_hz, profile: Profile {
            bass: OctaveWeights::new(bass)?, treble: OctaveWeights::new(treble)?,
            max_partial_hz,
        } })
    }

    pub fn concert_a_hz(self) -> f64 { self.concert_a_hz }
    pub fn profile(self) -> Profile { self.profile }

    /// Fixed WASM packet: A4, four bass weights, four treble weights, band Hz.
    /// An empty packet preserves tensions; band zero alone means no ceiling.
    /// The shared weight owner admits and normalizes the four partial pairs.
    pub fn from_values(values: &[f64]) -> Result<Option<Self>, String> {
        if values.is_empty() { return Ok(None); }
        if values.len() != 10 {
            return Err("tuning needs A4, four bass weights, four treble weights and a band".into());
        }
        Self::new(values[0], [values[1], values[2], values[3], values[4]],
            [values[5], values[6], values[7], values[8]],
            if values[9] == 0.0 { None } else { Some(values[9]) }).map(Some)
    }
}

#[derive(Clone, Debug)]
pub struct PreparationReport {
    pub imported_scale: bool,
    pub courses: usize,
    pub tuning: Option<Tuning>,
    /// Actual signed partial-pair mismatches from the same physical solver.
    pub octaves: Vec<OctaveReport>,
}
impl PreparationReport {
    /// Bounded numeric metadata for native and browser callers. All numbers
    /// originate in admitted controls or the solver's finite physical report.
    /// RMS is a local frequency-mismatch objective, not perceptual roughness.
    pub fn json(&self) -> String {
        let tuning = match self.tuning {
            None => "null".to_string(),
            Some(t) => {
                let band = t.profile.max_partial_hz.map_or_else(|| "null".into(), |x| x.to_string());
                let rms = self.octaves.iter().map(|o| o.weighted_rms_hz).fold(0.0_f64, f64::max);
                format!(concat!("{{\"policy\":\"pinned-string-octaves\",\"concertAHz\":{},",
                    "\"bassOctaveWeights\":{:?},\"trebleOctaveWeights\":{:?},",
                    "\"maxPartialHz\":{},\"octavePairs\":{},\"maxWeightedRmsHz\":{}}}"),
                    t.concert_a_hz, t.profile.bass.values(), t.profile.treble.values(),
                    band, self.octaves.len(), rms)
            }
        };
        format!("{{\"importedScale\":{},\"courses\":{},\"tuning\":{}}}",
            self.imported_scale, self.courses, tuning)
    }
}

/// Import an optional physical scale, then apply an optional explicit tuning.
/// Publish the courses only after BOTH operations succeed. Board coordinates,
/// full damping, bridge rows, acoustic surface, and source labels are untouched.
///
/// With neither option, preserve even the original course allocation. With an
/// imported scale but no tuning, retain its supplied tensions bit for bit. This
/// is the playback route for CSVs already prepared by `piano_tune` or an author.
/// Sparse scales obey the original reader and solver's explicit-neighbour rules.
/// Existing running engines are not modified: construct a new LivePiano or
/// PianoAudio from this prepared Description so all physical loads are rebuilt.
pub fn prepare(description: &mut Description, scale_csv: Option<&str>, tuning: Option<Tuning>)
    -> Result<PreparationReport, String> {
    prepare_with_materials(description, scale_csv, tuning, None)
}

/// Resolve the complete physical string scale before optional tension tuning.
/// Source sections, imported geometry and tuning publish together. A failed
/// material binding or later tuning proposal leaves the description untouched.
pub fn prepare_with_materials(description: &mut Description, scale_csv: Option<&str>,
    tuning: Option<Tuning>, materials: Option<&materials::StringMaterials>)
    -> Result<PreparationReport, String> {
    let mut report = PreparationReport { imported_scale: scale_csv.is_some(),
        courses: description.courses.len(), tuning, octaves: Vec::new() };
    if scale_csv.is_none() && tuning.is_none() && materials.is_none() { return Ok(report); }
    let mut courses = match scale_csv {
        Some(text) => {
            if text.len() > MAX_SCALE_BYTES {
                return Err("physical scale CSV exceeds the 256 KiB input budget".into());
            }
            geometry::read_scale(text)?
        }
        None => description.courses.clone(),
    };
    if let Some(materials) = materials { courses = materials.apply(&courses)?; }
    if let Some(t) = tuning {
        report.octaves = stretch_tuning::partials::tune(&mut courses, t.concert_a_hz, t.profile)?;
        // Retuning changes tension: the material strings face their tensile check again.
        if let Some(materials) = materials { materials.verify_courses(&courses)?; }
    }
    report.courses = courses.len();
    description.courses = courses;
    Ok(report)
}
