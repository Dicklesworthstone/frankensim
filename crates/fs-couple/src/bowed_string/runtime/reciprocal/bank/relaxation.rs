//! Per-layer bending spectra in the existing coupled string/body operator.
//!
//! For the original mass-normalized sine coordinate, integration of curvature
//! gives c_n^2 = DeltaE*I*(n*pi/L)^4/mu. Each selected layer/arm/plane/mode
//! retains its own first-order history. Tension is NOT relaxed, and neither
//! material identities nor relaxation times are averaged across bonded layers.
//! The common Maxwell owner supplies storage, restoring load and dissipation.
use super::*;
use crate::acoustic_realize::AcousticRealizeError;
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::section::relaxation::LayerBendingRelaxationProperties;
use fs_material::visco::GeneralizedMaxwell;
use fs_math::det;
use fs_phs::RelaxationBranch;

// Original specimens can already select an exclusive constitutive spectrum.
// They join the same owner before optional regional and layered additions.
mod source;

/// Explicit bending-modulus increment [Pa] / relaxation-time [s] selectors.
/// An unlisted string receives no new constitutive mechanism.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StringBendingRelaxationSelection {
    pub string: usize,
    pub properties: Vec<(String, String)>,
}

/// Independent source roles for every spatial layer of one original string.
/// A coating and core may use different arm counts and material property names.
#[derive(Clone, Debug, PartialEq)]
pub struct LayeredStringBendingRelaxationSelection {
    pub string: usize,
    /// Core first, then sleeves in the original section order. Empty branch
    /// lists explicitly select elastic layers; omitted layers are an error.
    pub layers: Vec<LayerBendingRelaxationProperties>,
}

/// One source-projected internal strain, in original component/modal coordinates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StringBendingMemoryChannel {
    pub string: usize,
    /// Includes the second transverse plane's offset, never a torsional mode.
    pub mode: usize,
    /// c_n multiplying mass-normalized displacement, in 1/s.
    pub coupling_per_s: f64,
    pub relaxation_time_s: f64,
}

/// Immutable source selection and its exact uniform-span curvature projection.
#[derive(Clone, Debug)]
pub struct StringBankBendingRelaxation {
    selections: Vec<LayeredStringBendingRelaxationSelection>,
    sections: Vec<ConcentricElasticSection>,
    spectra: Vec<Vec<GeneralizedMaxwell>>,
    bands: Vec<(f64, f64)>,
    channels: Vec<StringBendingMemoryChannel>,
    basis: Vec<(usize, Vec<ModalAcousticMode>)>,
    max_memory_states: usize,
}
impl StringBankBendingRelaxation {
    /// Normalized source roles, including explicitly elastic layers.
    #[must_use] pub fn layered_selections(&self) -> &[LayeredStringBendingRelaxationSelection] { &self.selections }
    /// Original material bundles and section integrals, in selection order.
    #[must_use] pub fn sections(&self) -> &[ConcentricElasticSection] { &self.sections }
    /// One original spectrum per layer, nested in selection order.
    #[must_use] pub fn spectra(&self) -> &[Vec<GeneralizedMaxwell>] { &self.spectra }
    /// Order: selection, layer, positive arm, transverse plane, sine mode.
    #[must_use] pub fn channels(&self) -> &[StringBendingMemoryChannel] { &self.channels }
    #[must_use] pub const fn max_memory_states(&self) -> usize { self.max_memory_states }

    fn bind(bank: &ReciprocalStringBank, selections: Vec<LayeredStringBendingRelaxationSelection>,
        max_memory_states: usize, gate: &CancelGate) -> Result<Self, ReciprocalStringError>
    {
        if selections.is_empty() || selections.len() > bank.cards.len()
            || !(1..=256).contains(&max_memory_states) {
            return Err(invalid("string bending memory needs bounded explicit selections and 1..=256 states"));
        }
        let mut sections = Vec::new();
        let mut spectra = Vec::new();
        let mut bands = Vec::new();
        let mut channels = Vec::new();
        let mut basis = Vec::new();
        for (index, selection) in selections.iter().enumerate() {
            poll(gate)?;
            let member = bank.source.strings.get(selection.string)
                .ok_or_else(|| invalid("bending memory names an absent original string"))?;
            if selections[..index].iter().any(|s| s.string == selection.string) {
                return Err(invalid("each original string may select bending memory only once"));
            }
            if member.source.bending_relaxation().is_some() {
                return Err(invalid("an original Prony string cannot select a second bending-memory owner"));
            }
            let section = section(&member.source)?.with_prony_bending(&selection.layers)
                .map_err(owner)?;
            let card = &bank.cards[selection.string];
            let planes: usize = if matches!(&member.source, MaterialStringSource::WithPolarizations(_)) { 2 } else { 1 };
            let component = &bank.system.components()[selection.string];
            if component.modes().len() != member.source.modes()
                || planes * card.mode_count > component.modes().len() {
                return Err(invalid("bending memory lost the complete original string basis"));
            }
            basis.push((selection.string, component.modes().to_vec()));
            for branch in section.branches() {
                poll(gate)?;
                if branch.relaxing_stiffness_n_m2 == 0.0 { continue; }
                if planes.checked_mul(card.mode_count)
                    .and_then(|n| n.checked_add(channels.len())).is_none_or(|n| n > max_memory_states) {
                    return Err(invalid("complete per-layer bending memory exceeds its declared state budget"));
                }
                let delta_ei = branch.relaxing_stiffness_n_m2;
                let scale = det::sqrt(delta_ei / card.linear_density_kg_m);
                if !delta_ei.is_finite() || delta_ei <= 0.0 || !scale.is_finite() || scale <= 0.0 {
                    return Err(invalid("string bending memory stiffness is unrepresentable"));
                }
                for plane in 0..planes { for mode in 0..card.mode_count {
                    poll(gate)?;
                    let k = (mode+1) as f64 * core::f64::consts::PI / card.length_m;
                    let coupling = scale * k * k;
                    if !coupling.is_finite() || coupling <= 0.0 {
                        return Err(invalid("string curvature projection lost a positive material channel"));
                    }
                    channels.push(StringBendingMemoryChannel { string: selection.string,
                        mode: plane*card.mode_count+mode, coupling_per_s: coupling,
                        relaxation_time_s: branch.relaxation_time_s });
                }}
            }
            let mut band = section.omega_band_rad_s();
            if let Some(winding_band) = winding_band(&member.source) {
                band = (band.0.max(winding_band.0), band.1.min(winding_band.1));
            }
            bands.push(band);
            sections.push(section.elastic().clone());
            spectra.push(section.spectra().to_vec());
        }
        Ok(Self { selections, sections, spectra, bands, channels, basis, max_memory_states })
    }

    pub(super) fn validate_band(&self, high: f64, dt: f64) -> Result<(), ReciprocalStringError> {
        if !high.is_finite() || high < 0.0 || !dt.is_finite() || dt <= 0.0
            || self.bands.iter().any(|&(lo,hi)| lo > 0.0 || high > hi)
            || self.channels.iter().any(|c| dt/c.relaxation_time_s > 2.0) {
            return Err(invalid("complete coupled band or clock exceeds original string bending-memory support"));
        }
        Ok(())
    }

    // New body/axial coordinates can extend the embedding, but changing the
    // string mass or basis needs a new cold source projection, never stale rows.
    pub(super) fn branches(&self, system: &FrictionModalSystem, gate: &CancelGate)
        -> Result<Vec<RelaxationBranch>, ReciprocalStringError>
    {
        poll(gate)?;
        for (index, modes) in &self.basis {
            if system.components().get(*index).is_none_or(|m| m.modes() != modes.as_slice()) {
                return Err(invalid("select string gas loading and basis before binding bending memory"));
            }
        }
        let dimension = system.mode_count().checked_mul(2)
            .ok_or_else(|| invalid("string material-memory embedding overflow"))?;
        self.channels.iter().map(|channel| {
            poll(gate)?;
            let start: usize = system.components()[..channel.string].iter().map(|m| m.modes().len()).sum();
            let mut projection = vec![0.0; dimension];
            projection[2*(start+channel.mode)] = channel.coupling_per_s;
            Ok(RelaxationBranch { projection, stiffness: 1.0,
                relaxation_time_s: channel.relaxation_time_s })
        }).collect()
    }
}

fn owner(error: impl core::fmt::Display) -> ReciprocalStringError {
    ReciprocalStringError::Material(AcousticRealizeError::Nonlinear(error.to_string()))
}
fn section(source: &MaterialStringSource) -> Result<ConcentricElasticSection, ReciprocalStringError> {
    match source {
        MaterialStringSource::WithPolarizations(p) => section(p.source()),
        MaterialStringSource::WithTorsion(t) => section(t.bending()),
        MaterialStringSource::Homogeneous(s) => ConcentricElasticSection::try_new(0.0,
            vec![ConcentricMaterialLayer { outer_radius_m: s.radius_m(), material: s.material().clone() }])
            .map_err(owner),
        MaterialStringSource::Concentric { section, .. } => Ok(section.elastic().clone()),
        MaterialStringSource::ConcentricProny { section, .. } => Ok(section.elastic().clone()),
        MaterialStringSource::Wound(s) => Ok(s.section().core().elastic().clone()),
    }
}
// Raising instantaneous core stiffness must not escape the wrapping material's
// support, even though the wrap contributes no bending-memory channel.
fn winding_band(source: &MaterialStringSource) -> Option<(f64, f64)> {
    match source {
        MaterialStringSource::WithPolarizations(p) => winding_band(p.source()),
        MaterialStringSource::WithTorsion(t) => winding_band(t.bending()),
        MaterialStringSource::Wound(s) => Some(s.section().omega_band_rad_s()),
        _ => None,
    }
}
fn validate_properties(properties: &[(String,String)]) -> Result<(), ReciprocalStringError> {
    if properties.is_empty() || properties.len()>8 {
        return Err(invalid("string memory requires 1..=8 explicit modulus/time pairs"));
    }
    let mut pairs = std::collections::BTreeSet::new();
    let mut moduli = std::collections::BTreeSet::new();
    let mut times = std::collections::BTreeSet::new();
    for (modulus,time) in properties {
        if modulus.is_empty() || time.is_empty() || modulus==time
            || !pairs.insert((modulus.as_str(),time.as_str()))
            || times.contains(modulus.as_str()) || moduli.contains(time.as_str()) {
            return Err(invalid("string memory needs distinct source pairs with consistent quantities"));
        }
        moduli.insert(modulus.as_str()); times.insert(time.as_str());
    }
    Ok(())
}

impl MaterialStringBankSource {
    /// Source sections and spectra survive multibow, forcing and PCM ownership.
    #[must_use]
    pub const fn string_bending_relaxation(&self) -> Option<&StringBankBendingRelaxation> {
        self.string_bending_relaxation.as_ref()
    }
}
impl ReciprocalStringBank {
    /// Add source-selected bending memory to homogeneous, bonded or wound strings.
    /// For the sliding close-wound reduction only the original core bends; its
    /// arms use the total core-plus-wrap mass and retain wrapping-band limits.
    /// Every material layer supplies its own modulus increments and time scales;
    /// both transverse planes share those data, not a fitted per-mode Q.
    /// The unchanged exact owner couples this storage to body, bridges and bows.
    /// Existing explicitly sourced Kelvin-Voigt bending viscosity is a parallel
    /// instantaneous dashpot; the selected spectrum must describe ADDITIONAL
    /// relaxing bending modulus, with the original Young modulus as equilibrium.
    ///
    /// Select body memory and gas loading first. Additional torque connections
    /// and axial coordinates retain this history through the shared rebuilding
    /// paths. Prepared mechanical-only propagation refuses rather than losing z.
    /// `initial_memory` follows channels() order: None is relaxed; an explicit
    /// zero vector is unrelaxed. Existing body history is never reset.
    /// Original exclusive Prony specimens already carry their histories and
    /// cannot be selected again. Selections on other strings compose after
    /// those original histories and any regional body history.
    ///
    /// Frozen uniform-span small-curvature bending only. No axial/torsional
    /// memory, finite-strain viscoelasticity, winding slip, heat evolution or
    /// measured instrument-fidelity claim is supplied by this selection.
    pub fn with_material_string_bending_relaxation(
        self, selections: Vec<StringBendingRelaxationSelection>,
        initial_memory: Option<Vec<f64>>, max_memory_states: usize, gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        poll(gate)?;
        if selections.is_empty() || selections.len() > self.cards.len()
            || !(1..=256).contains(&max_memory_states) {
            return Err(invalid("string bending memory needs bounded explicit selections and 1..=256 states"));
        }
        let mut layered = Vec::with_capacity(selections.len());
        for selection in selections {
            poll(gate)?;
            validate_properties(&selection.properties)?;
            let member = self.source.strings.get(selection.string)
                .ok_or_else(|| invalid("bending memory names an absent original string"))?;
            let layer_count = section(&member.source)?.layers().len();
            layered.push(LayeredStringBendingRelaxationSelection {
                string: selection.string,
                layers: vec![LayerBendingRelaxationProperties {
                    equilibrium_modulus: "young_modulus".into(), branches: selection.properties,
                }; layer_count],
            });
        }
        self.with_material_layered_string_bending_relaxation(layered, initial_memory,
            max_memory_states, gate)
    }

    /// Bind a separate source spectrum for each core/coating layer using the
    /// shared material-section owner. Layer-specific names and arm counts need
    /// no new string equation. Explicit empty spectra are elastic, but missing
    /// layer bindings or required property support refuse atomically.
    ///
    /// The same cold-state, memory-budget, complete-band and exact-runtime
    /// requirements as `with_material_string_bending_relaxation` apply. Every
    /// layer explicitly identifies its equilibrium modulus; a sampled storage
    /// modulus cannot silently replace it. No new thermal or axial memory is added.
    pub fn with_material_layered_string_bending_relaxation(
        mut self, selections: Vec<LayeredStringBendingRelaxationSelection>,
        initial_memory: Option<Vec<f64>>, max_memory_states: usize, gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        poll(gate)?;
        let body_count = self.source.material_chart_relaxation.as_ref()
            .map_or(0, |image| image.coupling_rows_per_s().len())
            +self.source.material_chart_tensor_relaxation.as_ref()
                .map_or(0,|image|image.coupling_rows_per_s().len());
        let existing_count = self.source.source_memory_count(&self.cards, gate)?
            .checked_add(body_count).ok_or_else(|| invalid("complete material history extent overflow"))?;
        if self.system.samples_rendered()!=0 || self.source.string_bending_relaxation.is_some()
            || self.system.material_relaxation_branches().len()!=existing_count {
            return Err(invalid("string bending memory requires one cold binding after declared body memory"));
        }
        let image = StringBankBendingRelaxation::bind(&self,selections,max_memory_states,gate)?;
        image.validate_band(self.system.linear_frequency_upper_bound_rad_s(),
            self.system.sample_period_s())?;
        let branches = image.branches(&self.system,gate)?;
        self.system = self.system.with_additional_material_relaxation(branches,initial_memory,gate)
            .map_err(ReciprocalStringError::Friction)?;
        self.source.string_bending_relaxation=Some(image);
        self.source.validate_body_memory_frequencies(&self.system,gate)?;
        self.source.validate_material_bridge_frequencies(&self.system,gate)?;
        poll(gate)?;
        Ok(self)
    }
}
