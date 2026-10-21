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
use fs_material::state_point::ResolvedMaterialStatePoint;
use fs_material::visco::GeneralizedMaxwell;
use fs_math::det;
use fs_phs::RelaxationBranch;
use fs_qty::{Density, Dims, Pressure, QuantitySpec};

/// Explicit bending-modulus increment [Pa] / relaxation-time [s] selectors.
/// An unlisted string receives no new constitutive mechanism.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct StringBendingRelaxationSelection {
    pub string: usize,
    pub properties: Vec<(String, String)>,
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
    selections: Vec<StringBendingRelaxationSelection>,
    sections: Vec<ConcentricElasticSection>,
    spectra: Vec<Vec<GeneralizedMaxwell>>,
    bands: Vec<(f64, f64)>,
    channels: Vec<StringBendingMemoryChannel>,
    basis: Vec<(usize, Vec<ModalAcousticMode>)>,
    max_memory_states: usize,
}
impl StringBankBendingRelaxation {
    #[must_use] pub fn selections(&self) -> &[StringBendingRelaxationSelection] { &self.selections }
    /// Original material bundles and section integrals, in selection order.
    #[must_use] pub fn sections(&self) -> &[ConcentricElasticSection] { &self.sections }
    /// One original spectrum per layer, nested in selection order.
    #[must_use] pub fn spectra(&self) -> &[Vec<GeneralizedMaxwell>] { &self.spectra }
    /// Order: selection, layer, positive arm, transverse plane, sine mode.
    #[must_use] pub fn channels(&self) -> &[StringBendingMemoryChannel] { &self.channels }
    #[must_use] pub const fn max_memory_states(&self) -> usize { self.max_memory_states }

    fn bind(bank: &ReciprocalStringBank, selections: Vec<StringBendingRelaxationSelection>,
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
            validate_properties(&selection.properties)?;
            let section = section(&member.source)?;
            let card = &bank.cards[selection.string];
            let planes = if matches!(&member.source, MaterialStringSource::WithPolarizations(_)) { 2 } else { 1 };
            let component = &bank.system.components()[selection.string];
            if component.modes().len() != member.source.modes()
                || planes * card.mode_count > component.modes().len() {
                return Err(invalid("bending memory lost the complete original string basis"));
            }
            basis.push((selection.string, component.modes().to_vec()));
            let mut layer_spectra = Vec::new();
            for (layer, geometry) in section.layers().iter().zip(section.layer_integrals()) {
                poll(gate)?;
                let mut band = (0.0_f64, f64::MAX);
                coefficient(&layer.material, "density", Density::DIMS, &mut band)?;
                let equilibrium = coefficient(&layer.material, "young_modulus", Pressure::DIMS, &mut band)?;
                let terms = selection.properties.iter().map(|(modulus, time)| Ok((
                    coefficient(&layer.material, modulus, Pressure::DIMS, &mut band)?,
                    coefficient(&layer.material, time, Dims([0,0,1,0,0,0]), &mut band)?,
                ))).collect::<Result<Vec<_>, ReciprocalStringError>>()?;
                let spectrum = GeneralizedMaxwell::new(equilibrium, terms).map_err(owner)?;
                if band.0 > band.1 || !band.0.is_finite() || !band.1.is_finite() || band.0 < 0.0 {
                    return Err(invalid("string bending spectrum has no common source frequency domain"));
                }
                let positive = spectrum.terms.iter().filter(|(e,_)| *e > 0.0).count();
                if positive.checked_mul(planes).and_then(|n| n.checked_mul(card.mode_count))
                    .and_then(|n| n.checked_add(channels.len())).is_none_or(|n| n > max_memory_states) {
                    return Err(invalid("complete per-layer bending memory exceeds its declared state budget"));
                }
                for &(modulus, tau) in &spectrum.terms {
                    if modulus == 0.0 { continue; }
                    let delta_ei = modulus * geometry.second_moment_m4();
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
                            relaxation_time_s: tau });
                    }}
                }
                bands.push(band);
                layer_spectra.push(spectrum);
            }
            sections.push(section);
            spectra.push(layer_spectra);
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
fn coefficient(source: &ResolvedMaterialStatePoint, key: &str, dims: Dims,
    band: &mut (f64,f64)) -> Result<f64, ReciprocalStringError>
{
    let p = source.property(key).ok_or_else(|| invalid("string memory lacks an explicit source coefficient"))?;
    if p.requirement().quantity() != QuantitySpec::dimensional(dims)
        || !p.is_constant_at_fixed_temperature() || !p.value_si().is_finite() {
        return Err(invalid("string memory requires dimensioned constant coefficients at the original material state"));
    }
    if let Some((lo,hi)) = p.angular_frequency_band_rad_s().map_err(owner)? {
        band.0=band.0.max(lo); band.1=band.1.min(hi);
    }
    Ok(p.value_si())
}

impl MaterialStringBankSource {
    /// Source sections and spectra survive multibow, forcing and PCM ownership.
    #[must_use]
    pub const fn string_bending_relaxation(&self) -> Option<&StringBankBendingRelaxation> {
        self.string_bending_relaxation.as_ref()
    }
}
impl ReciprocalStringBank {
    /// Add source-selected bending memory to homogeneous or bonded strings.
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
    ///
    /// Frozen uniform-span small-curvature bending only. No axial/torsional
    /// memory, finite-strain viscoelasticity, winding slip, heat evolution or
    /// measured instrument-fidelity claim is supplied by this selection.
    pub fn with_material_string_bending_relaxation(
        mut self, selections: Vec<StringBendingRelaxationSelection>,
        initial_memory: Option<Vec<f64>>, max_memory_states: usize, gate: &CancelGate,
    ) -> Result<Self, ReciprocalStringError> {
        poll(gate)?;
        let body_count = self.source.material_chart_relaxation.as_ref()
            .map_or(0, |image| image.coupling_rows_per_s().len());
        if self.system.samples_rendered()!=0 || self.source.string_bending_relaxation.is_some()
            || self.system.material_relaxation_branches().len()!=body_count {
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
