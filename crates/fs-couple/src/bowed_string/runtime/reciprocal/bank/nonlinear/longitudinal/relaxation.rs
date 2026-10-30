//! Linear axial hereditary stress on the original retained longitudinal modes.
//!
//! A bonded layer contributes `Delta E * A * k_n^2 / mu_solid` to the
//! mass-normalized axial stiffness. Its relaxation time and independent history
//! survive; bending's area moment and transverse fluid-added mass never enter
//! this projection. The existing section owner admits the original uniaxial
//! Young-modulus spectra; the existing shared Maxwell owner advances them.
use super::*;
use crate::bowed_string::runtime::reciprocal::bank::{
    ReciprocalStringError, invalid as bank_invalid,
};
use crate::render::schedule::force::coupled::contact::friction::FrictionModalSystem;
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::section::relaxation::{ConcentricPronySection, LayerRelaxationProperties};
use fs_phs::RelaxationBranch;

/// One spectrum per original physical layer, for one retained axial member.
#[derive(Clone, Debug, PartialEq)]
pub struct LongitudinalRelaxationSelection {
    pub string: usize,
    /// Core first, followed by bonded sleeves. An empty branch list explicitly
    /// selects an elastic layer, whose equilibrium claim still needs a band.
    pub layers: Vec<LayerRelaxationProperties>,
}

/// One independent layer/arm/axial-mode history in the original basis.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LongitudinalMemoryChannel {
    pub string: usize,
    pub layer: usize,
    pub arm: usize,
    pub axial_component: usize,
    pub mode: usize,
    /// One-based physical sine index, not its retained-vector position.
    pub sine_index: usize,
    /// `sqrt(Delta E * A / mu_solid) * k_n`, in inverse seconds.
    pub coupling_per_s: f64,
    pub relaxation_time_s: f64,
}

/// Immutable source and projection. Actual evolving histories remain in the
/// same coupled mechanical state as strings, contacts, body and cavity.
#[derive(Clone, Debug)]
pub struct StringBankLongitudinalRelaxation {
    selections: Vec<LongitudinalRelaxationSelection>,
    sections: Vec<ConcentricPronySection>,
    channels: Vec<LongitudinalMemoryChannel>,
    basis: Vec<(usize, Vec<ModalAcousticMode>)>,
    max_memory_states: usize,
}
impl StringBankLongitudinalRelaxation {
    pub fn selections(&self) -> &[LongitudinalRelaxationSelection] { &self.selections }
    /// Original uniaxial spectra, equilibrium sections and complete receipts.
    /// This projection uses each layer's AREA, not the bending moment carried
    /// by the common section-spectrum representation.
    pub fn sections(&self) -> &[ConcentricPronySection] { &self.sections }
    pub fn channels(&self) -> &[LongitudinalMemoryChannel] { &self.channels }
    pub const fn max_memory_states(&self) -> usize { self.max_memory_states }

    fn bind(bank: &ReciprocalStringBank, selections: Vec<LongitudinalRelaxationSelection>,
        max_memory_states: usize, gate: &CancelGate) -> Result<Self, ImpactError>
    {
        let axial = bank.source.longitudinal_waves.as_ref()
            .ok_or_else(|| invalid("axial memory requires retained original longitudinal waves"))?;
        if selections.is_empty() || selections.len() > axial.strings.len()
            || !(1..=256).contains(&max_memory_states) {
            return Err(invalid("axial memory needs distinct source selections and 1..=256 history states"));
        }
        let mut sections = Vec::with_capacity(selections.len());
        let mut channels = Vec::new();
        let mut basis = Vec::with_capacity(selections.len());
        for (index, selection) in selections.iter().enumerate() {
            if gate.is_requested() { return Err(ImpactError::Cancelled); }
            let member = bank.source.strings.get(selection.string)
                .ok_or_else(|| invalid("axial memory names an absent original string"))?;
            let law = &axial.strings[selection.string];
            if selections[..index].iter().any(|s| s.string == selection.string) {
                return Err(invalid("one axial member cannot select two material-memory owners"));
            }
            // Density, equilibrium E, every Delta E/tau and original layer
            // ordering are validated once by the shared material owner.
            let section = elastic_section(&member.source)?.with_prony_bending(&selection.layers)
                .map_err(|e| ImpactError::Owner(e.to_string()))?;
            // Homogeneous specimens compute pi*(r*r), while the annular owner
            // factors pi*(r-0)*(r+0). Preserve the original mass basis and admit
            // only its floating-point rounding, not a new equivalent section.
            let same = |a: f64, b: f64| (a - b).abs() <= 16.0 * f64::EPSILON * a.abs().max(b.abs());
            if !same(section.elastic().line_density_kg_m(), law.solid_linear_density_kg_m)
                || !same(section.elastic().axial_stiffness_n(), law.axial_rigidity_n) {
                return Err(invalid("axial memory must preserve its original solid inertia and equilibrium rigidity"));
            }
            let component = &bank.system.components()[law.axial_component];
            if component.modes().len() != law.model.axial_modes().len() {
                return Err(invalid("axial memory lost its complete original sine basis"));
            }
            basis.push((law.axial_component, component.modes().to_vec()));
            for (layer, (spectrum, geometry)) in section.spectra().iter()
                .zip(section.elastic().layer_integrals()).enumerate() {
                for (arm, &(delta, tau)) in spectrum.terms.iter().enumerate() {
                    if gate.is_requested() { return Err(ImpactError::Cancelled); }
                    if delta == 0.0 { continue; }
                    if channels.len().checked_add(law.model.axial_modes().len())
                        .is_none_or(|n| n > max_memory_states) {
                        return Err(invalid("complete per-layer axial histories exceed their declared budget"));
                    }
                    let delta_ea = delta * geometry.area_m2();
                    let scale = det::sqrt(delta_ea / law.solid_linear_density_kg_m);
                    if !delta_ea.is_finite() || delta_ea <= 0.0 || !scale.is_finite() || scale <= 0.0 {
                        return Err(invalid("axial material-memory rigidity is unrepresentable"));
                    }
                    for (mode, &sine_index) in law.model.axial_modes().iter().enumerate() {
                        let k = sine_index as f64 * core::f64::consts::PI / bank.cards[selection.string].length_m;
                        let coupling = scale * k;
                        if !coupling.is_finite() || coupling <= 0.0 {
                            return Err(invalid("axial strain projection lost a positive source channel"));
                        }
                        channels.push(LongitudinalMemoryChannel { string: selection.string,
                            layer, arm, axial_component: law.axial_component, mode, sine_index,
                            coupling_per_s: coupling, relaxation_time_s: tau });
                    }
                }
            }
            sections.push(section);
        }
        Ok(Self { selections, sections, channels, basis, max_memory_states })
    }

    pub(in crate::bowed_string::runtime::reciprocal::bank) fn validate_band(&self,
        high: f64, dt: f64) -> Result<(), ReciprocalStringError>
    {
        if !high.is_finite() || high < 0.0 || !dt.is_finite() || dt <= 0.0
            || self.sections.iter().any(|section| {
                let (lo, hi) = section.omega_band_rad_s(); lo > 0.0 || high > hi
            }) || self.sections.iter().flat_map(|section| section.spectra())
                .flat_map(|spectrum| &spectrum.terms).any(|&(_, tau)| dt / tau > 2.0) {
            return Err(bank_invalid("complete coupled band or clock exceeds original axial-memory support"));
        }
        Ok(())
    }

    pub(in crate::bowed_string::runtime::reciprocal::bank) fn branches(&self,
        system: &FrictionModalSystem, gate: &CancelGate) -> Result<Vec<RelaxationBranch>, ReciprocalStringError>
    {
        for (component, modes) in &self.basis {
            if system.components().get(*component).is_none_or(|c| c.modes() != modes.as_slice()) {
                return Err(bank_invalid("axial material history requires its unchanged original solid basis"));
            }
        }
        let dimension = system.mode_count().checked_mul(2)
            .ok_or_else(|| bank_invalid("axial history embedding overflow"))?;
        self.channels.iter().map(|channel| {
            if gate.is_requested() { return Err(ReciprocalStringError::Cancelled); }
            let start: usize = system.components()[..channel.axial_component].iter().map(|c| c.modes().len()).sum();
            let mut projection = vec![0.0; dimension];
            projection[2 * (start + channel.mode)] = channel.coupling_per_s;
            Ok(RelaxationBranch { projection, stiffness: 1.0, relaxation_time_s: channel.relaxation_time_s })
        }).collect()
    }
}

fn elastic_section(source: &MaterialStringSource) -> Result<ConcentricElasticSection, ImpactError> {
    match source {
        MaterialStringSource::WithPolarizations(p) => elastic_section(p.source()),
        MaterialStringSource::WithTorsion(t) => elastic_section(t.bending()),
        MaterialStringSource::Homogeneous(s) => ConcentricElasticSection::try_new(0.0,
            vec![ConcentricMaterialLayer { outer_radius_m: s.radius_m(), material: s.material().clone() }])
            .map_err(|e| ImpactError::Owner(e.to_string())),
        MaterialStringSource::Concentric { section, .. } => Ok(section.elastic().clone()),
        MaterialStringSource::ConcentricProny { section, .. } => Ok(section.elastic().clone()),
        MaterialStringSource::Wound(_) => Err(invalid("quasistatic winding spring has no resolved axial material history")),
    }
}

fn bank_error(error: ReciprocalStringError) -> ImpactError {
    match error {
        ReciprocalStringError::Cancelled => ImpactError::Cancelled,
        error => ImpactError::Owner(error.to_string()),
    }
}

impl ReciprocalStringBank {
    /// Attach source-derived axial Maxwell arms to the SAME retained waves,
    /// nonlinear transverse strain, reciprocal bearings, material body and sound.
    /// Every original layer/arm/sine coordinate gets an independent history.
    /// Existing bending, torsional, body and gas histories stay unchanged.
    ///
    /// Bind after longitudinal waves and other source memory, before playback.
    /// `None` initializes the added arms relaxed; explicit memory follows
    /// channels() order. This is linear u_x hereditary stress about a maintained
    /// equilibrium preload. It does not relax the prescribed tension or supply
    /// finite-slope hereditary y_x*y_xt stress, creep, self-heating or wrap slip.
    /// A separately selected axial Kelvin--Voigt law is an additional parallel
    /// mechanism, not a conversion of these spectra to damping ratios.
    ///
    /// The complete instantaneous band, original clock, exact-state capacity
    /// and explicit memory allowance are checked. No waveform fitting, extra
    /// stepping, missing-data fallback or frequency-dependent output filter exists.
    pub fn with_material_longitudinal_relaxation(mut self,
        selections: Vec<LongitudinalRelaxationSelection>, initial_memory: Option<Vec<f64>>,
        max_memory_states: usize, gate: &CancelGate) -> Result<Self, ImpactError>
    {
        if gate.is_requested() { return Err(ImpactError::Cancelled); }
        if self.system.samples_rendered() != 0 || self.source.longitudinal_waves.as_ref()
            .is_none_or(|axial| axial.material_relaxation.is_some()) {
            return Err(invalid("axial material memory requires one cold binding to retained longitudinal waves"));
        }
        let image = StringBankLongitudinalRelaxation::bind(&self, selections, max_memory_states, gate)?;
        let branches = image.branches(&self.system, gate).map_err(bank_error)?;
        self.system = self.system.with_additional_material_relaxation(branches, initial_memory, gate)
            .map_err(|e| ImpactError::Owner(e.to_string()))?;
        image.validate_band(self.system.linear_frequency_upper_bound_rad_s(), self.system.sample_period_s())
            .map_err(bank_error)?;
        self.source.longitudinal_waves.as_mut()
            .ok_or_else(|| invalid("axial material memory lost its original wave owner"))?
            .material_relaxation = Some(image);
        self.source.validate_body_memory_frequencies(&self.system, gate)
            .map_err(bank_error)?;
        self.source.validate_material_bridge_frequencies(&self.system, gate)
            .map_err(bank_error)?;
        if gate.is_requested() { return Err(ImpactError::Cancelled); }
        Ok(self)
    }
}
