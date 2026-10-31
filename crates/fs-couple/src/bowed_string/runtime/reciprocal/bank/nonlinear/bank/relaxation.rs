//! Original material spectra on Kirchhoff--Carrier's shared mean strain.
//!
//! Core arms use Delta E*A; a close-coiled winding uses its existing spring
//! geometry times Delta G. Only hereditary stress/storage is added: equilibrium
//! EA, preload, mass, bending, gas, contacts and the time owner remain unchanged.
use super::*;
use crate::bowed_string::runtime::reciprocal::bank::{
    ReciprocalStringError, invalid as bank_invalid,
};
use crate::render::schedule::force::coupled::contact::friction::FrictionModalSystem;
use crate::render::schedule::force::coupled::contact::friction::nonlinear::{
    ModalCoordinate, history::MappedModalHistoryBinding,
};
use fs_material::section::relaxation::{ConcentricPronySection, LayerRelaxationProperties};
use fs_material::section::winding::{CloseWoundSection, relaxation::CloseWoundRelaxation};
use fs_nlmodal::longitudinal::relaxation::{StrainRelaxation, StrainRelaxationArm};

mod history;
use history::MeanStrainHistory;

/// Explicit constitutive roles for one original member. Wound sections require
/// an explicit winding selection, even when its branch list is empty.
#[derive(Clone, Debug, PartialEq)]
pub struct StretchingRelaxationSelection {
    pub string: usize,
    pub layers: Vec<LayerRelaxationProperties>,
    pub winding: Option<LayerRelaxationProperties>,
}

/// Source identity and actual section integrals survive the mechanical lowering.
#[derive(Clone, Debug)]
pub struct StretchingRelaxationSource {
    pub core: ConcentricPronySection,
    pub winding: Option<CloseWoundRelaxation>,
}

/// One original nonzero arm, in selection/core-layer/arm then winding-arm order.
/// Native heat-channel selectors can be obtained from the retained binding and
/// channel indices; losses are computed by each actual arm, not partitioned.
#[derive(Clone, Debug)]
pub struct StretchingMemoryChannel {
    pub string: usize,
    pub layer: usize,
    pub arm: usize,
    pub winding: bool,
    pub relaxing_rigidity_n: f64,
    pub relaxation_time_s: f64,
    pub history_binding: usize,
    pub heat_channel: usize,
}

/// Source/projection only. Every accepted relaxed strain stays in the original
/// nonlinear mechanical owner, together with body, contact and fluid histories.
#[derive(Clone, Debug)]
pub struct StringBankStretchingRelaxation {
    selections: Vec<StretchingRelaxationSelection>,
    sources: Vec<StretchingRelaxationSource>,
    channels: Vec<StretchingMemoryChannel>,
}
impl StringBankStretchingRelaxation {
    pub fn selections(&self) -> &[StretchingRelaxationSelection] {
        &self.selections
    }
    pub fn sources(&self) -> &[StretchingRelaxationSource] {
        &self.sources
    }
    pub fn channels(&self) -> &[StretchingMemoryChannel] {
        &self.channels
    }

    /// Read the very same accepted arm strains used by the mechanical owner.
    /// Geometric strain/slope are unchanged; only actual axial stress/storage
    /// are added, without staging or advancing any material state.
    pub(super) fn add_observation(
        &self,
        string: usize,
        length_m: f64,
        memory: &[&[f64]],
        observation: &mut StringObservation,
    ) -> Result<(), ImpactError> {
        for channel in self
            .channels
            .iter()
            .filter(|channel| channel.string == string)
        {
            let z = memory
                .get(channel.history_binding)
                .and_then(|history| history.get(channel.heat_channel))
                .ok_or_else(|| {
                    invalid("mean-strain observation lost its original accepted material arm")
                })?;
            let strain = observation.additional_strain - z;
            observation.tension_n += channel.relaxing_rigidity_n * strain;
            observation.stretching_energy_j +=
                0.5 * channel.relaxing_rigidity_n * length_m * strain * strain;
        }
        if !observation.tension_n.is_finite()
            || observation.tension_n <= 0.0
            || !observation.stretching_energy_j.is_finite()
        {
            return Err(invalid(
                "mean-strain observation requires finite material storage and positive total tension",
            ));
        }
        Ok(())
    }

    pub(in crate::bowed_string::runtime::reciprocal::bank) fn validate_band(
        &self,
        high: f64,
        dt: f64,
    ) -> Result<(), ReciprocalStringError> {
        if !high.is_finite()
            || high < 0.0
            || !dt.is_finite()
            || dt <= 0.0
            || high * dt >= core::f64::consts::PI
        {
            return Err(bank_invalid(
                "mean-strain material history requires the original resolved mechanical clock",
            ));
        }
        for source in &self.sources {
            for (lo, hi) in std::iter::once(source.core.omega_band_rad_s())
                .chain(source.winding.as_ref().map(|w| w.omega_band_rad_s()))
            {
                if lo > 0.0 || hi < high {
                    return Err(bank_invalid(
                        "complete coupled band exceeds original core or winding relaxation support",
                    ));
                }
            }
        }
        Ok(())
    }
}

fn winding(source: &MaterialStringSource) -> Option<&CloseWoundSection> {
    match source {
        MaterialStringSource::WithPolarizations(p) => winding(p.source()),
        MaterialStringSource::WithTorsion(t) => winding(t.bending()),
        MaterialStringSource::Wound(s) => Some(s.section()),
        _ => None,
    }
}

pub(in crate::bowed_string::runtime::reciprocal::bank) fn bindings(
    source: &MaterialStringBankSource,
    cards: &[crate::bowed_string::BowedStringCard],
    system: &FrictionModalSystem,
    selections: Vec<StretchingRelaxationSelection>,
    initial_memory: Option<Vec<f64>>,
    maximum_states: usize,
    gate: &CancelGate,
) -> Result<
    (
        StringBankStretchingRelaxation,
        Vec<MappedModalHistoryBinding>,
    ),
    ImpactError,
> {
    if gate.is_requested() {
        return Err(ImpactError::Cancelled);
    }
    let stretching = source.geometric_stretching.as_ref().ok_or_else(|| {
        invalid("mean-strain material memory requires original geometric stretching")
    })?;
    if system.samples_rendered() != 0
        || stretching.material_relaxation.is_some()
        || source.longitudinal_waves.is_some()
        || selections.is_empty()
        || selections.len() > source.strings.len()
        || !(1..=256).contains(&maximum_states)
    {
        return Err(invalid(
            "mean-strain memory needs one bounded cold source binding after static initialization",
        ));
    }
    let mut sources = Vec::with_capacity(selections.len());
    let mut channels = Vec::new();
    let mut laws = Vec::new();
    let first_binding = system.nonlinear_material_memory().len();
    for (selection_index, selection) in selections.iter().enumerate() {
        if gate.is_requested() {
            return Err(ImpactError::Cancelled);
        }
        let member = source
            .strings
            .get(selection.string)
            .ok_or_else(|| invalid("mean-strain memory names an absent original member"))?;
        if selections[..selection_index]
            .iter()
            .any(|s| s.string == selection.string)
        {
            return Err(invalid("mean-strain memory repeats a member"));
        }
        let wound = winding(&member.source);
        let elastic = match wound {
            Some(section) => section.core().elastic().clone(),
            None => super::super::longitudinal::relaxation::elastic_section(&member.source)?,
        };
        let core = elastic
            .with_prony_bending(&selection.layers)
            .map_err(|e| ImpactError::Owner(e.to_string()))?;
        let winding = match (wound, &selection.winding) {
            (Some(section), Some(properties)) => Some(
                section
                    .clone()
                    .with_winding_relaxation(properties.clone())
                    .map_err(|e| ImpactError::Owner(e.to_string()))?,
            ),
            (None, None) => None,
            _ => {
                return Err(invalid(
                    "winding memory must explicitly match the original helical geometry",
                ));
            }
        };
        let equilibrium = core.elastic().axial_stiffness_n()
            + winding
                .as_ref()
                .map_or(0.0, |w| w.section().winding_axial_stiffness_n());
        let original = member.source.axial_rigidity_n()?;
        if (equilibrium - original).abs()
            > 16.0 * f64::EPSILON * equilibrium.abs().max(original.abs())
        {
            return Err(invalid(
                "material memory cannot replace original equilibrium axial rigidity",
            ));
        }
        let length_m = cards[selection.string].length_m;
        let mut arms = Vec::new();
        let mut add = |layer, arm, winding, rigidity: f64, tau| -> Result<(), ImpactError> {
            if rigidity == 0.0 {
                return Ok(());
            }
            let stiffness_j = rigidity * length_m;
            if !stiffness_j.is_finite() || stiffness_j <= 0.0 || channels.len() >= maximum_states {
                return Err(invalid(
                    "complete source mean-strain arms exceed storage or memory limits",
                ));
            }
            channels.push(StretchingMemoryChannel {
                string: selection.string,
                layer,
                arm,
                winding,
                relaxing_rigidity_n: rigidity,
                relaxation_time_s: tau,
                history_binding: first_binding + laws.len(),
                heat_channel: arms.len(),
            });
            arms.push(StrainRelaxationArm {
                stiffness_j,
                relaxation_time_s: tau,
            });
            Ok(())
        };
        for (layer, (spectrum, geometry)) in core
            .spectra()
            .iter()
            .zip(core.elastic().layer_integrals())
            .enumerate()
        {
            for (arm, &(delta, tau)) in spectrum.terms.iter().enumerate() {
                let rigidity = delta * geometry.area_m2();
                if delta > 0.0 && rigidity == 0.0 {
                    return Err(invalid("positive axial relaxation rigidity underflow"));
                }
                add(layer, arm, false, rigidity, tau)?;
            }
        }
        if let Some(winding) = &winding {
            for (arm, source) in winding.arms().iter().enumerate() {
                add(
                    core.elastic().layers().len(),
                    arm,
                    true,
                    source.relaxing_rigidity_n,
                    source.relaxation_time_s,
                )?;
            }
        }
        if !arms.is_empty() {
            let material =
                StrainRelaxation::new(0, arms).map_err(|e| ImpactError::Owner(e.to_string()))?;
            laws.push((
                selection.string,
                MeanStrainHistory {
                    potential: stretching.potentials[selection.string].clone(),
                    length_m,
                    material,
                },
            ));
        }
        sources.push(StretchingRelaxationSource { core, winding });
    }
    if initial_memory
        .as_ref()
        .is_some_and(|m| m.len() != channels.len() || m.iter().any(|v| !v.is_finite()))
    {
        return Err(invalid(
            "complete finite initial relaxed mean strains required",
        ));
    }
    let mut start = 0;
    let mut bindings = Vec::with_capacity(laws.len());
    for (string, law) in laws {
        let end = start + law.material.memory_count();
        let coordinates = (0..law.potential.mode_count())
            .map(|mode| ModalCoordinate {
                component: string,
                mode,
            })
            .collect();
        bindings.push(MappedModalHistoryBinding {
            coordinates,
            initial_memory: initial_memory.as_ref().map(|m| m[start..end].to_vec()),
            law: Box::new(law),
        });
        start = end;
    }
    if gate.is_requested() {
        return Err(ImpactError::Cancelled);
    }
    Ok((
        StringBankStretchingRelaxation {
            selections,
            sources,
            channels,
        },
        bindings,
    ))
}

pub(in crate::bowed_string::runtime::reciprocal::bank) fn retain_source(
    source: &mut MaterialStringBankSource,
    image: StringBankStretchingRelaxation,
) {
    source
        .geometric_stretching
        .as_mut()
        .expect("admitted original geometric law")
        .material_relaxation = Some(image);
}

impl ReciprocalStringBank {
    /// Bind source axial and optional winding-shear memory after geometric
    /// stretching and static initialization, before playback. Original source
    /// equilibrium EA stays in the existing potential exactly once. Histories
    /// act on BOTH transverse planes' common mean strain, using the original
    /// joint nonlinear/contact step and independently computed per-arm heat.
    ///
    /// This is the admitted quasistatic axial reduction. It adds no resolved
    /// winding/longitudinal inertia, interturn contact, friction or self-heating.
    pub fn with_material_stretching_relaxation(
        mut self,
        selections: Vec<StretchingRelaxationSelection>,
        initial_memory: Option<Vec<f64>>,
        maximum_states: usize,
        gate: &CancelGate,
    ) -> Result<Self, ImpactError> {
        let (image, bindings) = bindings(
            &self.source,
            &self.cards,
            &self.system,
            selections,
            initial_memory,
            maximum_states,
            gate,
        )?;
        if !bindings.is_empty() {
            self.system = self
                .system
                .with_mapped_material_histories(bindings, maximum_states, gate)
                .map_err(|e| ImpactError::Owner(e.to_string()))?;
        }
        retain_source(&mut self.source, image);
        self.source
            .validate_body_memory_frequencies(&self.system, gate)
            .map_err(|e| ImpactError::Owner(e.to_string()))?;
        self.source
            .validate_material_bridge_frequencies(&self.system, gate)
            .map_err(|e| ImpactError::Owner(e.to_string()))?;
        if gate.is_requested() {
            return Err(ImpactError::Cancelled);
        }
        Ok(self)
    }
}
