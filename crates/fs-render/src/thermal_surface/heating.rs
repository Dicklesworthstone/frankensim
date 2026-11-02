//! Resolved surface absorption on the original native thermal timeline.
//!
//! Each source supplies power already intercepted by its physical surface patch,
//! including any area, incidence cosine and quadrature weight. The retained
//! response partitions those watts once into reflection, transmission and true
//! absorption. Only absorption heats the explicitly assigned thermal cell.
//! Outgoing light remains the responsibility of its original transport owner.
//!
//! Responses are fixed or supplied by explicit native-temperature/enthalpy
//! substrate, film or diffuse laws. Coated sources assign the original coherent film and
//! substrate absorption to separate native cells. No phase mixing, within-film
//! deposition depth, emission or cooling law is inferred. Native caloric storage,
//! conduction and independently composed boundary owners remain authoritative.

mod response;
pub use response::SurfacePowerResponse;
/// Coherent coating absorption assigned to separate native thermal receivers.
pub mod coating;
pub use coating::{
    CoatedSurfaceIlluminationSource, CoatedSurfacePowerResponse, SurfaceHeatingReceiver,
    TransmittingSurfaceIlluminationSource, TransmittingSurfacePowerResponse,
};
/// Finite-film optical constants on the original temperature or enthalpy state.
pub mod coating_feedback;
pub use coating_feedback::SurfaceCoatingLaw;
/// Independently evolving film/substrate optics on the native thermal timeline.
pub mod coating_state;
use coating_state::CoatingFeedbackState;
pub use coating_state::{
    CoatingOpticsModel, StatefulCoatedSurfaceIlluminationSource, StatefulCoatingLayer,
    SubstrateOpticsModel,
};

/// Explicit conductor optical constants versus native temperature or enthalpy.
pub mod feedback;
mod state;
use state::SurfaceFeedbackState;
pub use state::{StatefulSurfaceIlluminationSource, SurfaceIncidence};

pub mod diffuse_feedback;
mod diffuse_state;
pub use diffuse_state::StatefulDiffuseSurfaceIlluminationSource;
use diffuse_state::SurfaceDiffuseFeedbackState;

use super::ThermalSurfaceError;
use crate::conductor::ConductorError;
use crate::dielectric::DielectricError;
use crate::volumes::heating::feedback::OpticalFeedbackError;
use crate::volumes::heating::illumination::schedule::SpectralPowerSchedule;
use crate::volumes::heating::{CheckedSum, OpticalHeatingError};
use fs_blake3::{ContentHash, DomainHasher};
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionObserver, DepositionSimulation, Event, MaterialId,
};
use std::collections::BTreeMap;
use std::sync::Arc;

/// One physical patch/source assignment. Several sources may heat the same cell.
/// All sources in one observer must use the same ordered wavelength packet.
#[derive(Clone, Debug)]
pub struct SurfaceIlluminationSource<const N: usize> {
    pub source_id: u64,
    pub thermal_cell: CellId,
    pub material_id: MaterialId,
    pub temperature_bounds_k: (f64, f64),
    pub response: SurfacePowerResponse<N>,
    pub power: SpectralPowerSchedule<N>,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct SurfaceHeatingBudget {
    pub max_sources: usize,
    pub max_total_power_segments: usize,
    /// Whole native advance, including dark intervals.
    pub max_source_evaluations: usize,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfaceHeatingTolerance {
    /// Absolute allowance for lane balance and final receiver reduction, W.
    /// This is arithmetic admission, not optical/material uncertainty.
    pub power_roundoff_w: f64,
}

/// Accepted spectral energy. Native external heat already includes thermal input.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SurfaceHeatingEnergy<const N: usize> {
    pub duration_s: f64,
    pub incident_j: [f64; N],
    pub reflected_j: [f64; N],
    pub transmitted_j: [f64; N],
    pub absorbed_j: [f64; N],
    pub residual_j: [f64; N],
    pub thermal_input_j: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum SurfaceHeatingError {
    Invalid(&'static str),
    Budget,
    Optics(ConductorError),
    Dielectric(DielectricError),
    Surface(ThermalSurfaceError),
    Profile(OpticalFeedbackError),
    Heating(OpticalHeatingError),
    Thermal(DepositionError),
    TemperatureDomain {
        thermal_cell: CellId,
        temperature_k: f64,
        minimum_k: f64,
        maximum_k: f64,
    },
}
impl From<ConductorError> for SurfaceHeatingError {
    fn from(value: ConductorError) -> Self {
        Self::Optics(value)
    }
}
impl From<DielectricError> for SurfaceHeatingError {
    fn from(value: DielectricError) -> Self {
        Self::Dielectric(value)
    }
}
impl From<ThermalSurfaceError> for SurfaceHeatingError {
    fn from(value: ThermalSurfaceError) -> Self {
        Self::Surface(value)
    }
}
impl From<OpticalFeedbackError> for SurfaceHeatingError {
    fn from(value: OpticalFeedbackError) -> Self {
        Self::Profile(value)
    }
}
impl From<OpticalHeatingError> for SurfaceHeatingError {
    fn from(value: OpticalHeatingError) -> Self {
        Self::Heating(value)
    }
}
impl From<DepositionError> for SurfaceHeatingError {
    fn from(value: DepositionError) -> Self {
        Self::Thermal(value)
    }
}
impl core::fmt::Display for SurfaceHeatingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "surface heating refused: {self:?}")
    }
}
impl std::error::Error for SurfaceHeatingError {}

#[derive(Clone, Copy, Debug)]
struct Receiver {
    material: MaterialId,
    minimum_k: f64,
    maximum_k: f64,
}

#[derive(Clone, Debug)]
struct PreparedSurfaceSource<const N: usize> {
    source_id: u64,
    /// Only a single-receiver source has an implicit destination. Coherent
    /// film partitions supply all actual receivers separately, including none
    /// for a bare transmitting interface.
    receiver: Option<SurfaceHeatingReceiver>,
    response: SurfacePowerResponse<N>,
    power: SpectralPowerSchedule<N>,
}

#[derive(Debug)]
struct Model<const N: usize> {
    identity: ContentHash,
    sources: Vec<PreparedSurfaceSource<N>>,
    receivers: BTreeMap<CellId, Receiver>,
    transitions_s: Vec<f64>,
    budget: SurfaceHeatingBudget,
    tolerance: SurfaceHeatingTolerance,
    /// Present only for sources with an explicit film/substrate partition.
    coating_receivers: BTreeMap<u64, Box<[coating::AbsorptionReceiver<N>]>>,
}

#[derive(Clone, Debug)]
struct Sums<const N: usize> {
    duration: CheckedSum,
    // Incident, reflected, transmitted, absorbed, residual; canonical lane order.
    channels: [[CheckedSum; N]; 5],
    thermal: CheckedSum,
}
impl<const N: usize> Default for Sums<N> {
    fn default() -> Self {
        Self {
            duration: CheckedSum::default(),
            channels: [[CheckedSum::default(); N]; 5],
            thermal: CheckedSum::default(),
        }
    }
}
impl<const N: usize> Sums<N> {
    fn value(&self) -> Result<SurfaceHeatingEnergy<N>, SurfaceHeatingError> {
        let mut values = [[0.0; N]; 5];
        for (sums, values) in self.channels.iter().zip(&mut values) {
            for (sum, value) in sums.iter().zip(values) {
                *value = sum.value()?;
            }
        }
        let [
            incident_j,
            reflected_j,
            transmitted_j,
            absorbed_j,
            residual_j,
        ] = values;
        Ok(SurfaceHeatingEnergy {
            duration_s: self.duration.value()?,
            incident_j,
            reflected_j,
            transmitted_j,
            absorbed_j,
            residual_j,
            thermal_input_j: self.thermal.value()?,
        })
    }
}

#[derive(Clone, Debug)]
struct Pending<const N: usize> {
    start_s: f64,
    duration_s: f64,
    end_s: f64,
    sums: Sums<N>,
    sources_w: Arc<BTreeMap<CellId, f64>>,
}

/// Surface absorption with exact source-profile intervals and native rollback.
/// Immutable source data and accepted maps are shared across candidate clones;
/// spectral integration and receiver reduction use the original compensated sum.
#[derive(Clone, Debug)]
pub struct SurfaceOpticalHeating<const N: usize> {
    model: Arc<Model<N>>,
    feedback: Option<SurfaceFeedbackState<N>>,
    diffuse_feedback: Option<SurfaceDiffuseFeedbackState<N>>,
    coating_feedback: Option<CoatingFeedbackState<N>>,
    time_s: f64,
    evaluations: usize,
    response_evaluations: usize,
    sums: Sums<N>,
    energy: SurfaceHeatingEnergy<N>,
    sources_w: Arc<BTreeMap<CellId, f64>>,
    pending: Option<Pending<N>>,
}

impl<const N: usize> SurfaceOpticalHeating<N> {
    /// Admit a nonempty bounded source set against an activated native state.
    /// Source IDs canonicalize reductions. Receiver material and temperature
    /// support are checked even for a currently dark or perfectly reflecting patch.
    pub fn try_new(
        sources: &[SurfaceIlluminationSource<N>],
        thermal: &DepositionSimulation<'_>,
        budget: SurfaceHeatingBudget,
        tolerance: SurfaceHeatingTolerance,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, SurfaceHeatingError> {
        checkpoint()?;
        if sources.len() > budget.max_sources {
            return Err(SurfaceHeatingError::Budget);
        }
        let mut prepared = Vec::new();
        prepared
            .try_reserve_exact(sources.len())
            .map_err(|_| SurfaceHeatingError::Budget)?;
        for source in sources {
            checkpoint()?;
            prepared.push(PreparedSurfaceSource {
                source_id: source.source_id,
                receiver: Some(SurfaceHeatingReceiver {
                    thermal_cell: source.thermal_cell,
                    material_id: source.material_id,
                    temperature_bounds_k: source.temperature_bounds_k,
                }),
                response: source.response,
                power: source.power.clone(),
            });
        }
        Self::try_new_prepared(prepared, thermal, budget, tolerance, checkpoint)
    }

    fn try_new_prepared(
        sources: Vec<PreparedSurfaceSource<N>>,
        thermal: &DepositionSimulation<'_>,
        budget: SurfaceHeatingBudget,
        tolerance: SurfaceHeatingTolerance,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, SurfaceHeatingError> {
        checkpoint()?;
        if N == 0
            || sources.is_empty()
            || !tolerance.power_roundoff_w.is_finite()
            || tolerance.power_roundoff_w < 0.0
        {
            return Err(invalid(
                "nonempty surface sources and finite nonnegative roundoff allowance required",
            ));
        }
        if sources.len() > budget.max_sources || budget.max_source_evaluations == 0 {
            return Err(SurfaceHeatingError::Budget);
        }
        let wavelengths = sources[0].response.wavelengths_nm().map(f64::to_bits);
        let mut segments = 0_usize;
        let mut receivers = BTreeMap::<CellId, Receiver>::new();
        for source in &sources {
            checkpoint()?;
            if source.response.wavelengths_nm().map(f64::to_bits) != wavelengths {
                return Err(invalid("one shared surface wavelength packet required"));
            }
            if let Some(binding) = source.receiver {
                let (lo, hi) = binding.temperature_bounds_k;
                if !lo.is_finite() || !hi.is_finite() || lo <= 0.0 || hi < lo {
                    return Err(invalid(
                        "finite positive surface temperature support required",
                    ));
                }
                let cell = thermal
                    .cells()
                    .get(&binding.thermal_cell)
                    .ok_or_else(|| invalid("surface receiver must be an active native cell"))?;
                if cell.material_id != binding.material_id {
                    return Err(invalid(
                        "surface receiver must retain its declared native material",
                    ));
                }
                let receiver = receivers.entry(binding.thermal_cell).or_insert(Receiver {
                    material: binding.material_id,
                    minimum_k: lo,
                    maximum_k: hi,
                });
                receiver.minimum_k = receiver.minimum_k.max(lo);
                receiver.maximum_k = receiver.maximum_k.min(hi);
            }
            segments = segments
                .checked_add(source.power.segment_count())
                .ok_or(SurfaceHeatingError::Budget)?;
            if segments > budget.max_total_power_segments {
                return Err(SurfaceHeatingError::Budget);
            }
        }
        let mut retained = sources;
        retained.sort_unstable_by_key(|source| source.source_id);
        if retained
            .windows(2)
            .any(|pair| pair[0].source_id == pair[1].source_id)
        {
            return Err(invalid("surface illumination source IDs must be unique"));
        }
        let mut transitions_s = Vec::new();
        transitions_s
            .try_reserve_exact(segments.checked_mul(2).ok_or(SurfaceHeatingError::Budget)?)
            .map_err(|_| SurfaceHeatingError::Budget)?;
        let explicit_receivers = retained.iter().any(|source| source.receiver.is_none());
        let mut hash = DomainHasher::new(if explicit_receivers {
            "org.frankensim.render.partitioned-surface-optical-heating.v1"
        } else {
            "org.frankensim.render.surface-optical-heating.v1"
        });
        hash.update(&(N as u64).to_le_bytes());
        hash.update(&(retained.len() as u64).to_le_bytes());
        for source in &retained {
            checkpoint()?;
            hash.update(&source.source_id.to_le_bytes());
            if explicit_receivers {
                hash.update(&[u8::from(source.receiver.is_some())]);
            }
            if let Some(receiver) = source.receiver {
                for id in [receiver.thermal_cell.0, receiver.material_id.0] {
                    hash.update(&id.to_le_bytes());
                }
                for value in [
                    receiver.temperature_bounds_k.0,
                    receiver.temperature_bounds_k.1,
                ] {
                    hash.update(&value.to_bits().to_le_bytes());
                }
            }
            hash.update(source.response.identity().as_bytes());
            hash.update(source.power.identity().as_bytes());
            for segment in source.power.segments() {
                checkpoint()?;
                transitions_s.extend([segment.start_s, segment.end_s]);
            }
        }
        transitions_s.sort_unstable_by(f64::total_cmp);
        transitions_s.dedup();
        let sums = Sums::default();
        let result = Self {
            model: Arc::new(Model {
                identity: hash.finalize(),
                sources: retained,
                receivers,
                transitions_s,
                budget,
                tolerance,
                coating_receivers: BTreeMap::new(),
            }),
            feedback: None,
            diffuse_feedback: None,
            coating_feedback: None,
            time_s: thermal.time_s(),
            evaluations: 0,
            response_evaluations: 0,
            energy: sums.value()?,
            sums,
            sources_w: Arc::new(BTreeMap::new()),
            pending: None,
        };
        result.validate_receivers(thermal, &mut checkpoint)?;
        checkpoint()?;
        Ok(result)
    }

    pub fn energy(&self) -> &SurfaceHeatingEnergy<N> {
        &self.energy
    }
    pub fn time_s(&self) -> f64 {
        self.time_s
    }
    pub fn wavelengths_nm(&self) -> &[f64; N] {
        self.model.sources[0].response.wavelengths_nm()
    }
    pub fn model_identity(&self) -> ContentHash {
        self.model.identity
    }
    pub fn source_evaluations(&self) -> usize {
        self.evaluations
    }
    /// Actual Fresnel packet rebuilds in this advance. Admission builds are
    /// excluded; unchanged constitutive values and dark packets reuse storage.
    pub fn response_evaluations(&self) -> usize {
        self.response_evaluations
    }
    pub fn source_count(&self) -> usize {
        self.model.sources.len()
    }
    /// Last accepted interval's combined native heat map; empty before any step.
    pub fn sources_w(&self) -> &BTreeMap<CellId, f64> {
        &self.sources_w
    }

    fn check_clock(&self, thermal: &DepositionSimulation<'_>) -> Result<(), SurfaceHeatingError> {
        if thermal.time_s() != self.time_s {
            return Err(invalid("surface and native thermal clocks differ"));
        }
        Ok(())
    }

    fn validate_receivers(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), SurfaceHeatingError> {
        for (&id, receiver) in &self.model.receivers {
            checkpoint()?;
            let cell = thermal
                .cells()
                .get(&id)
                .ok_or_else(|| invalid("surface receiver is no longer active"))?;
            if cell.material_id != receiver.material {
                return Err(invalid("surface receiver material changed"));
            }
            let t = thermal.temperature_k(id)?;
            if !t.is_finite() || !(receiver.minimum_k..=receiver.maximum_k).contains(&t) {
                return Err(SurfaceHeatingError::TemperatureDomain {
                    thermal_cell: id,
                    temperature_k: t,
                    minimum_k: receiver.minimum_k,
                    maximum_k: receiver.maximum_k,
                });
            }
        }
        if let Some(feedback) = &self.feedback {
            feedback.validate(thermal, checkpoint)?;
        }
        if let Some(feedback) = &self.diffuse_feedback {
            feedback.validate(thermal, checkpoint)?;
        }
        if let Some(feedback) = &self.coating_feedback {
            feedback.validate(thermal, checkpoint)?;
        }
        Ok(())
    }
}

impl<const N: usize> DepositionObserver for SurfaceOpticalHeating<N> {
    type Error = SurfaceHeatingError;

    fn begin(&mut self) -> Result<(), Self::Error> {
        self.evaluations = 0;
        self.response_evaluations = 0;
        self.pending = None;
        Ok(())
    }
    fn requires_step_interval(&self) -> bool {
        true
    }

    fn next_transition_time_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<f64>, Self::Error> {
        self.check_clock(thermal)?;
        checkpoint()?;
        let index = self
            .model
            .transitions_s
            .partition_point(|time| *time <= self.time_s);
        Ok(self.model.transitions_s.get(index).copied())
    }

    fn sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.check_clock(thermal)?;
        checkpoint()?;
        self.pending = None;
        Ok(BTreeMap::new())
    }

    fn step_sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        start_s: f64,
        duration_s: f64,
        end_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.check_clock(thermal)?;
        self.validate_receivers(thermal, checkpoint)?;
        if start_s != self.time_s {
            return Err(invalid(
                "surface source interval starts at the wrong native time",
            ));
        }
        if self.evaluations >= self.model.budget.max_source_evaluations {
            return Err(SurfaceHeatingError::Budget);
        }
        self.pending = None;
        if let Some(feedback) = &mut self.feedback {
            feedback.refresh(thermal, checkpoint)?;
        }
        if let Some(feedback) = &mut self.diffuse_feedback {
            feedback.refresh(thermal, checkpoint)?;
        }
        if let Some(feedback) = &mut self.coating_feedback {
            feedback.refresh(thermal, checkpoint)?;
        }
        let mut sums = self.sums.clone();
        let mut receiving = BTreeMap::<CellId, CheckedSum>::new();
        let mut absorbed = CheckedSum::default();
        let mut residuals = [CheckedSum::default(); N];
        for (index, source) in self.model.sources.iter().enumerate() {
            let mean = source
                .power
                .mean_power_w(start_s, duration_s, end_s, &mut *checkpoint)?;
            let illuminated = mean.iter().any(|watts| *watts > 0.0);
            let coating_packet = match &mut self.coating_feedback {
                Some(feedback) => feedback.response(
                    index,
                    illuminated,
                    &mut self.response_evaluations,
                    checkpoint,
                )?,
                None => None,
            };
            let response = if let Some(packet) = &coating_packet {
                packet.total
            } else {
                match &mut self.diffuse_feedback {
                    Some(feedback) => match feedback.response(
                        index,
                        illuminated,
                        &mut self.response_evaluations,
                        checkpoint,
                    )? {
                        Some(response) => response,
                        None => match &mut self.feedback {
                            Some(feedback) => feedback
                                .response(
                                    index,
                                    illuminated,
                                    &mut self.response_evaluations,
                                    checkpoint,
                                )?
                                .unwrap_or(source.response),
                            None => source.response,
                        },
                    },
                    None => match &mut self.feedback {
                        Some(feedback) => feedback
                            .response(
                                index,
                                illuminated,
                                &mut self.response_evaluations,
                                checkpoint,
                            )?
                            .unwrap_or(source.response),
                        None => source.response,
                    },
                }
            };
            let coating_receivers = coating_packet
                .as_ref()
                .map(|packet| packet.receivers.as_ref())
                .or_else(|| {
                    self.model
                        .coating_receivers
                        .get(&source.source_id)
                        .map(|receivers| receivers.as_ref())
                });
            for lane in 0..N {
                let incident = mean[lane];
                let reflected = physical_product(incident, response.reflected_fraction()[lane])?;
                let transmitted =
                    physical_product(incident, response.transmitted_fraction()[lane])?;
                let heat = physical_product(incident, response.absorbed_fraction()[lane])?;
                let mut residual = CheckedSum::default();
                for value in [incident, -reflected, -transmitted, -heat] {
                    residual.add(value)?;
                }
                let residual = residual.value()?;
                for (index, (channel, power)) in sums
                    .channels
                    .iter_mut()
                    .zip([incident, reflected, transmitted, heat, residual])
                    .enumerate()
                {
                    // The final channel is a signed arithmetic residual. It
                    // does not describe a physical positive power transfer.
                    let joules = if index < 4 {
                        physical_product(power, duration_s)?
                    } else {
                        power * duration_s
                    };
                    channel[lane].add(joules)?;
                }
                residuals[lane].add(residual)?;
                absorbed.add(heat)?;
                if let Some(receivers) = coating_receivers {
                    for receiver in receivers {
                        checkpoint()?;
                        let local_heat = physical_product(incident, receiver.absorbed[lane])?;
                        if local_heat > 0.0 {
                            receiving
                                .entry(receiver.cell)
                                .or_default()
                                .add(local_heat)?;
                        }
                    }
                } else if heat > 0.0 {
                    let receiver = source.receiver.ok_or_else(|| {
                        invalid("absorbing surface requires an explicit thermal receiver")
                    })?;
                    receiving
                        .entry(receiver.thermal_cell)
                        .or_default()
                        .add(heat)?;
                }
            }
        }
        let tolerance = self.model.tolerance.power_roundoff_w;
        for residual in residuals {
            if residual.value()?.abs() > tolerance {
                return Err(invalid(
                    "surface incident/reflected/transmitted/absorbed power is unbalanced",
                ));
            }
        }
        let mut sources_w = BTreeMap::new();
        let mut thermal_w = CheckedSum::default();
        for (id, sum) in receiving {
            checkpoint()?;
            let watts = sum.value()?;
            thermal_w.add(watts)?;
            sources_w.insert(id, watts);
        }
        let thermal_w = thermal_w.value()?;
        if (thermal_w - absorbed.value()?).abs() > tolerance {
            return Err(invalid(
                "surface absorption and native receiver power disagree",
            ));
        }
        sums.duration.add(duration_s)?;
        sums.thermal.add(physical_product(thermal_w, duration_s)?)?;
        sums.value()?;
        checkpoint()?;
        self.pending = Some(Pending {
            start_s,
            duration_s,
            end_s,
            sums,
            sources_w: Arc::new(sources_w.clone()),
        });
        self.evaluations += 1;
        Ok(sources_w)
    }

    fn event(&mut self, _: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        self.check_clock(thermal)?;
        self.validate_receivers(thermal, &mut || Ok(()))
    }

    fn step(
        &mut self,
        start_s: f64,
        duration_s: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        let pending = self
            .pending
            .as_ref()
            .ok_or_else(|| invalid("surface acceptance needs prepared interval power"))?;
        if start_s != self.time_s
            || pending.start_s != start_s
            || pending.duration_s != duration_s
            || pending.end_s != thermal.time_s()
        {
            return Err(invalid(
                "accepted native interval differs from surface heating proposal",
            ));
        }
        self.validate_receivers(thermal, checkpoint)?;
        let energy = pending.sums.value()?;
        checkpoint()?;
        let pending = self.pending.take().expect("validated surface interval");
        self.sums = pending.sums;
        self.energy = energy;
        self.sources_w = pending.sources_w;
        self.time_s = thermal.time_s();
        Ok(())
    }
}

fn invalid(message: &'static str) -> SurfaceHeatingError {
    SurfaceHeatingError::Invalid(message)
}

/// A nonzero physical channel cannot silently become zero before native heat
/// sees it. Existing checked reductions still own nonfinite/overflow refusal.
fn physical_product(a: f64, b: f64) -> Result<f64, SurfaceHeatingError> {
    let value = a * b;
    if a > 0.0 && b > 0.0 && value == 0.0 {
        Err(invalid(
            "positive surface power or energy is below the representable range",
        ))
    } else {
        Ok(value)
    }
}
