//! Prescribed photon illumination on the original conservative heat timeline.
//!
//! Independently named sources carry incident group powers, not invented nuclear
//! inventories. Each final thermal interval integrates its affine power history
//! analytically, applies the original immutable photon response, and credits
//! absorbed joules once. The existing thermal owner retains calorics, conduction,
//! boundaries and atomic rollback. Transport is steady; optional original
//! temperature tables refresh its material kernel at native step boundaries.
//! Photon flight time and evolving density/geometry are not inferred.

use super::photon::{photon_call, photon_material_input};
use super::radiation::bind_receivers;
use super::*;
use crate::endf::decay::spectra::DecayRadiation;
use crate::radiation::grouped::GroupedRadiationResponse;
use crate::radiation::transport::refinement::{
    PhotonSlabRefinement, RefinedPhotonSlabResponse, refine_slab_response_with_sources,
};
use crate::radiation::transport::{
    PhotonOrdinate, PhotonSlabLimits, PhotonSlabResponse, PhotonSlabSource,
    prepare_slab_response_with_sources,
};
use fs_math::dd::Dd;

mod affine;
use affine::{affine_energy, interpolate};

mod feedback;
pub use feedback::{PhotonFeedbackBudget, PhotonFeedbackWork, TemperaturePhotonSlabLayer};
use feedback::{PhotonFeedbackCandidate, PhotonTemperatureFeedback};

const MAX_SCHEDULE_VALUES: usize = 1_048_576;

/// One knot of total group powers already entering the declared source geometry.
#[derive(Clone, Debug, PartialEq)]
pub struct PhotonPowerKnot {
    /// Absolute time [s]. The first knot is zero; later knots strictly increase.
    pub time_s: f64,
    /// Nonnegative group-integrated watts in the original energy-grid order.
    pub group_power_w: Vec<f64>,
}

/// Immutable, continuous piecewise-affine photon power on one explicit grid.
/// The final knot is held afterward; a finite pulse therefore ends at zero.
/// This is already incident energy per second, with no extra area/cosine factor.
#[derive(Clone, Debug)]
pub struct PhotonPowerSchedule {
    data_id: Arc<str>,
    radiation: DecayRadiation,
    edges_j: Arc<[f64]>,
    knots: Arc<[PhotonPowerKnot]>,
    values: usize,
}

impl PhotonPowerSchedule {
    /// Admit the complete source spectrum, absolute-time knots and scalar budget.
    /// No group rebinning, isotope inference or spectral normalization is done.
    pub fn new(
        data_id: String,
        radiation: DecayRadiation,
        group_edges_j: Vec<f64>,
        knots: Vec<PhotonPowerKnot>,
        max_values: usize,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Self, FissionThermalError> {
        checkpoint()?;
        let groups = group_edges_j
            .len()
            .checked_sub(1)
            .filter(|&n| n > 0)
            .ok_or(FissionThermalError::Invalid(
                "photon power needs an energy grid",
            ))?;
        let values = groups
            .checked_mul(knots.len())
            .filter(|&n| n > 0 && n <= max_values && n <= MAX_SCHEDULE_VALUES)
            .ok_or(FissionThermalError::Invalid(
                "photon power schedule exceeds scalar allowance",
            ))?;
        if data_id.trim().is_empty()
            || !matches!(radiation, DecayRadiation::Gamma | DecayRadiation::XRay)
            || !group_edges_j.iter().copied().all(nonnegative)
            || !group_edges_j.windows(2).all(|w| w[0] < w[1])
            || knots[0].time_s != 0.0
        {
            return Err(FissionThermalError::Invalid(
                "invalid photon power source or grid",
            ));
        }
        let mut previous = None;
        for knot in &knots {
            checkpoint()?;
            if !nonnegative(knot.time_s)
                || previous.is_some_and(|time| knot.time_s <= time)
                || knot.group_power_w.len() != groups
            {
                return Err(FissionThermalError::Invalid("invalid photon power knot"));
            }
            for chunk in knot.group_power_w.chunks(256) {
                checkpoint()?;
                if !chunk.iter().copied().all(nonnegative) {
                    return Err(FissionThermalError::Invalid(
                        "photon group powers must be nonnegative",
                    ));
                }
            }
            previous = Some(knot.time_s);
        }
        Ok(Self {
            data_id: data_id.into(),
            radiation,
            edges_j: group_edges_j.into(),
            knots: knots.into(),
            values,
        })
    }

    /// Authored source identity, independent of the transport response identity.
    pub fn data_id(&self) -> &str {
        &self.data_id
    }
    /// Original immutable power knots; no resampling occurs during preparation.
    pub fn knots(&self) -> &[PhotonPowerKnot] {
        &self.knots
    }
    /// Original photon energy edges [J].
    pub fn group_edges_j(&self) -> &[f64] {
        &self.edges_j
    }
    /// Original gamma or X-ray kind.
    pub fn radiation(&self) -> DecayRadiation {
        self.radiation
    }

    fn next_after(&self, time_s: f64) -> Option<f64> {
        self.knots
            .get(self.knots.partition_point(|k| k.time_s <= time_s))
            .map(|k| k.time_s)
    }

    // The thermal owner's final duration is authoritative; exact end_s chooses
    // the power endpoint, including native event-time rounding.
    fn interval(
        &self,
        start_s: f64,
        duration_s: f64,
        end_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Vec<f64>, f64), FissionThermalError> {
        let index = self.knots.partition_point(|k| k.time_s <= start_s) - 1;
        let left = &self.knots[index];
        let right = self.knots.get(index + 1);
        if right.is_some_and(|k| k.time_s < end_s) {
            return Err(FissionThermalError::Invalid(
                "photon interval crosses a power knot",
            ));
        }
        let mut energy = Vec::with_capacity(left.group_power_w.len());
        let mut roundoff = Sum::default();
        for (group, &power) in left.group_power_w.iter().enumerate() {
            if group % 64 == 0 {
                checkpoint()?;
            }
            let (p0, p1) = if let Some(right) = right {
                let q = right.group_power_w[group];
                (
                    interpolate(power, q, left.time_s, right.time_s, start_s),
                    interpolate(power, q, left.time_s, right.time_s, end_s),
                )
            } else {
                (Dd::from_f64(power), Dd::from_f64(power))
            };
            let (value, error) =
                affine_energy(p0, p1, duration_s).ok_or(FissionThermalError::Arithmetic)?;
            energy.push(value);
            roundoff.add(error)?;
        }
        Ok((energy, roundoff.nonnegative_upper_bound()?))
    }
}

/// One independently identified source geometry and its prescribed spectrum.
#[derive(Clone, Debug)]
pub struct PhotonIllumination {
    /// Boundary beam or explicit uniform layer emission, with independent ID.
    pub source: PhotonSlabSource,
    /// Total incident group-power history, already integrated over beam area.
    pub power: PhotonPowerSchedule,
}

/// Whole-advance source-storage and numerical-energy allowances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhotonHeatingBudget {
    /// Maximum total group-power values across every bound schedule.
    pub max_schedule_values: usize,
    /// Aggregate observed schedule, response and thermal-credit roundoff [J].
    /// This does not bound transport discretization or thermal timestep error.
    pub energy_tolerance_j: f64,
}

/// Separate prescribed source and material heat accounts for one advance.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct ExternalPhotonHeatingAudit {
    /// Integrated incident photon energy [J].
    pub incident_j: f64,
    /// Energy absorbed in the original material cells [J].
    pub absorbed_j: f64,
    /// Uncollided and scattered energy escaping the transport domain [J].
    pub escaped_j: f64,
    /// Energy represented by the mean rates credited to native thermal cells [J].
    pub thermal_credit_j: f64,
    /// Observed power-schedule integration and export rounding [J].
    pub schedule_roundoff_j: f64,
    /// Original transport preparation/application and aggregate export defects [J].
    pub transport_roundoff_j: f64,
    /// Absorbed-energy to mean-power conversion and credit export defects [J].
    pub coupling_roundoff_j: f64,
    /// Signed absorbed + escaped - incident account residual [J].
    pub residual_j: f64,
    /// Number of accepted common source/thermal intervals.
    pub substeps: usize,
}

#[derive(Clone, Debug, Default)]
struct PhotonAccount {
    incident: Sum,
    absorbed: Sum,
    escaped: Sum,
    credit: Sum,
    schedule_roundoff: Sum,
    transport_roundoff: Sum,
    coupling_roundoff: Sum,
    substeps: usize,
}

impl PhotonAccount {
    fn snapshot(&self) -> Result<ExternalPhotonHeatingAudit, FissionThermalError> {
        let incident_j = self.incident.value()?;
        let absorbed_j = self.absorbed.value()?;
        let escaped_j = self.escaped.value()?;
        let residual_j = finite(crate::sum(&[absorbed_j, escaped_j, -incident_j]))?;
        let mut transport = self.transport_roundoff;
        transport.add(self.incident.export_roundoff()?)?;
        transport.add(self.absorbed.export_roundoff()?)?;
        transport.add(self.escaped.export_roundoff()?)?;
        transport.add(residual_j.abs())?;
        let mut coupling = self.coupling_roundoff;
        coupling.add(self.credit.export_roundoff()?)?;
        Ok(ExternalPhotonHeatingAudit {
            incident_j,
            absorbed_j,
            escaped_j,
            residual_j,
            thermal_credit_j: self.credit.value()?,
            schedule_roundoff_j: self.schedule_roundoff.nonnegative_upper_bound()?,
            transport_roundoff_j: transport.nonnegative_upper_bound()?,
            coupling_roundoff_j: coupling.nonnegative_upper_bound()?,
            substeps: self.substeps,
        })
    }
}

#[derive(Clone, Debug)]
struct PhotonPending {
    start_s: f64,
    duration_s: f64,
    end_s: f64,
    account: PhotonAccount,
    audit: ExternalPhotonHeatingAudit,
    transport: Option<PhotonFeedbackCandidate>,
}

/// External photon heating without an artificial nuclear inventory or cell emitter.
/// Immutable response/schedules are shared by trial clones. Accepted state is
/// local, and composes with nuclear or other heat through `CoupledObservers`.
#[derive(Clone, Debug)]
pub struct ExternalPhotonHeating {
    time_s: f64,
    response: Arc<GroupedRadiationResponse>,
    sources: Arc<[(u64, PhotonPowerSchedule)]>,
    receivers: Arc<BTreeMap<CellId, ThermalBinding>>,
    budget: PhotonHeatingBudget,
    account: PhotonAccount,
    audit: ExternalPhotonHeatingAudit,
    pending: Option<PhotonPending>,
    feedback: Option<PhotonTemperatureFeedback>,
}

impl ExternalPhotonHeating {
    /// Bind an explicit fixed response to every original receiver and source.
    /// Spectrum kind/grid must agree exactly. IDs need not name thermal cells.
    /// All response destinations, including zero-heat cells, remain bound.
    pub fn new(
        thermal: &DepositionSimulation<'_>,
        response: Arc<GroupedRadiationResponse>,
        mut sources: Vec<(u64, PhotonPowerSchedule)>,
        receivers: Vec<RadiationReceiver>,
        budget: PhotonHeatingBudget,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<Self, FissionThermalError> {
        checkpoint()?;
        admit_schedules(&sources, budget, &mut checkpoint)?;
        if receivers.is_empty() {
            return Err(FissionThermalError::Invalid(
                "photon heating needs original material receivers",
            ));
        }
        sources.sort_by_key(|source| source.0);
        if !response.source_ids().eq(sources.iter().map(|s| s.0)) {
            return Err(FissionThermalError::Invalid(
                "photon source coverage must be exact and unique",
            ));
        }
        admit_spectra(
            &sources,
            response.radiation(),
            response.group_edges_j(),
            &mut checkpoint,
        )?;
        let receivers =
            bind_receivers(thermal, receivers, response.receiver_ids(), &mut checkpoint)?;
        checkpoint()?;
        Ok(Self {
            time_s: thermal.time_s(),
            response,
            sources: sources.into(),
            receivers: Arc::new(receivers),
            budget,
            account: PhotonAccount::default(),
            audit: ExternalPhotonHeatingAudit::default(),
            pending: None,
            feedback: None,
        })
    }

    /// Derive boundary/internal-source transport from original material cells.
    /// Native mass/volume set density and volume/area sets layer depth. Transport
    /// is prepared once; no thermal emissivity supplies photon absorption.
    #[allow(clippy::too_many_arguments)]
    pub fn from_slab(
        thermal: &DepositionSimulation<'_>,
        data_id: String,
        geometry: PhotonSlabGeometry,
        layers: Vec<PhotonSlabLayer>,
        sources: Vec<PhotonIllumination>,
        ordinates: Vec<PhotonOrdinate>,
        limits: PhotonSlabLimits,
        budget: PhotonHeatingBudget,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, PhotonSlabResponse), FissionThermalError> {
        checkpoint()?;
        if sources.len() > limits.response.max_sources {
            return Err(RadiationError::Budget("photon source count").into());
        }
        let (geometries, powers) = split_sources(sources, budget, &mut checkpoint)?;
        let input = photon_material_input(
            thermal,
            geometry,
            layers,
            powers.iter().map(|s| s.0).collect(),
            false,
            ordinates.len(),
            limits,
            &mut checkpoint,
        )?;
        admit_spectra(&powers, input.radiation, &input.edges_j, &mut checkpoint)?;
        let prepared = photon_call(&mut checkpoint, |poll| {
            prepare_slab_response_with_sources(
                data_id,
                input.radiation,
                input.edges_j,
                input.cells,
                geometries,
                ordinates,
                limits,
                poll,
            )
        })?;
        let heating = Self::new(
            thermal,
            Arc::clone(&prepared.response),
            powers,
            input.receivers,
            budget,
            checkpoint,
        )?;
        Ok((heating, prepared))
    }

    /// Refine photon deposition while preserving original thermal receiver cells.
    /// Beam energy remains attached to its original source through subdivision.
    #[allow(clippy::too_many_arguments)]
    pub fn from_refined_slab(
        thermal: &DepositionSimulation<'_>,
        data_id: String,
        geometry: PhotonSlabGeometry,
        layers: Vec<PhotonSlabLayer>,
        sources: Vec<PhotonIllumination>,
        ordinates: Vec<PhotonOrdinate>,
        limits: PhotonSlabLimits,
        refinement: PhotonSlabRefinement,
        budget: PhotonHeatingBudget,
        mut checkpoint: impl FnMut() -> Result<(), FissionError>,
    ) -> Result<(Self, RefinedPhotonSlabResponse), FissionThermalError> {
        checkpoint()?;
        if sources.len() > limits.response.max_sources {
            return Err(RadiationError::Budget("photon source count").into());
        }
        let (geometries, powers) = split_sources(sources, budget, &mut checkpoint)?;
        let input = photon_material_input(
            thermal,
            geometry,
            layers,
            powers.iter().map(|s| s.0).collect(),
            false,
            ordinates.len(),
            limits,
            &mut checkpoint,
        )?;
        admit_spectra(&powers, input.radiation, &input.edges_j, &mut checkpoint)?;
        let prepared = photon_call(&mut checkpoint, |poll| {
            refine_slab_response_with_sources(
                data_id,
                input.radiation,
                input.edges_j,
                input.cells,
                geometries,
                ordinates,
                limits,
                refinement,
                poll,
            )
        })?;
        let heating = Self::new(
            thermal,
            Arc::clone(&prepared.transport.response),
            powers,
            input.receivers,
            budget,
            checkpoint,
        )?;
        Ok((heating, prepared))
    }

    /// Accepted source/thermal clock [s].
    pub fn time_s(&self) -> f64 {
        self.time_s
    }
    /// Accepted account for the last advance, separate from other heat sources.
    pub fn audit(&self) -> ExternalPhotonHeatingAudit {
        self.audit
    }
    /// Most recently applied immutable transport response. For temperature
    /// feedback, `prepared_temperatures_k` identifies its native start state;
    /// the endpoint is validated but its unused response is prepared lazily.
    pub fn response(&self) -> &Arc<GroupedRadiationResponse> {
        &self.response
    }
    /// Independently named original source schedules in canonical source order.
    pub fn sources(&self) -> &[(u64, PhotonPowerSchedule)] {
        &self.sources
    }

    /// Advance incident histories and original thermal storage atomically.
    pub fn advance_to<'thermal>(
        &self,
        thermal: &DepositionSimulation<'thermal>,
        target_s: f64,
        config: StepConfig,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'thermal>, Self, EnergyAudit), FissionThermalError> {
        self.validate(thermal, self.time_s, &mut checkpoint)?;
        let result = thermal.advance_observed(self, target_s, config, &mut checkpoint)?;
        result
            .1
            .validate(&result.0, result.1.time_s, &mut checkpoint)?;
        result.1.admit_audit(result.1.audit)?;
        Ok(result)
    }

    fn validate(
        &self,
        thermal: &DepositionSimulation<'_>,
        expected_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), FissionThermalError> {
        checkpoint()?;
        if thermal.time_s() != expected_s {
            return Err(FissionThermalError::ClockMismatch {
                thermal_time_s: thermal.time_s(),
                fission_time_s: expected_s,
            });
        }
        for (&id, binding) in self.receivers.iter() {
            checkpoint()?;
            binding.validate(thermal, id)?;
            thermal.temperature_k(id)?;
        }
        if let Some(feedback) = &self.feedback {
            feedback.validate_temperature_support(thermal, checkpoint)?;
        }
        Ok(())
    }

    fn admit_audit(&self, audit: ExternalPhotonHeatingAudit) -> Result<(), FissionThermalError> {
        let mut defects = Sum::default();
        for error in [
            audit.schedule_roundoff_j,
            audit.transport_roundoff_j,
            audit.coupling_roundoff_j,
        ] {
            defects.add(error)?;
        }
        let defect_j = defects.nonnegative_upper_bound()?;
        if defect_j > self.budget.energy_tolerance_j {
            return Err(FissionThermalError::EnergyBudget {
                defect_j,
                tolerance_j: self.budget.energy_tolerance_j,
            });
        }
        Ok(())
    }
}

impl DepositionObserver for ExternalPhotonHeating {
    type Error = FissionThermalError;
    fn begin(&mut self) -> Result<(), Self::Error> {
        if self.pending.is_some() {
            return Err(FissionThermalError::Invalid("uncommitted photon interval"));
        }
        self.account = PhotonAccount::default();
        self.audit = ExternalPhotonHeatingAudit::default();
        if let Some(feedback) = &mut self.feedback {
            feedback.begin();
        }
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
        self.validate(thermal, self.time_s, checkpoint)?;
        let mut next: Option<f64> = None;
        for (_, source) in self.sources.iter() {
            checkpoint()?;
            if let Some(time) = source.next_after(self.time_s) {
                next = Some(next.map_or(time, |old| old.min(time)));
            }
        }
        Ok(next)
    }
    fn step_sources(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        start_s: f64,
        duration_s: f64,
        end_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.validate(thermal, self.time_s, checkpoint)?;
        if self.pending.is_some()
            || start_s != self.time_s
            || !duration_s.is_finite()
            || duration_s <= 0.0
            || !end_s.is_finite()
            || end_s <= start_s
        {
            return Err(FissionThermalError::Invalid(
                "invalid photon/thermal interval",
            ));
        }
        let mut account = self.account.clone();
        let mut energies = Vec::with_capacity(self.sources.len());
        for (id, source) in self.sources.iter() {
            checkpoint()?;
            let (energy, roundoff) = source.interval(start_s, duration_s, end_s, checkpoint)?;
            account.schedule_roundoff.add(roundoff)?;
            energies.push((*id, energy));
        }
        let borrowed: Vec<_> = energies
            .iter()
            .map(|(id, values)| (*id, values.as_slice()))
            .collect();
        // A response refresh belongs to this tentative interval. Publishing it
        // before the thermal endpoint is accepted would leak a failed trial's
        // cache and preparation budget into later source calls.
        let illuminated = energies
            .iter()
            .any(|(_, groups)| groups.iter().any(|&energy| energy > 0.0));
        let transport = if illuminated {
            self.feedback
                .as_ref()
                .map(|feedback| feedback.candidate(thermal, checkpoint))
                .transpose()?
                .flatten()
        } else {
            // Dark intervals still validate original material support, but
            // spend no transport solve on photons that were never supplied.
            None
        };
        let response = transport
            .as_ref()
            .map_or(&self.response, |prepared| &prepared.response);
        let mut checkpoint_error = None;
        let deposition = response.apply(&borrowed, self.budget.energy_tolerance_j, || {
            checkpoint().map_err(|error| {
                checkpoint_error = Some(error);
                RadiationError::Cancelled
            })
        });
        if let Some(error) = checkpoint_error {
            return Err(error.into());
        }
        let deposition = deposition?;
        account.incident.add(deposition.emitted_j)?;
        account.absorbed.add(deposition.deposited_j)?;
        account.escaped.add(deposition.escaped_j)?;
        account.transport_roundoff.add(deposition.roundoff_j)?;
        let mut watts = BTreeMap::new();
        for (id, energy_j) in deposition.receiver_energy_j {
            checkpoint()?;
            let power = finite(energy_j / duration_s)?;
            if energy_j > 0.0 && power == 0.0 {
                return Err(FissionThermalError::Arithmetic);
            }
            let credit = finite(power * duration_s)?;
            account.credit.add(credit)?;
            account.coupling_roundoff.add((credit - energy_j).abs())?;
            watts.insert(CellId(id), power);
        }
        account.substeps = account
            .substeps
            .checked_add(1)
            .ok_or(FissionThermalError::Arithmetic)?;
        let audit = account.snapshot()?;
        self.admit_audit(audit)?;
        checkpoint()?;
        self.pending = Some(PhotonPending {
            start_s,
            duration_s,
            end_s,
            account,
            audit,
            transport,
        });
        Ok(watts)
    }
    fn event(&mut self, _: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        self.validate(thermal, self.time_s, &mut || Ok(()))
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
            .ok_or(FissionThermalError::Invalid("no photon candidate interval"))?;
        if pending.start_s != start_s || pending.duration_s != duration_s {
            return Err(FissionThermalError::Invalid(
                "thermal step changed the photon interval",
            ));
        }
        self.validate(thermal, pending.end_s, checkpoint)?;
        self.admit_audit(pending.audit)?;
        checkpoint()?;
        let pending = self.pending.take().expect("validated photon interval");
        self.time_s = pending.end_s;
        self.account = pending.account;
        self.audit = pending.audit;
        if let Some(transport) = pending.transport {
            self.response = transport.response;
            self.feedback = Some(transport.feedback);
        }
        Ok(())
    }
}

fn admit_spectra(
    sources: &[(u64, PhotonPowerSchedule)],
    radiation: DecayRadiation,
    edges_j: &[f64],
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<(), FissionThermalError> {
    for (_, source) in sources {
        checkpoint()?;
        if source.radiation != radiation || source.group_edges_j() != edges_j {
            return Err(FissionThermalError::Invalid(
                "photon source and response spectra disagree",
            ));
        }
    }
    Ok(())
}

fn admit_schedules(
    sources: &[(u64, PhotonPowerSchedule)],
    budget: PhotonHeatingBudget,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<(), FissionThermalError> {
    if sources.is_empty() || !nonnegative(budget.energy_tolerance_j) {
        return Err(FissionThermalError::Invalid(
            "photon sources and a finite nonnegative allowance required",
        ));
    }
    let mut values = 0_usize;
    for (index, (_, source)) in sources.iter().enumerate() {
        if index % 64 == 0 {
            checkpoint()?;
        }
        values = values
            .checked_add(source.values)
            .filter(|&n| n <= budget.max_schedule_values && n <= MAX_SCHEDULE_VALUES)
            .ok_or(FissionThermalError::Invalid(
                "photon source schedules exceed aggregate allowance",
            ))?;
    }
    Ok(())
}

fn split_sources(
    sources: Vec<PhotonIllumination>,
    budget: PhotonHeatingBudget,
    checkpoint: &mut impl FnMut() -> Result<(), FissionError>,
) -> Result<(Vec<PhotonSlabSource>, Vec<(u64, PhotonPowerSchedule)>), FissionThermalError> {
    checkpoint()?;
    // Each schedule has at least one value; bound allocation before splitting.
    if sources.is_empty()
        || sources.len() > budget.max_schedule_values
        || sources.len() > MAX_SCHEDULE_VALUES
    {
        return Err(FissionThermalError::Invalid(
            "photon source count exceeds allowance",
        ));
    }
    let mut geometries = Vec::with_capacity(sources.len());
    let mut powers = Vec::with_capacity(sources.len());
    for source in sources {
        checkpoint()?;
        geometries.push(source.source);
        powers.push((source.source.source_id, source.power));
    }
    admit_schedules(&powers, budget, checkpoint)?;
    Ok((geometries, powers))
}

#[cfg(test)]
mod tests;
