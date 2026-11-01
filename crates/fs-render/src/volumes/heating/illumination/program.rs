//! Finite-duration source power on the existing native optical/thermal clock.

use super::schedule::SpectralPowerSchedule;
use super::{
    CheckedSum, EnergySums, IlluminationBudget, IncidentIlluminationSource, IncidentScatteringBeam,
    OpticalFeedbackError, OpticalHeatingError, ScatteringDepositionConfig,
    ScatteringHeatingTolerance, TemperatureIlluminationHeating, TemperatureOpticalField, invalid,
};
use fs_blake3::DomainHasher;
use fs_geom::{Point3, Vec3};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use std::collections::BTreeMap;
use std::sync::Arc;

/// One prescribed spectral pulse/ramp source or weighted source-quadrature node.
/// Wavelengths come from the shared optical field; powers already include any
/// spatial/angular quadrature weight. The source is dark outside its segments.
#[derive(Clone, Debug)]
pub struct ScheduledIlluminationSource<const N: usize> {
    pub source_id: u64,
    pub origin: Point3,
    pub unit_direction: Vec3,
    pub power: SpectralPowerSchedule<N>,
}

pub(super) struct Program<const N: usize> {
    profiles: Vec<SpectralPowerSchedule<N>>,
    pub(super) transitions_s: Vec<f64>,
}

#[derive(Clone)]
struct Pending<const N: usize> {
    start_s: f64,
    duration_s: f64,
    end_s: f64,
    sources_w: Arc<BTreeMap<CellId, f64>>,
    sums: EnergySums<N>,
}

#[derive(Clone)]
pub(super) struct ProgramState<const N: usize> {
    mean_power_w: Vec<[f64; N]>,
    accepted_sources_w: Option<Arc<BTreeMap<CellId, f64>>>,
    pending: Option<Pending<N>>,
}

impl<const N: usize> ProgramState<N> {
    pub(super) fn clear_pending(&mut self) {
        self.pending = None;
    }
    pub(super) fn sources_w(&self) -> Option<&BTreeMap<CellId, f64>> {
        self.pending
            .as_ref()
            .map(|pending| pending.sources_w.as_ref())
            .or(self.accepted_sources_w.as_deref())
    }
}

impl<'g, const N: usize> TemperatureIlluminationHeating<'g, N> {
    /// Admit bounded, piecewise-linear spectral source powers without tracing.
    /// Native steps end at every source transition and use its exact interval
    /// mean power. Optical coefficients stay frozen over each native step;
    /// changing material state still requires temporal convergence assessment.
    ///
    /// Transport and uncertainty admission use each lane's maximum power over
    /// the complete schedule. These simultaneous envelopes conservatively bound
    /// every actual profile, even when different sources peak at different times.
    /// Power-only changes reuse that response without reseeding or retracing.
    #[allow(clippy::too_many_arguments)]
    pub fn from_scheduled_field(
        field: Arc<TemperatureOpticalField<'g, N>>,
        sources: &[ScheduledIlluminationSource<N>],
        thermal: &DepositionSimulation<'_>,
        budget: IlluminationBudget,
        config: ScatteringDepositionConfig,
        tolerance: ScatteringHeatingTolerance,
        max_total_power_segments: usize,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, OpticalFeedbackError> {
        checkpoint()?;
        if sources.len() > budget.max_sources {
            return Err(OpticalFeedbackError::Budget);
        }
        let mut segments = 0_usize;
        for source in sources {
            checkpoint()?;
            segments = segments
                .checked_add(source.power.segment_count())
                .ok_or(OpticalFeedbackError::Budget)?;
            if segments > max_total_power_segments {
                return Err(OpticalFeedbackError::Budget);
            }
        }
        let transition_capacity = segments
            .checked_mul(2)
            .ok_or(OpticalFeedbackError::Budget)?;
        let mut ordered = Vec::new();
        ordered
            .try_reserve_exact(sources.len())
            .map_err(|_| OpticalFeedbackError::Budget)?;
        ordered.extend(sources.iter());
        ordered.sort_unstable_by_key(|source| source.source_id);
        let mut envelopes = Vec::new();
        let mut profiles = Vec::new();
        let mut transitions_s = Vec::new();
        let mut mean_power_w = Vec::new();
        envelopes
            .try_reserve_exact(sources.len())
            .map_err(|_| OpticalFeedbackError::Budget)?;
        profiles
            .try_reserve_exact(sources.len())
            .map_err(|_| OpticalFeedbackError::Budget)?;
        mean_power_w
            .try_reserve_exact(sources.len())
            .map_err(|_| OpticalFeedbackError::Budget)?;
        mean_power_w.resize(sources.len(), [0.0; N]);
        transitions_s
            .try_reserve_exact(transition_capacity)
            .map_err(|_| OpticalFeedbackError::Budget)?;
        for source in ordered {
            checkpoint()?;
            envelopes.push(IncidentIlluminationSource {
                source_id: source.source_id,
                beam: IncidentScatteringBeam {
                    origin: source.origin,
                    unit_direction: source.unit_direction,
                    wavelengths_nm: *field.wavelengths_nm(),
                    incident_power_w: source.power.peak_power_w(),
                },
            });
            for segment in source.power.segments() {
                checkpoint()?;
                transitions_s.extend([segment.start_s, segment.end_s]);
            }
            profiles.push(source.power.clone());
        }
        transitions_s.sort_unstable_by(f64::total_cmp);
        transitions_s.dedup();
        let mut result = Self::from_field(
            field,
            &envelopes,
            thermal,
            budget,
            config,
            tolerance,
            &mut checkpoint,
        )?;
        let model = Arc::get_mut(&mut result.model).expect("new illumination model is unshared");
        let mut identity = DomainHasher::new("org.frankensim.render.scheduled-illumination.v1");
        identity.update(model.identity.as_bytes());
        for profile in &profiles {
            checkpoint()?;
            identity.update(profile.identity().as_bytes());
        }
        model.identity = identity.finalize();
        model.program = Some(Program {
            profiles,
            transitions_s,
        });
        result.program_state = Some(ProgramState {
            mean_power_w,
            accepted_sources_w: None,
            pending: None,
        });
        checkpoint()?;
        Ok(result)
    }

    pub(super) fn prepare_scheduled_interval(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        start_s: f64,
        duration_s: f64,
        end_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, OpticalFeedbackError> {
        self.check_clock(thermal)?;
        checkpoint()?;
        if start_s != self.time_s
            || !duration_s.is_finite()
            || duration_s <= 0.0
            || !end_s.is_finite()
            || end_s <= start_s
            || (end_s != start_s + duration_s && duration_s != end_s - start_s)
        {
            return Err(invalid("matching scheduled illumination interval required"));
        }
        if self.evaluations >= self.budget.max_source_evaluations {
            return Err(OpticalFeedbackError::Budget);
        }
        let program = self
            .model
            .program
            .as_ref()
            .expect("scheduled illumination model");
        let state = self
            .program_state
            .as_mut()
            .expect("scheduled illumination state");
        state.clear_pending();
        let mut illuminated = false;
        for (retained, profile) in state.mean_power_w.iter_mut().zip(&program.profiles) {
            let mean = profile.mean_power_w(start_s, duration_s, end_s, &mut *checkpoint)?;
            illuminated |= mean.iter().any(|power| *power > 0.0);
            *retained = mean;
        }
        if illuminated {
            self.prepare_transport(thermal, checkpoint)?;
        } else {
            // Dark intervals preserve the last transport for later reuse, but
            // still enforce the actual native material/source support domain.
            self.model
                .field
                .validate_temperatures(thermal, checkpoint)?;
        }
        let mut sums = self.sums.clone();
        sums.add_duration(duration_s)?;
        let mut receiving = BTreeMap::<CellId, CheckedSum>::new();
        let mut absorbed = CheckedSum::default();
        if illuminated {
            let prepared = self.current.as_ref().expect("prepared illumination");
            let state = self
                .program_state
                .as_ref()
                .expect("scheduled illumination state");
            for ((heating, source), mean) in prepared
                .sources
                .iter()
                .zip(&self.model.sources)
                .zip(&state.mean_power_w)
            {
                checkpoint()?;
                let mut factors = [0.0; N];
                for lane in 0..N {
                    let peak = source.beam.incident_power_w[lane];
                    let power = mean[lane];
                    if power > peak || (power > 0.0 && peak == 0.0) {
                        return Err(invalid("scheduled power exceeds its admitted envelope"));
                    }
                    if power > 0.0 {
                        factors[lane] = power / peak;
                        if factors[lane] == 0.0 {
                            return Err(invalid("scheduled power ratio is not representable"));
                        }
                    }
                    absorbed.add(heating.deposition().absorbed.mean_w[lane] * factors[lane])?;
                }
                for cell in &heating.deposition().cells {
                    checkpoint()?;
                    let receiver = self
                        .model
                        .field
                        .receiver(cell.cell_index)
                        .ok_or(OpticalHeatingError::UnmappedCell(cell.cell_index))?;
                    let mut heat = CheckedSum::default();
                    for lane in 0..N {
                        heat.add(cell.absorbed.mean_w[lane] * factors[lane])?;
                    }
                    let power = heat.value()?;
                    if power > 0.0 {
                        receiving.entry(receiver).or_default().add(power)?;
                    }
                }
                sums.add_source_radiation(&heating.scaled_radiation_energy(duration_s, &factors)?)?;
            }
        }
        let mut sources_w = BTreeMap::new();
        let mut thermal_input = CheckedSum::default();
        for (receiver, power) in receiving {
            checkpoint()?;
            let power = power.value()?;
            thermal_input.add(power)?;
            sources_w.insert(receiver, power);
        }
        let thermal_w = thermal_input.value()?;
        if (thermal_w - absorbed.value()?).abs() > self.model.tolerance.power_roundoff_w {
            return Err(invalid(
                "scheduled cell heat disagrees with scaled absorption",
            ));
        }
        sums.add_thermal_input(thermal_w * duration_s)?;
        let result = sources_w.clone();
        checkpoint()?;
        self.program_state
            .as_mut()
            .expect("scheduled illumination state")
            .pending = Some(Pending {
            start_s,
            duration_s,
            end_s,
            sources_w: Arc::new(sources_w),
            sums,
        });
        self.evaluations += 1;
        Ok(result)
    }

    pub(super) fn accept_scheduled_interval(
        &mut self,
        start_s: f64,
        duration_s: f64,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), OpticalFeedbackError> {
        let state = self
            .program_state
            .as_mut()
            .expect("scheduled illumination state");
        let pending = state
            .pending
            .as_ref()
            .ok_or_else(|| invalid("native step needs prepared scheduled illumination"))?;
        if start_s != self.time_s
            || pending.start_s != start_s
            || pending.duration_s != duration_s
            || pending.end_s != thermal.time_s()
        {
            return Err(invalid(
                "accepted native interval differs from scheduled illumination",
            ));
        }
        self.model
            .field
            .validate_temperatures(thermal, checkpoint)?;
        let energy = pending.sums.value()?;
        checkpoint()?;
        let pending = state.pending.take().expect("validated scheduled interval");
        self.sums = pending.sums;
        self.energy = energy;
        self.time_s = thermal.time_s();
        state.accepted_sources_w = Some(pending.sources_w);
        Ok(())
    }
}
