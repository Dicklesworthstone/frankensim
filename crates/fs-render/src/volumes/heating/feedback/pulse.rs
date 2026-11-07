//! Prescribed pulse/ramp illumination on the original nonlinear thermal owner.
//!
//! Source transitions clip native intervals. Linear and two-photon optical
//! coefficients are frozen at the actual native start state; the complete
//! affine power waveform is integrated analytically through the cached cell
//! path. A nonlinear source is never evaluated at a mean intensity. This
//! removes source-waveform quadrature error, not coupled temperature-step error.

use super::*;
use crate::volumes::heating::illumination::schedule::{
    SpectralPowerSchedule, SpectralPowerSegment,
};

mod moments;

#[derive(Clone, Debug)]
struct Pending<const N: usize> {
    start_s: f64,
    duration_s: f64,
    end_s: f64,
    end_power_w: [f64; N],
    power: OpticalPowerBudget<N>,
    sources_w: BTreeMap<CellId, f64>,
    sums: EnergySums<N>,
    energy: OpticalEnergyBudget<N>,
}

#[derive(Clone, Debug)]
pub(super) struct PulseState<const N: usize> {
    profile: SpectralPowerSchedule<N>,
    pending: Option<Pending<N>>,
}

struct Interval<const N: usize> {
    start_w: [f64; N],
    end_w: [f64; N],
    mean_w: [f64; N],
}

impl<const N: usize> TemperatureOpticalHeating<N> {
    /// Bind prescribed affine ramps, pulse jumps and dark gaps before this
    /// observer's first advance. The existing `IncidentBeam` powers are the
    /// per-lane admission ceilings; every schedule endpoint must fit them.
    /// Geometry, wavelength, density and uniform beam area stay prescribed.
    ///
    /// The original `SpectralPowerSchedule` owns the waveform and transitions.
    /// Native intervals stop at every transition and receive the analytically
    /// integrated nonlinear cell heat, including downstream shielding. Both
    /// sides of endpoint jumps obey the actual local irradiance/temperature
    /// domains. The original source-evaluation limit covers each preparation
    /// and each nonlinear endpoint-domain pass.
    ///
    /// Use the existing `advance_refined` to assess temperature feedback error.
    /// Prompt two-photon thermalization still excludes carrier storage,
    /// saturation, nonlinear refraction, diffraction and ultrafast relaxation.
    pub fn with_power_schedule(
        mut self,
        thermal: &DepositionSimulation<'_>,
        profile: SpectralPowerSchedule<N>,
        max_power_segments: usize,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, OpticalFeedbackError> {
        checkpoint()?;
        self.check_clock(thermal)?;
        self.check_pulse_pending()?;
        if self.energy.duration_s != 0.0 {
            return Err(invalid(
                "bind the optical power schedule before the first advance",
            ));
        }
        if profile.segment_count() > max_power_segments {
            return Err(OpticalFeedbackError::Budget);
        }
        if profile
            .peak_power_w()
            .iter()
            .zip(self.model.beam.incident_power_w)
            .any(|(&peak, ceiling)| peak > ceiling)
        {
            return Err(invalid(
                "optical pulse exceeds the original incident-power ceiling",
            ));
        }
        let present = power_at(&profile, self.time_s)?;
        let (power, sources_w) = self
            .model
            .evaluate_at_power(thermal, present, &mut checkpoint)?;
        checkpoint()?;
        self.power = power;
        self.sources_w = sources_w;
        self.pulse = Some(PulseState {
            profile,
            pending: None,
        });
        Ok(self)
    }

    pub(super) fn check_pulse_pending(&self) -> Result<(), OpticalFeedbackError> {
        if self
            .pulse
            .as_ref()
            .is_some_and(|state| state.pending.is_some())
        {
            Err(invalid("uncommitted optical pulse interval"))
        } else {
            Ok(())
        }
    }

    pub(super) fn next_pulse_transition(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Option<f64>, OpticalFeedbackError> {
        let Some(state) = &self.pulse else {
            return Ok(None);
        };
        checkpoint()?;
        self.check_clock(thermal)?;
        state.profile.next_transition_time_s(self.time_s)
    }

    pub(super) fn prepare_pulse_interval(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        start_s: f64,
        duration_s: f64,
        end_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, OpticalFeedbackError> {
        checkpoint()?;
        self.check_clock(thermal)?;
        self.check_pulse_pending()?;
        if start_s != self.time_s {
            return Err(invalid(
                "optical pulse starts at a different accepted clock",
            ));
        }
        self.admit_evaluation()?;
        let state = self.pulse.as_ref().expect("scheduled optical source");
        let interval = interval(&state.profile, start_s, duration_s, end_s, checkpoint)?;
        let (power, sources_w) =
            self.model
                .evaluate_pulse(thermal, &interval, state.profile.identity(), checkpoint)?;
        let mut sums = self.sums.clone();
        sums.add(&power, duration_s)?;
        let energy = sums.value()?;
        checkpoint()?;
        let result = sources_w.clone();
        self.pulse
            .as_mut()
            .expect("scheduled optical source")
            .pending = Some(Pending {
            start_s,
            duration_s,
            end_s,
            end_power_w: interval.end_w,
            power,
            sources_w,
            sums,
            energy,
        });
        Ok(result)
    }

    pub(super) fn accept_pulse_interval(
        &mut self,
        start_s: f64,
        duration_s: f64,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), OpticalFeedbackError> {
        let state = self.pulse.as_ref().expect("scheduled optical source");
        let pending = state
            .pending
            .as_ref()
            .ok_or_else(|| invalid("native step needs a prepared optical pulse interval"))?;
        if start_s != self.time_s
            || pending.start_s != start_s
            || pending.duration_s != duration_s
            || pending.end_s != thermal.time_s()
        {
            return Err(invalid(
                "accepted native interval differs from prepared optical pulse",
            ));
        }
        let left = pending.end_power_w;
        let right = power_at(&state.profile, pending.end_s)?;
        if self.model.two_photon.is_some() {
            self.admit_evaluation()?;
        }
        self.model
            .check_domains_at_power(thermal, left, checkpoint)?;
        if self.model.two_photon.is_some() && left != right {
            self.admit_evaluation()?;
            self.model
                .check_domains_at_power(thermal, right, checkpoint)?;
        }
        checkpoint()?;
        let pending = self
            .pulse
            .as_mut()
            .expect("scheduled optical source")
            .pending
            .take()
            .expect("validated optical pulse interval");
        self.power = pending.power;
        self.sources_w = pending.sources_w;
        self.sums = pending.sums;
        self.energy = pending.energy;
        self.time_s = pending.end_s;
        Ok(())
    }

    pub(super) fn check_pulse_event(
        &mut self,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), OpticalFeedbackError> {
        self.check_pulse_pending()?;
        let state = self.pulse.as_ref().expect("scheduled optical source");
        let power = power_at(&state.profile, self.time_s)?;
        if self.model.two_photon.is_some() {
            self.admit_evaluation()?;
        }
        self.model
            .check_domains_at_power(thermal, power, &mut || Ok(()))
    }
}

impl<const N: usize> FeedbackModel<N> {
    fn evaluate_pulse(
        &self,
        thermal: &DepositionSimulation<'_>,
        interval: &Interval<N>,
        profile_identity: ContentHash,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(OpticalPowerBudget<N>, BTreeMap<CellId, f64>), OpticalFeedbackError> {
        // A frozen linear transport is linear in source power. Its exact mean
        // uses the original kernel once, with no nonlinear moment work.
        if self.two_photon.is_none() {
            return self.evaluate_at_power(thermal, interval.mean_w, checkpoint);
        }
        // Preserve the exact original cell arithmetic for rectangular levels.
        if interval.start_w == interval.end_w {
            return self.evaluate_at_power(thermal, interval.start_w, checkpoint);
        }
        let mut hash = DomainHasher::new("org.frankensim.render.temperature-pulse-state.v1");
        hash.update(self.identity.as_bytes());
        hash.update(profile_identity.as_bytes());
        for power in interval.start_w.into_iter().chain(interval.end_w) {
            hash.update(&power.to_bits().to_le_bytes());
        }
        let mut peak: [f64; N] =
            std::array::from_fn(|lane| interval.start_w[lane].max(interval.end_w[lane]));
        let fraction = std::array::from_fn::<_, N, _>(|lane| {
            if peak[lane] == 0.0 {
                1.0
            } else {
                interval.start_w[lane].min(interval.end_w[lane]) / peak[lane]
            }
        });
        let mut shape = [1.0; N];
        let mut transmitted = interval.mean_w;
        let mut absorbed = [CheckedSum::default(); N];
        let mut sources = BTreeMap::<CellId, CheckedSum>::new();
        for (index, segment) in self.path.iter().enumerate() {
            checkpoint()?;
            let temperature = thermal.temperature_k(segment.thermal_cell)?;
            let (coefficients, enthalpy) =
                self.tables[segment.table].at(thermal, segment, temperature)?;
            hash.update(&temperature.to_bits().to_le_bytes());
            if let Some(h) = enthalpy {
                hash.update(&h.to_bits().to_le_bytes());
            }
            let mut heat = CheckedSum::default();
            for lane in 0..N {
                let linear = absorb_beam(
                    peak[lane],
                    segment.path.density_kg_m3,
                    segment.path.length_m,
                    coefficients[lane],
                );
                // Passive cell transport is monotone in its incident power:
                // the local peak admits every interior point of this ramp.
                let full = if let Some(nonlinear) = &self.two_photon {
                    nonlinear.absorb(index, segment, temperature, peak[lane], coefficients[lane])?
                } else {
                    linear
                };
                let mean = moments::mean_cell(
                    fraction[lane],
                    shape[lane],
                    full.0,
                    full.1,
                    linear.0,
                    linear.1,
                )
                .map_err(invalid)?;
                absorbed[lane].add(mean.absorbed_w)?;
                heat.add(mean.absorbed_w)?;
                transmitted[lane] = mean.transmitted_w;
                shape[lane] = mean.next_shape;
                peak[lane] = full.1;
            }
            sources
                .entry(segment.thermal_cell)
                .or_default()
                .add(heat.value()?)?;
        }
        self.finish_power(
            hash.finalize(),
            interval.mean_w,
            transmitted,
            absorbed,
            sources,
            checkpoint,
        )
    }
}

fn interval<const N: usize>(
    profile: &SpectralPowerSchedule<N>,
    start_s: f64,
    duration_s: f64,
    end_s: f64,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<Interval<N>, OpticalFeedbackError> {
    if profile
        .next_transition_time_s(start_s)?
        .is_some_and(|time| time < end_s)
    {
        return Err(invalid(
            "native optical interval crosses a pulse transition",
        ));
    }
    let mean_w = profile.mean_power_w(start_s, duration_s, end_s, checkpoint)?;
    let segments = profile.segments();
    let index = segments.partition_point(|s| s.end_s <= start_s);
    let Some(segment) = segments.get(index).filter(|s| s.start_s <= start_s) else {
        return Ok(Interval {
            start_w: [0.0; N],
            end_w: [0.0; N],
            mean_w,
        });
    };
    Ok(Interval {
        start_w: segment_power(segment, start_s),
        end_w: segment_power(segment, end_s),
        mean_w,
    })
}

fn power_at<const N: usize>(
    profile: &SpectralPowerSchedule<N>,
    time_s: f64,
) -> Result<[f64; N], OpticalFeedbackError> {
    if !nonnegative(time_s) {
        return Err(invalid("finite nonnegative pulse time required"));
    }
    let segments = profile.segments();
    let index = segments.partition_point(|s| s.end_s <= time_s);
    Ok(segments
        .get(index)
        .filter(|s| s.start_s <= time_s)
        .map_or([0.0; N], |s| segment_power(s, time_s)))
}

fn segment_power<const N: usize>(segment: &SpectralPowerSegment<N>, time_s: f64) -> [f64; N] {
    let fraction = (time_s - segment.start_s) / (segment.end_s - segment.start_s);
    std::array::from_fn(|lane| {
        let a = segment.start_power_w[lane];
        let b = segment.end_power_w[lane];
        if fraction == 0.0 {
            a
        } else if fraction == 1.0 {
            b
        } else if a <= b {
            a + (b - a) * fraction
        } else {
            b + (a - b) * (1.0 - fraction)
        }
    })
}
