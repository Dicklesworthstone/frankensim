//! Source-mesh pressure amplitudes through the existing pressure event selector.
use super::{GestureCompileError, ScheduledControl, ScheduledTwoSidedPressure, invalid};
use crate::bernoulli_aperture::performance::pressure_candidate;
use crate::pcm_wav::observation::PressureRenderer;
use crate::render::{ControlDelta, RenderError};

/// Assignment to one blade's previously bound dimensionless pressure trace.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialPressureEvent {
    /// Apply before this MECHANICAL sample, within the finite render window.
    pub sample: u64,
    /// Signed reference pressure [Pa]. Local pressure is psi(x)*pressure_pa.
    /// Zero releases this load without resetting material or acoustic state.
    pub pressure_pa: f64,
}

pub(super) struct PressureSchedule {
    pub(super) events: Vec<MaterialPressureEvent>,
    controls: Vec<ScheduledControl>,
    next: usize,
    held: f64,
    pub(super) work: u64,
}
impl PressureSchedule {
    fn compile(mut events: Vec<MaterialPressureEvent>, samples: u64,
        max_events: usize, max_work: u64) -> Result<Self, GestureCompileError>
    {
        let work = 2 * events.len() as u128;
        if events.len() > max_events || work > u128::from(max_work)
            || events.iter().any(|e| e.sample >= samples || !e.pressure_pa.is_finite()) {
            return Err(invalid("spatial pressure assignments require finite pascals, valid samples and shared capacity"));
        }
        // The shared pressure selector already implements stable same-sample
        // last-write semantics and immutable candidate evaluation.
        events.sort_by_key(|e| e.sample);
        let controls = events.iter().map(|e| ScheduledControl {
            sample: e.sample,
            delta: ControlDelta::SetBlowingPressure { voice: 0, pressure_pa: e.pressure_pa },
        }).collect();
        Ok(Self { events, controls, next: 0, held: 0.0, work: work as u64 })
    }
    pub(super) fn candidate(&self, sample: u64) -> f64 {
        pressure_candidate(&self.controls, self.next, self.held, sample).1
    }
    pub(super) fn accept(&mut self, sample: u64) {
        (self.next, self.held) = pressure_candidate(&self.controls, self.next, self.held, sample);
    }
}

impl ScheduledTwoSidedPressure {
    /// Schedule an additional spatial pressure source on one original material
    /// blade. Its shape/overlap is already bound to the specimen; these are Pa,
    /// independent of N-valued force inputs and the ideal series pressure source.
    /// Equal-time events use stable last-write order on the mechanical clock.
    /// Attachment and accepted controls share the existing whole-phrase budgets.
    pub fn with_material_pressure_controls(mut self, blade: usize,
        events: Vec<MaterialPressureEvent>) -> Result<Self, GestureCompileError>
    {
        if blade >= 2 || self.poisoned || self.source.poisoned
            || self.source.mechanical_samples_rendered() != 0
            || self.spatial_pressures[blade].is_some()
            || self.source.held_spatial_pressures_pa[blade] != 0.0 {
            return Err(invalid("spatial pressure controls require one cold attachment per unloaded material port"));
        }
        let mut probe = [0.0; 2];
        probe[blade] = 1.0;
        self.source.system().valve().validate_pressures(probe)
            .map_err(|e| GestureCompileError::Render(RenderError::Voice(e)))?;
        let (assignments, work) = self.remaining_control_budgets()?;
        let samples = self.config.samples.checked_mul(self.source.mechanics_per_output() as u64)
            .ok_or_else(|| invalid("spatial pressure mechanical window overflow"))?;
        self.spatial_pressures[blade] = Some(PressureSchedule::compile(events, samples, assignments, work)?);
        Ok(self)
    }

    /// Set unscheduled spatial pressure amplitudes between callbacks. A bound
    /// pressure program owns its whole phrase and refuses manual overrides.
    pub fn set_spatial_pressures(&mut self, pressures_pa: [f64; 2]) -> Result<(), RenderError> {
        self.validate_sample_count(0)?;
        if self.spatial_pressures.iter().any(Option::is_some) {
            return Err(RenderError::Control { what: "spatial pressure programs cannot be overridden manually" });
        }
        self.source.set_spatial_pressures(pressures_pa)
    }
    /// Original raw events committed with completely ingested mechanical steps.
    pub fn applied_material_pressure_controls(&self, blade: usize) -> &[MaterialPressureEvent] {
        self.spatial_pressures.get(blade).and_then(Option::as_ref)
            .map_or(&[], |p| &p.events[..p.next])
    }
    /// Original raw events still pending in stable mechanical-time order.
    pub fn pending_material_pressure_controls(&self, blade: usize) -> &[MaterialPressureEvent] {
        self.spatial_pressures.get(blade).and_then(Option::as_ref)
            .map_or(&[], |p| &p.events[p.next..])
    }
    /// Last accepted scheduled amplitudes [Pa], or explicit unscheduled loads.
    pub fn held_material_pressures_pa(&self) -> [f64; 2] {
        std::array::from_fn(|i| self.spatial_pressures[i].as_ref()
            .map_or(self.source.held_spatial_pressures_pa[i], |p| p.held))
    }
}
