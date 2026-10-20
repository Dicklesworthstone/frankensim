//! Physical pressure and jet-direction assignments on one mechanical clock.
//! Gesture interpolation remains in fs-scenario. The existing pressure compiler
//! handles pressure; the optional angle track samples the same ramp owner.
use super::{FlueOutletConfig, FlueVoice, invalid, source_error};
use crate::flue_loop::FlueControl;
use crate::render::ControlDelta;
use crate::render::schedule::{
    GestureCompileError, PressureGestureBinding, compile_pressure_gestures,
};
use fs_scenario::gesture::{GestureSchedule, GestureTarget};

#[derive(Clone, Copy, Debug)]
pub(super) struct ScheduledFlueControl {
    pub sample: u64,
    pub delta: FlueControl,
}

pub(super) fn compile(
    voice: &FlueVoice,
    schedule: &GestureSchedule,
    config: &FlueOutletConfig,
    mechanics_rate_hz: u32,
    mechanical_samples: u64,
) -> Result<Vec<ScheduledFlueControl>, GestureCompileError> {
    let tracks = schedule.tracks();
    if !(1..=2).contains(&tracks.len()) {
        return Err(invalid("flue needs one blowing-pressure track and at most one jet-angle track"));
    }
    let mut pressure = None;
    let mut angle = None;
    for track in tracks {
        match track.target {
            GestureTarget::BlowingPressure if pressure.is_none() => pressure = Some(track),
            GestureTarget::JetAngle if angle.is_none() => angle = Some(track),
            _ => return Err(invalid("flue tracks must address pressure and optional jet angle exactly once")),
        }
    }
    let pressure = pressure.ok_or_else(|| invalid("flue needs an explicitly typed blowing-pressure track"))?;
    let rate = schedule.control_rate_hz;
    if rate == 0 || rate > mechanics_rate_hz || mechanical_samples == 0 {
        return Err(invalid("flue controls need a finite window and a positive rate no faster than mechanics"));
    }
    let last_tick = u128::from(mechanical_samples - 1) * u128::from(rate)
        / u128::from(mechanics_rate_hz);
    let ticks = last_tick + 1;
    if ticks * tracks.len() as u128 > config.max_controls as u128
        || tracks.iter().any(|track| track.events.iter()
            .any(|event| event.time_s > last_tick as f64 / f64::from(rate))) {
        return Err(invalid("flue commands exceed the observed window or total control-assignment budget"));
    }
    // Bound the actual source sampler's track/event scans before cloning or
    // emitting any controls. Coalescing a held input does not erase this cost.
    let per_tick: u128 = tracks.iter()
        .map(|track| tracks.len() as u128 + track.events.len() as u128 + 1).sum();
    let required = ticks.checked_mul(per_tick)
        .ok_or_else(|| invalid("flue gesture sampling work overflows"))?;
    if required > u128::from(config.max_compile_work) {
        return Err(GestureCompileError::WorkBudget { required, allowed: config.max_compile_work });
    }
    let pressure_schedule = GestureSchedule::try_new(rate, vec![pressure.clone()])
        .map_err(GestureCompileError::Gesture)?;
    let pressure_controls = compile_pressure_gestures(&pressure_schedule,
        &[PressureGestureBinding { track: pressure.id.clone(), voice: 0 }],
        mechanics_rate_hz, mechanical_samples, config.max_compile_work)?;
    let mut controls = Vec::new();
    controls.try_reserve_exact(pressure_controls.len())
        .map_err(|_| invalid("flue controls exceed available allocation"))?;
    for event in pressure_controls {
        let ControlDelta::SetBlowingPressure { voice: 0, pressure_pa } = event.delta else {
            return Err(invalid("flue pressure compiler produced an incompatible control"));
        };
        voice.validate_blowing_pressure(pressure_pa)
            .map_err(|e| GestureCompileError::Render(source_error(e)))?;
        controls.push(ScheduledFlueControl {
            sample: event.sample, delta: FlueControl::SetBlowingPressure(pressure_pa),
        });
    }
    if let Some(angle) = angle {
        let ticks = u64::try_from(ticks).map_err(|_| invalid("flue control clock exceeds u64"))?;
        let mut previous = None;
        for tick in 0..ticks {
            let value = schedule.sample(&angle.id, tick).map_err(GestureCompileError::Gesture)?;
            voice.validate_jet_angle(value)
                .map_err(|e| GestureCompileError::Render(source_error(e)))?;
            if previous != Some(value.to_bits()) {
                let sample = (u128::from(tick) * u128::from(mechanics_rate_hz))
                    .div_ceil(u128::from(rate));
                let sample = u64::try_from(sample)
                    .map_err(|_| invalid("flue control sample exceeds u64"))?;
                controls.try_reserve(1)
                    .map_err(|_| invalid("flue controls exceed available allocation"))?;
                controls.push(ScheduledFlueControl { sample, delta: FlueControl::SetJetAngle(value) });
                previous = Some(value.to_bits());
            }
        }
        // Stable ordering: pressure precedes angle at the same physical sample.
        // Both are applied before the source step; neither advances a history.
        controls.sort_by_key(|event| event.sample);
    }
    Ok(controls)
}
