//! Physical pressure, jet direction and closing forces on one mechanical clock.
//! All interpolation remains in fs-scenario, including interruptible force ramps.
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
    if !(1..=66).contains(&tracks.len()) {
        return Err(invalid("flue needs bounded pressure, angle and physical boundary-force tracks"));
    }
    let mut pressure = None;
    let mut angle = None;
    let mut forces = Vec::new();
    let source_paths = voice.network().and_then(|network| network.source_control_paths());
    for track in tracks {
        match track.target {
            GestureTarget::BlowingPressure if pressure.is_none() => pressure = Some(track),
            GestureTarget::JetAngle if angle.is_none() => angle = Some(track),
            GestureTarget::Parameter(ref path) if source_paths.is_some_and(|source| path == &source.pressure)
                && pressure.is_none() => pressure = Some(track),
            GestureTarget::Parameter(ref path) if source_paths.is_some_and(|source| path == &source.jet_angle)
                && angle.is_none() => angle = Some(track),
            GestureTarget::Parameter(ref path) => {
                let node = voice.network().and_then(|network| network.boundary_node(path))
                    .ok_or_else(|| invalid("flue parameter track has no bound physical boundary consumer"))?;
                if forces.iter().any(|(_, existing)| *existing == node) {
                    return Err(invalid("one physical boundary can have only one closing-force track"));
                }
                forces.push((track, node));
            }
            _ => return Err(invalid("flue tracks must address pressure, optional jet angle or bound closing force")),
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
    // Lower an already-bound pressure parameter to the original scalar pressure
    // compiler. The public schedule remains all-v1 or all-resolved-v2; neither
    // the grammar nor the shared interpolation owner is changed by this adapter.
    let mut pressure_track = pressure.clone();
    pressure_track.target = GestureTarget::BlowingPressure;
    let pressure_schedule = GestureSchedule::try_new(rate, vec![pressure_track])
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
    let scalar_tracks = angle.into_iter().map(|track| (track, None))
        .chain(forces.into_iter().map(|(track, node)| (track, Some(node))));
    for (track, boundary) in scalar_tracks {
        let ticks = u64::try_from(ticks).map_err(|_| invalid("flue control clock exceeds u64"))?;
        let mut previous = None;
        for tick in 0..ticks {
            let value = schedule.sample(&track.id, tick).map_err(GestureCompileError::Gesture)?;
            let delta = if let Some(node) = boundary {
                voice.validate_boundary_closing_force(node, value)
                    .map_err(|e| GestureCompileError::Render(source_error(e)))?;
                FlueControl::SetBoundaryClosingForce { node, force_n: value }
            } else {
                voice.validate_jet_angle(value)
                    .map_err(|e| GestureCompileError::Render(source_error(e)))?;
                FlueControl::SetJetAngle(value)
            };
            if previous != Some(value.to_bits()) {
                let sample = (u128::from(tick) * u128::from(mechanics_rate_hz))
                    .div_ceil(u128::from(rate));
                let sample = u64::try_from(sample)
                    .map_err(|_| invalid("flue control sample exceeds u64"))?;
                controls.try_reserve(1)
                    .map_err(|_| invalid("flue controls exceed available allocation"))?;
                controls.push(ScheduledFlueControl { sample, delta });
                previous = Some(value.to_bits());
            }
        }
    }
    // Stable ordering: all assignments precede the same physical source step;
    // a force/angle change alone does not advance or reset any participant.
    controls.sort_by_key(|event| event.sample);
    Ok(controls)
}
