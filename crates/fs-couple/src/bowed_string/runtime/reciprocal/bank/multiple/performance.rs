//! Scheduled simultaneous bow contacts on the original source-bound string bank.
//! Geometry, materials, contact admission and joint mechanics remain in the
//! parent owner. All due source-track controls act before one mechanical step.

use super::MultiBowStringBank;
use crate::pcm_wav::observation::PressureRenderer;
use crate::render::{RenderError, schedule::GestureCompileError};
use crate::render::schedule::force::coupled::contact::friction::FrictionDrive;
use crate::render::schedule::force::coupled::render::friction::{
    FrictionGestureConfig, compile_friction_track_admitted,
};
use fs_scenario::gesture::GestureSchedule;

fn source_invalid(what: &'static str) -> GestureCompileError { GestureCompileError::Invalid { what } }
fn render_invalid(what: &'static str) -> RenderError { RenderError::Control { what } }

impl MultiBowStringBank {
    /// Bind one typed path per actual contact, in the source's contact order.
    /// Return the v2 source for provenance; all physical steps remain unchanged.
    pub fn into_scheduled_paths(
        self,
        source: &fs_scenario::gesture::ResolvedSchedule,
        paths: &[fs_scenario::assembly::ParameterPath],
        config: MultiBowPerformanceConfig,
    ) -> Result<(MultiBowPerformance, crate::gesture_path::BoundPaths), GestureCompileError> {
        use crate::gesture_path::{PathBinding, PathInput, bind_paths};
        if paths.len() != self.source.bows().len() {
            return Err(source_invalid(
                "one typed path is required per physical bow contact",
            ));
        }
        let bindings: Vec<_> = paths
            .iter()
            .zip(self.source.bows())
            .map(|(path, bow)| PathBinding {
                path: path.clone(),
                input: PathInput::Bow {
                    member: bow.location.string,
                },
            })
            .collect();
        let bound = bind_paths(source, &bindings).map_err(GestureCompileError::Gesture)?;
        let tracks = paths
            .iter()
            .map(|p| source.track(p).map(|t| t.id.clone()))
            .collect::<Result<Vec<_>, _>>()
            .map_err(GestureCompileError::Gesture)?;
        let performance = self.into_scheduled_admitted(source.schedule().clone(), tracks, config,Some(paths))?;
        Ok((performance, bound))
    }
    /// Compile a complete finite performance on the unchanged physical clock.
    ///
    /// Exactly one distinct track per contact, in `source().bows()` order, must name its
    /// actual member and fixed station. Every track uses the common friction
    /// compiler, including interrupted ramps and nondivisor control rates.
    /// Events due on the same sample all act before ONE joint friction step.
    /// Source input and numerical limits remain inspectable. This starts only
    /// a local performance timeline; preexisting vibration is never reset.
    ///
    /// No force/gain normalization is applied. Existing body/cavity pressure
    /// approximations are inherited. Include release and ringdown inside the
    /// finite horizon; rendering beyond it refuses rather than silently stopping
    /// friction or extending a held drive. Cancellation belongs to the existing
    /// PressureRenderer stream writer at callback boundaries.
    ///
    /// # Errors
    /// Incomplete/wrong/duplicate tracks, station/clock mismatch, source inputs,
    /// total compilation/event budgets, unrepresentable horizon or allocation.
    pub fn into_scheduled(
        self,
        schedule: GestureSchedule,
        tracks: Vec<String>,
        config: MultiBowPerformanceConfig,
    ) -> Result<MultiBowPerformance, GestureCompileError> {
        self.into_scheduled_admitted(schedule,tracks,config,None)
    }

    fn into_scheduled_admitted(self,schedule:GestureSchedule,tracks:Vec<String>,config:MultiBowPerformanceConfig,
        paths:Option<&[fs_scenario::assembly::ParameterPath]>)->Result<MultiBowPerformance,GestureCompileError> {
        let rate = self.source.bank().options().sample_rate_hz;
        let base_sample = self.system.samples_rendered();
        if config.samples == 0 || config.max_block == 0 || tracks.len() != self.source.bows().len()
            || self.system.sample_period_s().to_bits() != f64::from(rate).recip().to_bits()
            || base_sample.checked_add(config.samples).is_none() {
            return Err(source_invalid("multi-bow performance needs complete contacts and a bounded horizon on its physical clock"));
        }
        let mut events = Vec::new();
        let mut remaining_work = config.max_work;
        for (port, (name, bow)) in tracks.iter().zip(self.source.bows()).enumerate() {
            if tracks[..port].contains(name) {
                return Err(source_invalid("each physical bow requires its own explicit track binding"));
            }
            let source_track = schedule.tracks().iter().find(|t| &t.id == name);
            let path=paths.and_then(|p|p.get(port));
            let accepts=match path {
                Some(path)=>source_track.is_some_and(|t|crate::gesture_path::matches_input(t,crate::gesture_path::PathInput::Bow {member:bow.location.string},Some(path))),
                None=>source_track.and_then(|t|fs_scenario::gesture::path::legacy_member(&t.target,"bow"))==Some(bow.location.string),
            };
            if !accepts {
                return Err(source_invalid("bow track must address the actual material string at that contact"));
            }
            let compiled = compile_friction_track_admitted(&schedule, name, FrictionGestureConfig {
                sample_rate_hz: rate, samples: config.samples, max_block: config.max_block,
                max_work: remaining_work, max_events: config.max_events-events.len(),
                station_fraction: bow.location.station_fraction,
            },path)?;
            remaining_work -= compiled.required_work;
            events.try_reserve(compiled.events.len())
                .map_err(|_| source_invalid("cannot reserve complete bow controls"))?;
            events.extend(compiled.events.into_iter().map(|e| BowDriveEvent { sample: e.sample, port, drive: e.drive }));
        }
        events.sort_unstable_by_key(|e| (e.sample, e.port));
        let count = self.source.bows().len();
        let forces = vec![0.0; self.system.mode_count()];
        Ok(MultiBowPerformance { bank: self, schedule, tracks, config, events, next: 0,
            base_sample, forces, drives: vec![FrictionDrive { speed_m_s: 0.0, normal_force_n: 0.0 };count],
            candidate: vec![FrictionDrive { speed_m_s: 0.0, normal_force_n: 0.0 };count], poisoned: false,
            string_forces: None, stop_motion: None })
    }
}

/// Finite performance/storage admission. The sample rate comes from the bank.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct MultiBowPerformanceConfig {
    /// Complete local output window, including any desired ringdown.
    pub samples: u64,
    /// Largest host callback, preserving the existing output cancellation bound.
    pub max_block: usize,
    /// Shared source-track/event visit allowance across ALL contacts.
    pub max_work: u64,
    /// Shared maximum emitted controls across ALL contacts (also bounds sorting).
    pub max_events: usize,
}
/// One complete drive assignment for a physical contact, not an independent voice.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BowDriveEvent {
    /// Sample relative to performance start, not the original absolute clock.
    pub sample: u64,
    /// Index in the preserved contact list.
    pub port: usize,
    /// Held physical speed and normal force.
    pub drive: FrictionDrive,
}

/// Shared PressureRenderer output of one source-retaining multi-contact system.
/// A failed callback poisons the wrapper; the mechanics retain their last
/// complete joint step, not a falsely resumable partially written audio block.
pub struct MultiBowPerformance {
    bank: MultiBowStringBank,
    schedule: GestureSchedule,
    tracks: Vec<String>,
    config: MultiBowPerformanceConfig,
    events: Vec<BowDriveEvent>,
    next: usize,
    base_sample: u64,
    forces: Vec<f64>,
    drives: Vec<FrictionDrive>,
    candidate: Vec<FrictionDrive>,
    poisoned: bool,
    string_forces: Option<force::BoundStringForces>,
    stop_motion: Option<motion::BoundStopMotion>,
}
impl MultiBowPerformance {
    /// Original material bank, actual contact maps and accepted joint physics.
    #[must_use] pub const fn bank(&self) -> &MultiBowStringBank { &self.bank }
    /// Original authoring schedule, with no rewritten ramps or time labels.
    #[must_use] pub const fn schedule(&self) -> &GestureSchedule { &self.schedule }
    /// Complete explicit track-to-contact binding, in physical contact order.
    #[must_use] pub fn tracks(&self) -> &[String] { &self.tracks }
    /// Actual finite window and aggregate numerical budgets.
    #[must_use] pub const fn config(&self) -> MultiBowPerformanceConfig { self.config }
    /// Compiled physical drives; inspecting them cannot advance the system.
    #[must_use] pub fn events(&self) -> &[BowDriveEvent] { &self.events }
    /// Remaining finite output samples on this performance's local clock.
    #[must_use] pub fn remaining_samples(&self) -> u64 { self.config.samples-self.samples_rendered() }
}
impl PressureRenderer for MultiBowPerformance {
    fn samples_rendered(&self) -> u64 { self.bank.system.samples_rendered()-self.base_sample }
    fn max_block_len(&self) -> usize { self.config.max_block }
    fn validate_sample_rate(&self, rate: u32) -> Result<(), RenderError> {
        self.validate_sample_count(0)?;
        if rate != self.bank.source.bank().options().sample_rate_hz {
            return Err(render_invalid("multi-bow output must use the original physical sample rate"));
        }
        Ok(())
    }
    fn validate_sample_count(&self, samples: u64) -> Result<(), RenderError> {
        if self.poisoned { return Err(RenderError::Poisoned); }
        if samples > self.remaining_samples() { return Err(render_invalid("multi-bow request exceeds its finite performance window")); }
        Ok(())
    }
    fn block(&mut self, output: &mut [f64]) -> Result<(), RenderError> {
        let count = u64::try_from(output.len()).map_err(|_| render_invalid("multi-bow output size exceeds its clock"))?;
        self.validate_sample_count(count)?;
        if output.is_empty() { return Err(RenderError::EmptyBlock); }
        if output.len() > self.config.max_block { return Err(render_invalid("multi-bow callback exceeds its admitted capacity")); }
        for slot in output {
            let sample = self.samples_rendered();
            let mut next = self.next;
            self.candidate.copy_from_slice(&self.drives);
            while let Some(event) = self.events.get(next) {
                if event.sample != sample { break; }
                self.candidate[event.port] = event.drive; next += 1;
            }
            let mut next_force = None;
            let forces = if let Some(loads) = self.string_forces.as_mut() {
                if loads.projected.events().get(loads.next).is_some_and(|e| e.sample == sample) {
                    let mut cursor = loads.next;
                    loads.candidate.copy_from_slice(&self.forces);
                    while let Some(event) = loads.projected.events().get(cursor) {
                        if event.sample != sample { break; }
                        loads.candidate[event.mode] = event.force_n_per_sqrt_kg;
                        cursor += 1;
                    }
                    next_force = Some(cursor);
                    &loads.candidate
                } else { &self.forces }
            } else { &self.forces };
            let mut next_motion = None;
            let result = if let Some(motion) = self.stop_motion.as_mut() {
                next_motion = Some(motion.prepare(sample + 1));
                self.bank.system.step_with_obstacle_offsets(forces, &self.candidate, &motion.candidate)
            } else { self.bank.system.step(forces, &self.candidate) };
            match result {
                Ok(frame) => {
                    *slot = frame.observer_pressure_pa;
                    self.drives.copy_from_slice(&self.candidate); self.next = next;
                    if let (Some(loads), Some(cursor)) = (self.string_forces.as_mut(), next_force) {
                        self.forces.copy_from_slice(&loads.candidate); loads.next = cursor;
                    }
                    if let (Some(motion), Some(cursor)) = (self.stop_motion.as_mut(), next_motion) {
                        motion.next = cursor;
                    }
                },
                Err(error) => { self.poisoned = true; return Err(RenderError::Friction(error)); }
            }
        }
        Ok(())
    }
}

/// Physical load-and-release excitation of the same material instrument.
pub mod force;

/// Prescribed normal support trajectories on the same physical sample clock.
pub mod motion;
