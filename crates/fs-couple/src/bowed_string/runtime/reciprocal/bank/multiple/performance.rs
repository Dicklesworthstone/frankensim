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
        let performance = self.into_scheduled_admitted(source.schedule().clone(), tracks, config,Some(paths),false)?;
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
        self.into_scheduled_admitted(schedule,tracks,config,None,false)
    }

    /// Drive a finite ribbon from ONE original centre-station gesture track.
    ///
    /// Its speed acts at every cell; any prescribed total normal force is divided
    /// by the explicit spatial quadrature. If cells use solved normal contact,
    /// the original track must supply zero normal force; each cell then receives
    /// its own simultaneous normal reaction without quadrature rescaling.
    /// Every cell keeps an independent friction
    /// regime in the joint solve. The original schedule is retained verbatim,
    /// including release, reversal and interrupted ramps. No per-cell or audio
    /// normalization is applied beyond distributing the declared total load.
    /// Use `into_scheduled` for independently prescribed cell drives instead.
    pub fn into_scheduled_ribbon(
        self, schedule: GestureSchedule, track: String, config: MultiBowPerformanceConfig,
    ) -> Result<MultiBowPerformance, GestureCompileError> {
        let ribbon = self.source.ribbon().ok_or_else(|| source_invalid(
            "ribbon scheduling requires a source-bound finite-width contact"))?;
        if schedule.tracks().len() != 1 || ribbon.stations().len() != self.source.bows().len() {
            return Err(source_invalid("one original gesture track must drive the complete ribbon"));
        }
        let tracks = vec![track; ribbon.stations().len()];
        self.into_scheduled_admitted(schedule,tracks,config,None,true)
    }

    fn into_scheduled_admitted(self,schedule:GestureSchedule,tracks:Vec<String>,config:MultiBowPerformanceConfig,
        paths:Option<&[fs_scenario::assembly::ParameterPath]>,ribbon_projection:bool)->Result<MultiBowPerformance,GestureCompileError> {
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
            if !ribbon_projection && tracks[..port].contains(name) {
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
            let (station_fraction,load_fraction) = if ribbon_projection {
                let ribbon = self.source.ribbon().ok_or_else(|| source_invalid("missing ribbon source"))?;
                (ribbon.centre().station_fraction,ribbon.load_fraction())
            } else { (bow.location.station_fraction,1.0) };
            let compiled = compile_friction_track_admitted(&schedule, name, FrictionGestureConfig {
                sample_rate_hz: rate, samples: config.samples, max_block: config.max_block,
                max_work: remaining_work, max_events: config.max_events-events.len(),
                station_fraction,
            },path)?;
            if self.system.contact_normal_load(port).is_some()
                && compiled.events.iter().any(|e| e.drive.normal_force_n != 0.0) {
                return Err(source_invalid("contact-loaded bow tracks must not prescribe an additional normal force"));
            }
            if self.system.material_shear_layer(port).is_some()
                && self.system.contact_normal_load(port).is_none()
                && compiled.events.iter().any(|e| e.drive.normal_force_n <= 0.0) {
                return Err(source_invalid("retained material shear layers require positive normal load throughout playback"));
            }
            remaining_work -= compiled.required_work;
            events.try_reserve(compiled.events.len())
                .map_err(|_| source_invalid("cannot reserve complete bow controls"))?;
            events.extend(compiled.events.into_iter().map(|e| BowDriveEvent { sample: e.sample, port,
                drive: FrictionDrive { speed_m_s:e.drive.speed_m_s,
                    normal_force_n:e.drive.normal_force_n*load_fraction } }));
        }
        events.sort_unstable_by_key(|e| (e.sample, e.port));
        let count = self.source.bows().len();
        let forces = vec![0.0; self.system.mode_count()];
        Ok(MultiBowPerformance { bank: self, schedule, tracks, config, events, next: 0,
            base_sample, forces, drives: vec![FrictionDrive { speed_m_s: 0.0, normal_force_n: 0.0 };count],
            candidate: vec![FrictionDrive { speed_m_s: 0.0, normal_force_n: 0.0 };count], poisoned: false,
            string_forces: None, stop_motion: None, stop_station_motion: None, bow_station_motion: None,
            fiber_support: None, shell_observer: None })
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
    stop_station_motion: Option<stations::BoundStopStations>,
    bow_station_motion: Option<bow_stations::BoundBowStations>,
    fiber_support: Option<fiber_motion::BoundFiberSupport>,
    shell_observer: Option<observation::BoundShellObserver>,
}
impl MultiBowPerformance {
    /// Add fiber support acceleration and material contact travel to the SAME
    /// performance clock. Bind after optional bow-station paths. Schedule speeds
    /// denote support endpoint velocities; acceleration limits apply to every
    /// compiled control transition. A control tick activates a target at its
    /// unchanged sample boundary; constant acceleration reaches that target at
    /// the following mechanical endpoint (one physical interval later).
    pub fn with_fiber_support_motion(mut self, initial_velocity_m_s: f64,
        limits: super::ribbon::motion::FiberSupportMotionLimits,
        gate: &fs_exec::CancelGate) -> Result<Self, RenderError> {
        self.validate_sample_count(0)?;
        self.bank = self.bank.with_fiber_support_motion(initial_velocity_m_s, limits, gate)?;
        Ok(self)
    }
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
        self.block_with_step(output, |bank, forces, drives, motion| {
            let frame = if let Some(bows) = motion.bow_stations {
                bank.step_with_bow_and_stop_stations_under_gate(forces, drives, bows,
                    motion.stop_stations, motion.stop_offsets_m, &fs_exec::CancelGate::new_clock_free())
            } else if let Some(stations) = motion.stop_stations {
                bank.step_with_stop_stations(forces, drives, stations, motion.stop_offsets_m)
            } else if let Some(offsets) = motion.stop_offsets_m {
                bank.step_with_stop_offsets(forces, drives, offsets).map_err(RenderError::Friction)
            } else {
                bank.step(forces, drives).map_err(RenderError::Friction)
            }?;
            Ok(frame.observer_pressure_pa)
        })
    }
}

// Prepared physical endpoint controls, borrowed from the one canonical sampler.
// Every step owner must receive them; a specialized thermal step cannot be
// bypassed when a finger presses or migrates during the same performance.
struct StepMotion<'a> {
    bow_stations: Option<&'a [f64]>,
    stop_stations: Option<&'a [f64]>,
    stop_offsets_m: Option<&'a [f64]>,
}

impl MultiBowPerformance {
    // One schedule, force, motion and observer loop for ordinary and thermally
    // coupled rendering. A specialized step must commit the SAME bank once.
    fn block_with_step(&mut self, output: &mut [f64],
        mut step: impl FnMut(&mut MultiBowStringBank, &[f64], &[FrictionDrive], StepMotion<'_>)
            -> Result<f64, RenderError>) -> Result<(), RenderError> {
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
            let drives = if let Some(support) = self.fiber_support.as_mut() {
                if let Err(error) = support.prepare(sample+1, &self.candidate) {
                    self.poisoned = true;
                    return Err(error);
                }
                support.drives.as_slice()
            } else { self.candidate.as_slice() };
            let next_motion = self.stop_motion.as_mut().map(|motion| motion.prepare(sample+1));
            let next_stations = self.stop_station_motion.as_mut().map(|motion| motion.prepare(sample+1));
            let next_bow_stations = self.bow_station_motion.as_mut().map(|motion| motion.prepare(sample+1));
            if let (Some(bows), Some(stops)) = (&self.bow_station_motion, &mut self.stop_station_motion) {
                // A mapped normal and its tangent have one physical location.
                // The independent stop program retains source placeholders at
                // these indices; the bow program owns the complete patch path.
                for binding in &self.bank.source.bow_contact_loads {
                    stops.candidate[binding.stop] = bows.candidate[binding.bow];
                }
            }
            let motion = StepMotion {
                bow_stations: self.bow_station_motion.as_ref().map(|motion| motion.candidate.as_slice()),
                stop_stations: self.stop_station_motion.as_ref().map(|motion| motion.candidate.as_slice()),
                stop_offsets_m: self.stop_motion.as_ref().map(|motion| motion.candidate.as_slice()),
            };
            let result = step(&mut self.bank, forces, drives, motion);
            match result {
                Ok(compact_pressure) => {
                    if let Some(support) = self.fiber_support.as_mut() {
                        support.commit();
                    }
                    *slot = if let Some(observer) = self.shell_observer.as_mut() {
                        match observer.pressure(&self.bank.system) {
                            Ok(pressure) => pressure,
                            Err(error) => { self.poisoned = true; return Err(error); }
                        }
                    } else { compact_pressure };
                    self.drives.copy_from_slice(&self.candidate); self.next = next;
                    if let (Some(loads), Some(cursor)) = (self.string_forces.as_mut(), next_force) {
                        self.forces.copy_from_slice(&loads.candidate); loads.next = cursor;
                    }
                    if let (Some(motion), Some(cursor)) = (self.stop_motion.as_mut(), next_motion) { motion.next = cursor; }
                    if let (Some(motion), Some(cursor)) = (self.stop_station_motion.as_mut(), next_stations) { motion.next = cursor; }
                    if let (Some(motion), Some(cursor)) = (self.bow_station_motion.as_mut(), next_bow_stations) { motion.next = cursor; }
                },
                Err(error) => { self.poisoned = true; return Err(error); }
            }
        }
        Ok(())
    }
}

/// Physical load-and-release excitation of the same material instrument.
pub mod force;

/// Normal-support trajectories on the same physical sample clock.
pub mod motion;

/// Longitudinal stop trajectories using source-derived maps on that same clock.
pub mod stations;

/// Prescribed bow-position trajectories on the same physical sample clock.
pub mod bow_stations;

/// Material-fiber support acceleration and laboratory energy conversion.
pub mod fiber_motion;

/// Source-bound finite-distance closed-shell pressure on the same accepted clock.
pub mod observation;


/// Source-bound contact heat on the original scheduled pressure/PCM route.
pub mod thermal;
