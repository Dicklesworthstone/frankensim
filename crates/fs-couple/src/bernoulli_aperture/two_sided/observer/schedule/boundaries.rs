//! Original boundary programs on the two-sided source's mechanical clock.
use super::{GestureCompileError, ScheduledTwoSidedPressure, TwoSidedPressure, invalid};
use crate::bernoulli_aperture::network::boundaries::{BoundaryDrive, BoundaryNetworkFrame};
use crate::bernoulli_aperture::performance::AperturePerformanceConfig;
use crate::bernoulli_aperture::performance::boundaries::{
    BoundaryControl, BoundaryProgram, BoundarySchedule,
};
use crate::bernoulli_aperture::two_sided::ValveSide;
use fs_scenario::assembly::{ParameterKey, ParameterPath};
use fs_scenario::entity::EntityDeclaration;
use fs_scenario::gesture::{GestureSchedule, GestureTarget};

/// Physical graph endpoint and its declared closing-force consumer.
#[derive(Clone, Debug)]
pub struct BoundaryForceBinding {
    /// Original node on the selected upstream or downstream graph.
    pub node: usize,
    /// Force in newtons; a positive value closes the aperture.
    pub closing_force: ParameterPath,
}

/// Explicit initial reservoirs/forces and typed force identities for one side.
#[derive(Clone, Debug)]
pub struct BoundaryGestureProgram {
    /// One complete initial drive per attached node in ascending order.
    pub initial: Vec<BoundaryDrive>,
    /// One force path per initial node; an absent track keeps the initial force.
    pub bindings: Vec<BoundaryForceBinding>,
}

/// Typed source pressure consumed by `control_source ASSEMBLY PART`.
#[must_use]
pub fn source_pressure_path(assembly: &str, part: &str) -> ParameterPath {
    let assembly = EntityDeclaration::assembly(assembly).identity();
    ParameterPath {
        entity: EntityDeclaration::part(assembly, part).identity(),
        parameter: ParameterKey::new("blowing_pressure", fs_scenario::gesture::path::PRESSURE),
    }
}

impl ScheduledTwoSidedPressure {
    /// Attach one complete program to an already keyed acoustic side. Events
    /// address mechanical samples, including substeps inside one audio frame.
    /// All keys retain their original material/contact and chimney histories.
    /// Zero force releases the actuator without clearing any physical state.
    ///
    /// Initial values and changes on BOTH sides share max_controls with pressure
    /// and source-plate forces. Boundary admission and source-force projection
    /// share max_compile_work, independently of the existing pressure compiler.
    /// Attachment order cannot renew either budget. A keyed side without an
    /// explicit program refuses playback before any physical sample advances.
    pub fn with_boundary_controls(
        mut self,
        side: ValveSide,
        program: BoundaryProgram,
    ) -> Result<Self, GestureCompileError> {
        if self.poisoned || self.source.poisoned
            || self.source.mechanical_samples_rendered() != 0
            || self.boundaries[side.index()].is_some() {
            return Err(invalid("boundary controls require one cold attachment per keyed acoustic side"));
        }
        let (assignments, work) = self.remaining_control_budgets()?;
        let boundaries = self.source.system.boundary_apertures(side)
            .ok_or_else(|| invalid("boundary controls require the side's original physical apertures"))?;
        let samples = self.config.samples.checked_mul(self.source.mechanics_per_output() as u64)
            .ok_or_else(|| invalid("boundary control mechanical window overflow"))?;
        self.boundaries[side.index()] = Some(BoundarySchedule::compile(
            boundaries, program, samples, assignments, work,
        )?);
        Ok(self)
    }

    pub(super) fn remaining_control_budgets(&self) -> Result<(usize, u64), GestureCompileError> {
        let mut assignments = self.config.max_controls.checked_sub(self.controls.len())
            .ok_or_else(|| invalid("pressure assignments exceed the combined storage budget"))?;
        let mut work = self.config.max_compile_work.checked_sub(self.gesture_work)
            .ok_or_else(|| invalid("typed gesture sampling exceeds auxiliary compilation work"))?;
        for force in self.plate_forces.iter().flatten() {
            assignments = assignments.checked_sub(force.events.len())
                .ok_or_else(|| invalid("source forces exceed the combined control storage budget"))?;
            work = work.checked_sub(force.work)
                .ok_or_else(|| invalid("source forces exceed the auxiliary compilation work budget"))?;
        }
        for pressure in self.spatial_pressures.iter().flatten() {
            assignments = assignments.checked_sub(pressure.events.len())
                .ok_or_else(|| invalid("spatial pressure exceeds the combined control storage budget"))?;
            work = work.checked_sub(pressure.work)
                .ok_or_else(|| invalid("spatial pressure exceeds auxiliary compilation work"))?;
        }
        for boundary in self.boundaries.iter().flatten() {
            assignments = assignments.checked_sub(boundary.assignment_count())
                .ok_or_else(|| invalid("boundaries exceed the combined control storage budget"))?;
            work = work.checked_sub(boundary.work)
                .ok_or_else(|| invalid("boundaries exceed the auxiliary compilation work budget"))?;
        }
        Ok((assignments, work))
    }

    /// Compile one all-typed pressure/key schedule with the existing fs-scenario
    /// interpolation owner. Positive key force closes the actual material slit.
    /// Every event executes before its addressed mechanical solve, including
    /// substeps inside one audio frame. No callback allocates control buffers.
    ///
    /// Declared keys without a force track retain their explicit initial drive.
    /// Unknown consumers, duplicate paths, negative closing force, late events
    /// and total sampling/storage exhaustion refuse before any source advances.
    pub fn with_typed_boundary_schedule(
        source: TwoSidedPressure,
        schedule: GestureSchedule,
        config: AperturePerformanceConfig,
        pressure_path: ParameterPath,
        programs: [Option<BoundaryGestureProgram>; 2],
    ) -> Result<Self, GestureCompileError> {
        let tracks = schedule.tracks();
        if !(1..=129).contains(&tracks.len())
            || pressure_path.parameter.name != "blowing_pressure"
            || pressure_path.parameter.dims != fs_scenario::gesture::path::PRESSURE
        {
            return Err(invalid(
                "typed reed controls need a declared pressure path and bounded source/key tracks",
            ));
        }
        let mechanics = source.mechanics_rate_hz();
        let samples = config
            .samples
            .checked_mul(source.mechanics_per_output() as u64)
            .ok_or_else(|| invalid("typed reed mechanical window overflow"))?;
        let rate = schedule.control_rate_hz;
        if rate == 0 || rate > mechanics || samples == 0 {
            return Err(invalid(
                "typed controls need a finite window and a clock no faster than mechanics",
            ));
        }
        let last_tick = u128::from(samples - 1) * u128::from(rate) / u128::from(mechanics);
        let ticks = last_tick + 1;
        if ticks * tracks.len() as u128 > config.max_controls as u128
            || tracks.iter().any(|t| {
                t.events
                    .iter()
                    .any(|event| event.time_s > last_tick as f64 / f64::from(rate))
            })
        {
            return Err(invalid(
                "typed source/key assignments exceed their finite window or storage allowance",
            ));
        }
        let per_tick: u128 = tracks
            .iter()
            .map(|t| tracks.len() as u128 + t.events.len() as u128 + 1)
            .sum();
        let required = ticks
            .checked_mul(per_tick)
            .ok_or_else(|| invalid("typed source/key compilation work overflow"))?;
        if required > u128::from(config.max_compile_work) {
            return Err(GestureCompileError::WorkBudget {
                required,
                allowed: config.max_compile_work,
            });
        }
        let mut paths = Vec::new();
        for program in programs.iter().flatten() {
            if program.initial.len() != program.bindings.len() || program.initial.len() > 64 {
                return Err(invalid(
                    "each physical key needs one initial drive and one typed force binding",
                ));
            }
            for (drive, binding) in program.initial.iter().zip(&program.bindings) {
                if drive.node != binding.node
                    || drive.opening_force_n > 0.0
                    || binding.closing_force.parameter.name != "closing_force"
                    || binding.closing_force.parameter.dims != fs_scenario::gesture::path::FORCE
                    || paths.contains(&&binding.closing_force)
                {
                    return Err(invalid(
                        "key force bindings must be unique, dimensioned and in complete node order",
                    ));
                }
                paths.push(&binding.closing_force);
            }
        }
        let mut pressure = None;
        let mut force_tracks = Vec::new();
        for track in tracks {
            let GestureTarget::Parameter(path) = &track.target else {
                return Err(invalid(
                    "typed source/key schedules cannot mix legacy aliases",
                ));
            };
            if path == &pressure_path {
                if pressure.replace(track).is_some() {
                    return Err(invalid(
                        "one physical pressure consumer cannot have two tracks",
                    ));
                }
                continue;
            }
            let address = programs.iter().enumerate().find_map(|(side, program)| {
                program.as_ref().and_then(|p| {
                    p.bindings
                        .iter()
                        .position(|b| &b.closing_force == path)
                        .map(|key| (side, key))
                })
            });
            let (side, key) = address
                .ok_or_else(|| invalid("typed control has no declared physical key consumer"))?;
            if force_tracks.iter().any(|(_, s, k)| *s == side && *k == key) {
                return Err(invalid(
                    "one physical key cannot have two closing-force tracks",
                ));
            }
            force_tracks.push((track, side, key));
        }
        let mut pressure = pressure
            .ok_or_else(|| invalid("typed source/key schedule requires pressure"))?
            .clone();
        pressure.target = GestureTarget::BlowingPressure;
        let pressure_schedule =
            GestureSchedule::try_new(rate, vec![pressure]).map_err(GestureCompileError::Gesture)?;
        let mut compiled: [Option<BoundaryProgram>; 2] = programs.each_ref().map(|p| {
            p.as_ref().map(|p| BoundaryProgram {
                initial: p.initial.clone(),
                controls: Vec::new(),
            })
        });
        let tick_count =
            u64::try_from(ticks).map_err(|_| invalid("typed control tick count exceeds u64"))?;
        for (track, side, key) in force_tracks {
            let program = compiled[side].as_mut().expect("bound physical side");
            let initial = program.initial[key];
            let mut previous = None;
            for tick in 0..tick_count {
                let force = schedule
                    .sample(&track.id, tick)
                    .map_err(GestureCompileError::Gesture)?;
                if !force.is_finite() || force < 0.0 {
                    return Err(invalid(
                        "physical key closing force must be finite and nonnegative",
                    ));
                }
                if previous != Some(force.to_bits()) {
                    let sample =
                        (u128::from(tick) * u128::from(mechanics)).div_ceil(u128::from(rate));
                    program
                        .controls
                        .try_reserve(1)
                        .map_err(|_| invalid("cannot reserve bounded key controls"))?;
                    program.controls.push(BoundaryControl {
                        sample: u64::try_from(sample)
                            .map_err(|_| invalid("key control sample exceeds u64"))?,
                        drive: BoundaryDrive {
                            opening_force_n: -force,
                            ..initial
                        },
                    });
                    previous = Some(force.to_bits());
                }
            }
        }
        let mut renderer = Self::new(source, pressure_schedule, config)?;
        renderer.schedule = schedule;
        renderer.gesture_work = required as u64;
        for (side, program) in [ValveSide::Upstream, ValveSide::Downstream]
            .into_iter()
            .zip(compiled)
        {
            if let Some(program) = program {
                renderer = renderer.with_boundary_controls(side, program)?;
            }
        }
        Ok(renderer)
    }

    /// Last accepted held drives, or explicit initial drives before playback.
    #[must_use]
    pub fn held_boundary_drives(&self, side: ValveSide) -> &[BoundaryDrive] {
        self.boundaries[side.index()].as_ref().map_or(&[], BoundarySchedule::held)
    }
    /// Events published with successfully observed mechanical samples only.
    #[must_use]
    pub fn applied_boundary_controls(&self, side: ValveSide) -> &[BoundaryControl] {
        self.boundaries[side.index()].as_ref().map_or(&[], BoundarySchedule::applied)
    }
    /// Remaining events in mechanical (sample, node) order.
    #[must_use]
    pub fn pending_boundary_controls(&self, side: ValveSide) -> &[BoundaryControl] {
        self.boundaries[side.index()].as_ref().map_or(&[], BoundarySchedule::pending)
    }
    /// Last key/graph work account ingested by the original pressure receiver.
    #[must_use]
    pub fn last_boundary_frame(&self, side: ValveSide) -> Option<&BoundaryNetworkFrame> {
        self.source.last_frame()?.boundaries[side.index()].as_ref()
    }
}
