//! Source-owned directional friction and physical heat on finite caloric bodies.
//!
//! The original associated kinetic ellipse supplies traction, actual slip and
//! heat. Every thermal trial refreshes its SAME four selected coefficient
//! claims at the original two body temperatures and declared pressure,
//! sliding-speed and path coordinates. The caller explicitly supplies the
//! effective material direction and heat allocation; neither is inferred from
//! material names. [`dynamics`] closes the same contact with finite inertia.

pub mod dynamics;

use super::*;
use fs_contact::interface_binding::thermal_friction::elliptic::{
    EllipticFrictionMode, SourceThermalEllipticCoulomb, SourceThermalEllipticState,
};
use fs_contact::interface_binding::thermal_friction::{
    FrictionQueryAxisKind, SourceThermalCoulomb,
};
use fs_matdb::QueryPoint;
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState};
use fs_material::thermal::interface::ThermalInterfacePoint;
use fs_qty::Temperature;
use fs_thermal::enthalpy::CellId;
use fs_tribo::{HeatPartition, WorkLedger};
use std::sync::Arc;

/// Original ordered thermal owners and an explicit fixed nominal contact patch.
#[derive(Clone, Copy, Debug)]
pub struct DirectionalSlidingInterface<'a> {
    pub point: ThermalInterfacePoint<'a>,
    pub caloric: [&'a EquilibriumEnthalpyPhaseCurve; 2],
    pub cells: [CellId; 2],
    pub nominal_area_m2: f64,
    /// Fixed coordinates for body motion and the effective material direction.
    pub tangent_frame_id: &'a str,
    /// Explicit fractions of actual dissipative work, including an outside share.
    pub partition: HeatPartition,
    /// Required when the partition exports work outside the finite bodies.
    pub outside_destination: Option<&'a str>,
}

/// Pinned source owner, accepted native directional history and physical work.
#[derive(Clone)]
pub struct DirectionalSlidingContact<'a> {
    data: DirectionalSlidingInterface<'a>,
    law: Arc<SourceThermalEllipticCoulomb>,
    state: SourceThermalEllipticState,
    work: WorkLedger,
    distance_m: f64,
    time_s: f64,
}

/// Original native response and finite-thermal accounts for one candidate.
pub struct DirectionalSlidingReport<S> {
    pub friction: TangentialStep<SourceThermalEllipticState>,
    pub phases: [EquilibriumPhaseState; 2],
    pub pressure_pa: f64,
    pub partition: DissipationStep,
    pub thermal: EnergyAudit,
    pub source: S,
    pub nonlinear: NewtonReport,
    pub friction_trials: usize,
    pub maximum_cell_residual_j: f64,
    /// Source-query path increment, independently matched by a motion consumer.
    pub source_distance_m: f64,
    pub drive_work_j: f64,
    pub deposited_heat_j: f64,
    pub external_heat_j: f64,
    pub total_energy_residual_j: f64,
}

/// One atomic thermal/contact proposal; callers cannot publish a half-step.
pub struct DirectionalSlidingCandidate<'a, 't, S> {
    thermal: DepositionSimulation<'t>,
    contact: DirectionalSlidingContact<'a>,
    report: DirectionalSlidingReport<S>,
}
impl<'a, 't, S> DirectionalSlidingCandidate<'a, 't, S> {
    pub const fn thermal(&self) -> &DepositionSimulation<'t> {
        &self.thermal
    }
    pub const fn contact(&self) -> &DirectionalSlidingContact<'a> {
        &self.contact
    }
    pub const fn report(&self) -> &DirectionalSlidingReport<S> {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        DepositionSimulation<'t>,
        DirectionalSlidingContact<'a>,
        DirectionalSlidingReport<S>,
    ) {
        (self.thermal, self.contact, self.report)
    }
}

impl<'a> DirectionalSlidingContact<'a> {
    /// Bind actual solid caloric states to the already selected original claims.
    /// Source-axis declarations identify their physical owners. An unowned
    /// evolving coordinate refuses; a fixed coordinate never silently changes.
    pub fn new(
        thermal: &DepositionSimulation<'_>,
        data: DirectionalSlidingInterface<'a>,
        law: SourceThermalEllipticCoulomb,
        initial_normal_load_n: f64,
        initial_speed_mps: f64,
        tangential: TangentialSolveConfig,
    ) -> Result<Self, SlidingError> {
        if data.cells[0] == data.cells[1]
            || !data.nominal_area_m2.is_finite()
            || data.nominal_area_m2 <= 0.0
            || data.tangent_frame_id.trim().is_empty()
            || data.partition.other() > 0.0
                && data.outside_destination.is_none_or(|s| s.trim().is_empty())
        {
            return Err(SlidingError::Invalid(
                "directional contact needs distinct cells, a finite patch/frame and named outside work",
            ));
        }
        for source in law.principal_sources() {
            if source.interface().card_identity() != data.point.interface.content_hash()
                || !data
                    .point
                    .is_thermal_axis(&source.axes().temperature, Temperature::DIMS)
            {
                return Err(SlidingError::Invalid(
                    "directional sources must retain the original interface and an owned contact temperature",
                ));
            }
            for axis in source.query_axes() {
                if axis.kind == FrictionQueryAxisKind::Temperature
                    && !data.point.is_thermal_axis(&axis.name, Temperature::DIMS)
                {
                    return Err(SlidingError::Invalid(
                        "directional evolving temperature must belong to an original caloric surface or declared mean",
                    ));
                }
            }
        }
        let mut contact = Self {
            data,
            state: law.initial_state(),
            law: Arc::new(law),
            work: WorkLedger::default(),
            distance_m: 0.0,
            time_s: thermal.time_s(),
        };
        let (phases, pressure, query) =
            contact.sample(thermal, initial_normal_load_n, initial_speed_mps, 0.0)?;
        let _ = phases;
        // Zero work establishes the actual initial receipts without advancing
        // either physical clock or replacing the initially pinned claims.
        let input = contact.input(
            1.0,
            initial_normal_load_n,
            [0.0; 2],
            0.0,
            pressure,
            &query,
            tangential,
        )?;
        contact.state = contact
            .law
            .advance_at_query(
                &contact.state,
                &input,
                &query,
                EllipticFrictionMode::StaticOrSliding,
            )?
            .state;
        Ok(contact)
    }

    pub const fn binding(&self) -> &DirectionalSlidingInterface<'a> {
        &self.data
    }
    pub fn source_law(&self) -> &SourceThermalEllipticCoulomb {
        &self.law
    }
    pub const fn state(&self) -> &SourceThermalEllipticState {
        &self.state
    }
    pub const fn work(&self) -> &WorkLedger {
        &self.work
    }
    pub const fn accumulated_distance_m(&self) -> f64 {
        self.distance_m
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }

    fn motion_dependent(&self) -> bool {
        self.law.principal_sources()[0]
            .query_axes()
            .iter()
            .any(|axis| {
                matches!(
                    axis.kind,
                    FrictionQueryAxisKind::SlidingSpeed
                        | FrictionQueryAxisKind::AccumulatedDistance
                )
            })
    }

    fn sample<E>(
        &self,
        thermal: &DepositionSimulation<'_>,
        normal: f64,
        speed: f64,
        distance: f64,
    ) -> Result<([EquilibriumPhaseState; 2], f64, QueryPoint), SlidingError<E>> {
        if !normal.is_finite()
            || normal < 0.0
            || !speed.is_finite()
            || speed < 0.0
            || !distance.is_finite()
            || distance < 0.0
        {
            return Err(SlidingError::Invalid(
                "directional load, speed and path must be finite and nonnegative",
            ));
        }
        let phase = |side: usize| {
            let cell = thermal
                .cells()
                .get(&self.data.cells[side])
                .ok_or(SlidingError::Invalid("inactive directional thermal body"))?;
            let curve = self.data.caloric[side];
            if !thermal.uses_material_law(cell.material_id, curve) {
                return Err(SlidingError::Invalid(
                    "directional contact must retain original caloric objects",
                ));
            }
            let state = curve.state_at_specific_enthalpy(cell.specific_h_j_kg)?;
            if state.liquid_mass_fraction() != 0.0 {
                return Err(SlidingError::Invalid(
                    "dry directional friction does not admit partially melted or liquid surfaces",
                ));
            }
            Ok::<_, SlidingError<E>>(state)
        };
        let phases = [phase(0)?, phase(1)?];
        let pressure = number(normal / self.data.nominal_area_m2)?;
        if normal > 0.0 && pressure == 0.0 {
            return Err(SlidingError::Invalid(
                "directional nominal pressure underflow",
            ));
        }
        let mut query = self.data.point.at(phases)?;
        let source: &SourceThermalCoulomb = &self.law.principal_sources()[0];
        if let Some(axis) = &source.axes().nominal_pressure {
            query = query
                .with(axis, pressure)
                .map_err(|e| ThermalMaterialError::Coordinate(e.to_string()))?;
        }
        for axis in source.query_axes() {
            let value = match axis.kind {
                FrictionQueryAxisKind::Temperature => continue,
                FrictionQueryAxisKind::NominalPressure => pressure,
                FrictionQueryAxisKind::SlidingSpeed => speed,
                FrictionQueryAxisKind::AccumulatedDistance => distance,
            };
            query = query
                .with(&axis.name, value)
                .map_err(|e| ThermalMaterialError::Coordinate(e.to_string()))?;
        }
        Ok((phases, pressure, query))
    }

    #[allow(clippy::too_many_arguments)]
    fn input<E>(
        &self,
        dt: f64,
        normal: f64,
        increment: [f64; 2],
        compliance: f64,
        pressure: f64,
        query: &QueryPoint,
        solve: TangentialSolveConfig,
    ) -> Result<TangentialStepInput, SlidingError<E>> {
        let temperature = query
            .axes()
            .get(&self.law.principal_sources()[0].axes().temperature)
            .copied()
            .ok_or(SlidingError::Invalid(
                "missing original directional contact-temperature coordinate",
            ))?;
        Ok(TangentialStepInput {
            normal_n: normal,
            slip_increment_m: increment,
            dt_s: dt,
            compliance_m_per_n: compliance,
            contact_temperature_k: Some(temperature),
            nominal_pressure_pa: self.law.principal_sources()[0]
                .axes()
                .nominal_pressure
                .as_ref()
                .map(|_| pressure),
            solve,
        })
    }

    /// Prescribed physical slip through the original source and enthalpy owners.
    /// Cut at thermal source/topology events first. Other heat and conduction
    /// are supplied by the callback; it must not include friction a second time.
    #[allow(clippy::too_many_arguments)]
    pub fn step_candidate<'t, S, E>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        drive: SlidingDrive,
        controls: SlidingSolve,
        energy_guess_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'t>) -> Result<StepCoefficients<S>, E>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<DirectionalSlidingCandidate<'a, 't, S>, SlidingError<E>> {
        let distance = number(det::hypot(
            drive.relative_increment_m[0],
            drive.relative_increment_m[1],
        ))?;
        self.step_candidate_with_motion(
            thermal,
            target_s,
            drive,
            0.0,
            distance,
            EllipticFrictionMode::StaticOrSliding,
            controls,
            energy_guess_j,
            resolve,
            cancelled,
        )
    }

    /// Mechanical compliance is an original work-conjugate body input. A
    /// compliant consumer must match source_distance_m to the solved actual
    /// slip before accepting this private contact/thermal proposal.
    #[allow(clippy::too_many_arguments, clippy::too_many_lines)]
    fn step_candidate_with_motion<'t, S, E>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        drive: SlidingDrive,
        compliance: f64,
        source_distance_m: f64,
        mode: EllipticFrictionMode,
        controls: SlidingSolve,
        energy_guess_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'t>) -> Result<StepCoefficients<S>, E>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<DirectionalSlidingCandidate<'a, 't, S>, SlidingError<E>> {
        let cancel = RefCell::new(cancelled);
        let checkpoint = || {
            if (*cancel.borrow_mut())() {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        };
        checkpoint()?;
        controls.newton.validate()?;
        if self.time_s != thermal.time_s()
            || controls.maximum_cells == 0
            || thermal.cells().len() > controls.maximum_cells
            || energy_guess_j.len() != thermal.cells().len()
            || energy_guess_j.iter().any(|v| !v.is_finite())
            || controls.maximum_friction_trials == 0
            || controls.maximum_newton_iterations == 0
            || !compliance.is_finite()
            || compliance < 0.0
            || !source_distance_m.is_finite()
            || source_distance_m < 0.0
            || [
                controls.cell_energy_tolerance_j,
                controls.thermal_energy_tolerance_j,
                controls.total_energy_tolerance_j,
            ]
            .iter()
            .any(|v| !v.is_finite() || *v < 0.0)
        {
            return Err(SlidingError::Invalid(
                "invalid directional clocks, work controls or thermal unknowns",
            ));
        }
        let step = BackwardEulerStep::new(thermal, target_s, checkpoint)?;
        let dt = step.duration_s();
        let speed = number(source_distance_m / dt)?;
        let distance = number(self.distance_m + source_distance_m)?;
        if source_distance_m > 0.0 && (speed == 0.0 || distance == self.distance_m) {
            return Err(SlidingError::Invalid(
                "directional source motion vanished in speed or accumulated path",
            ));
        }
        let trials = Cell::new(0usize);
        let evaluate = |x: &[f64]| {
            step.evaluate(
                x,
                |state| {
                    checkpoint()?;
                    if trials.get() >= controls.maximum_friction_trials {
                        return Err(SlidingError::FrictionBudget);
                    }
                    trials.set(trials.get() + 1);
                    let (phases, pressure, query) =
                        self.sample(state, drive.normal_load_n, speed, distance)?;
                    let input = self.input(
                        dt,
                        drive.normal_load_n,
                        drive.relative_increment_m,
                        compliance,
                        pressure,
                        &query,
                        controls.tangential,
                    )?;
                    let friction = self
                        .law
                        .advance_at_query(&self.state, &input, &query, mode)?;
                    let work = number(
                        number(friction.traction_n[0] * friction.slip_increment_m[0])?
                            + number(friction.traction_n[1] * friction.slip_increment_m[1])?,
                    )?;
                    let compatibility = det::hypot(
                        friction.compatibility_residual_m[0],
                        friction.compatibility_residual_m[1],
                    );
                    if !compatibility.is_finite()
                        || number(friction.stored_change_j + friction.dissipated_j - work)?.abs()
                            > controls.total_energy_tolerance_j
                    {
                        return Err(SlidingError::Invalid(
                            "native directional physical work does not close",
                        ));
                    }
                    let partition = friction.dissipation_step(self.data.partition)?;
                    let StepCoefficients {
                        mut sources_w,
                        conductances_w_k,
                        payload,
                    } = resolve(state).map_err(SlidingError::Source)?;
                    let external_heat = number(
                        dt * state.external_step_power_w(
                            &sources_w,
                            &conductances_w_k,
                            checkpoint,
                        )?,
                    )?;
                    let mut deposited = 0.0;
                    for (id, heat) in self
                        .data
                        .cells
                        .into_iter()
                        .zip([partition.surface_a_heat_j(), partition.surface_b_heat_j()])
                    {
                        let power = number(heat / dt)?;
                        if heat != 0.0 && power == 0.0 {
                            return Err(SlidingError::Invalid(
                                "directional heat vanished during power conversion",
                            ));
                        }
                        deposited = number(deposited + number(dt * power)?)?;
                        let source = sources_w.entry(id).or_insert(0.0);
                        *source = number(*source + power)?;
                    }
                    if number(deposited + partition.other_work_j() - partition.total_work_j())?
                        .abs()
                        > controls.total_energy_tolerance_j
                    {
                        return Err(SlidingError::Invalid(
                            "directional heat destinations do not close native dissipation",
                        ));
                    }
                    Ok::<_, SlidingError<E>>(StepCoefficients {
                        sources_w,
                        conductances_w_k,
                        payload: (
                            phases,
                            pressure,
                            friction,
                            partition,
                            payload,
                            work,
                            deposited,
                            external_heat,
                        ),
                    })
                },
                checkpoint,
            )
        };
        let problem = FallibleFiniteDifference::new(
            step.cell_ids().len(),
            controls.difference_step_j,
            |x: &[f64]| evaluate(x).map(|t| t.residual_j().to_vec()),
        )
        .map_err(numerical)?;
        let construction =
            NewtonKrylovState::new(&problem, energy_guess_j.to_vec(), controls.newton);
        if let Some(error) = problem.take_failure() {
            return Err(numerical(error));
        }
        let mut solver = construction?;
        let result =
            solver.run_controlled(&problem, controls.maximum_newton_iterations, &mut || {
                (*cancel.borrow_mut())() || problem.has_failed()
            });
        if let Some(error) = problem.take_failure() {
            return Err(numerical(error));
        }
        checkpoint()?;
        if result.progress != SolverRunProgress::Complete {
            return Err(SlidingError::Interrupted(result.progress));
        }
        if !result.report.converged {
            return Err(SlidingError::NotConverged(result.report));
        }
        let trial = evaluate(&solver.x)?;
        let maximum_cell_residual_j = trial.maximum_cell_residual_j();
        let (
            next,
            (phases, pressure, friction, partition, source, work, deposited, external_heat),
            thermal_audit,
        ) = trial.accept(
            controls.cell_energy_tolerance_j,
            controls.thermal_energy_tolerance_j,
            checkpoint,
        )?;
        let total = number(
            thermal_audit.stored_energy_change_j
                + friction.stored_change_j
                + partition.other_work_j()
                - work
                - external_heat,
        )?;
        if total.abs() > controls.total_energy_tolerance_j {
            return Err(SlidingError::Invalid(
                "directional finite heat, outside work and physical slip do not balance",
            ));
        }
        let actual_distance = number(det::hypot(
            friction.slip_increment_m[0],
            friction.slip_increment_m[1],
        ))?;
        let mut contact = self.clone();
        contact.distance_m = number(self.distance_m + actual_distance)?;
        if actual_distance > 0.0 && contact.distance_m == self.distance_m {
            return Err(SlidingError::Invalid(
                "actual directional path increment is unrepresentable",
            ));
        }
        contact.work.record(partition)?;
        contact.state = friction.state.clone();
        contact.time_s = target_s;
        checkpoint()?;
        Ok(DirectionalSlidingCandidate {
            thermal: next,
            contact,
            report: DirectionalSlidingReport {
                friction,
                phases,
                pressure_pa: pressure,
                partition,
                thermal: thermal_audit,
                source,
                nonlinear: result.report,
                friction_trials: trials.get(),
                maximum_cell_residual_j,
                source_distance_m,
                drive_work_j: work,
                deposited_heat_j: deposited,
                external_heat_j: external_heat,
                total_energy_residual_j: total,
            },
        })
    }
}
