//! A material field actuator sliding against a fixed, thermally finite guide.
//!
//! The original field descriptor owns circuit, body/support inertia and elastic
//! storage. The original contact/enthalpy solve owns traction, physical friction
//! heat, both temperatures and slip history. Two displacement residuals identify
//! their SAME accepted motion and material query; no second inertia integrator,
//! prescribed drag coefficient, post-step heat credit or history advance occurs.
//! Electrical bulk/contact losses and friction heat enter the same thermal trial.
//!
//! Surface A is one complete source-bound field body and surface B is fixed.
//! Both tangents use the contact's existing first axis, positive toward increasing
//! field gap. Normal force/area are held and do no work. Frozen guide geometry,
//! body mass and optional conservative support keep their original restrictions.
//! Static reaction at rest and uninterrupted sliding are admitted. A step that
//! crosses zero velocity refuses: split at that event before continuing.

use super::*;
use crate::circuit_enthalpy::{CircuitEnthalpyError, ResistorCell};
use crate::circuit_enthalpy::contacts::{self, ContactHeatCells, ContactThermalCoefficients};
use crate::material_circuit::{MaterialCircuitSolve};
use crate::material_circuit::motion::{FieldMotionCandidate as NativeCandidate, FieldMotionInput,
    FieldMotionStep as NativeStep, MaterialFieldMotion};
use fs_tribo::FrictionRegime;
use std::collections::BTreeMap;

/// Original domain failures, including the caller's unchanged material error.
pub type Error<E> = SlidingError<CircuitEnthalpyError<E>>;

#[derive(Clone, Copy)]
pub struct Input<'a> {
    pub motion: &'a MaterialFieldMotion,
    pub electrical_ports: &'a [f64],
    /// Exterior tangent force, excluding the contact reaction [N].
    pub external_force_n: f64,
    pub normal_load_n: f64,
    pub resistor_cells: &'a [ResistorCell],
    pub contact_cells: &'a [ContactHeatCells],
}

/// Independent contact/thermal, descriptor and displacement/work gates.
#[derive(Clone, Copy, Debug)]
pub struct Controls {
    pub sliding: SlidingSolve,
    pub circuit: MaterialCircuitSolve,
    pub maximum_motion_trials: usize,
    /// All descriptor evaluations across ALL discarded outer/thermal trials.
    pub maximum_field_trials: usize,
    pub maximum_newton_iterations: usize,
    pub difference_step_m: f64,
    pub displacement_tolerance_m: f64,
    pub energy_tolerance_j: f64,
    pub newton: NewtonKrylovConfig,
}

pub struct Report<S> {
    /// Internal friction transfer, thermal stores and original source sample.
    /// Its external_heat_j includes electrical heat; external_heat_j below does not.
    pub sliding: SlidingReport<S>,
    pub field: NativeStep,
    pub nonlinear: NewtonReport,
    pub motion_trials: usize,
    pub field_trials: usize,
    /// Actual gap change minus native slip, and queried slip minus native slip [m].
    pub displacement_residual_m: [f64; 2],
    pub deposited_joule_j: f64,
    /// Independently evaluated authored/callback heat, excluding BOTH losses.
    pub external_heat_j: f64,
    pub external_mechanical_work_j: f64,
    pub contact_transfer_residual_j: f64,
    pub total_energy_residual_j: f64,
}

pub struct Candidate<'a, 't, S> {
    thermal: DepositionSimulation<'t>,
    contact: MaterialSlidingContact<'a>,
    motion: MaterialFieldMotion,
    report: Report<S>,
}
impl<'a, 't, S> Candidate<'a, 't, S> {
    pub const fn thermal(&self) -> &DepositionSimulation<'t> { &self.thermal }
    pub const fn contact(&self) -> &MaterialSlidingContact<'a> { &self.contact }
    pub const fn motion(&self) -> &MaterialFieldMotion { &self.motion }
    pub const fn report(&self) -> &Report<S> { &self.report }
    pub fn into_parts(self) -> (DepositionSimulation<'t>, MaterialSlidingContact<'a>,
        MaterialFieldMotion, Report<S>) {
        (self.thermal, self.contact, self.motion, self.report)
    }
}

struct Trial<S> {
    field: NativeCandidate,
    source: S,
    external_heat_j: f64,
    joule_j: f64,
}
fn finite<E>(x: f64) -> Result<f64, CircuitEnthalpyError<E>> {
    if x.is_finite() { Ok(x) }
    else { Err(CircuitEnthalpyError::Invalid("unrepresentable field/contact transfer")) }
}

impl<'a> MaterialSlidingContact<'a> {
    /// Close actual field displacement, dry traction, charge/flux and finite heat.
    /// The material callback is the original electrothermal source type and must
    /// supply only OTHER heat sources. Native Joule and friction heat are added
    /// here once, with explicit original branch/contact destinations.
    ///
    /// Guesses are [released displacement, material-query displacement] in m.
    /// Inertial compliance dt²/(2m) is an algebraic contact-solve coordinate, not
    /// an additional elastic interface or a second body response. Both native
    /// solved slip and actual field displacement must agree before publication.
    /// Positive/negative static force and sliding are supported; no unresolved
    /// stop/reversal, source fallback or fitted friction heat is accepted.
    #[allow(clippy::too_many_arguments)]
    pub fn field_motion_candidate<'t, S, E>(
        &self, thermal: &DepositionSimulation<'t>, target_s: f64, input: Input<'_>,
        controls: Controls, displacement_guess_m: [f64; 2], energy_guess_j: &[f64],
        resolve: impl Fn(&DepositionSimulation<'t>) -> Result<ContactThermalCoefficients<S>, E>,
        cancelled: impl FnMut() -> bool,
    ) -> Result<Candidate<'a, 't, S>, Error<E>> {
        let cancel = RefCell::new(cancelled);
        let checkpoint = || if (*cancel.borrow_mut())() { Err(DepositionError::Cancelled) } else { Ok(()) };
        checkpoint()?;
        controls.newton.validate()?;
        let dt = number(target_s-thermal.time_s())?;
        if self.time_s() != thermal.time_s() || dt <= 0.0
            || !input.external_force_n.is_finite()
            || displacement_guess_m.iter().any(|x| !x.is_finite())
            || controls.maximum_motion_trials == 0 || controls.maximum_field_trials == 0
            || controls.maximum_newton_iterations == 0
            || [controls.displacement_tolerance_m,controls.energy_tolerance_j]
                .iter().any(|x| !x.is_finite() || *x <= 0.0) {
            return Err(SlidingError::Invalid("invalid field/contact clocks, guesses or work controls"));
        }
        let regions = input.motion.body().regions();
        let cell = thermal.cells().get(&self.binding().cells[0])
            .ok_or(SlidingError::Invalid("moving contact thermal body is absent"))?;
        if regions.len() != 1
            || regions[0].material().card_identity() != self.binding().point.materials[0].content_hash()
            || regions[0].inertia().mass() != cell.mass_kg
            || regions[0].geometry().volume_m3() != cell.geometry.volume_m3 {
            return Err(SlidingError::Invalid("surface A must be the same complete single-region field body, source, mass and reference volume"));
        }
        let compliance = number(0.5*dt*(dt/input.motion.generalized_mass_kg()))?;
        if compliance <= 0.0 { return Err(SlidingError::Invalid("field contact inertial compliance is unrepresentable")); }
        let mut destinations = BTreeMap::new();
        if input.resistor_cells.len() != input.motion.circuit().resistors().len() {
            return Err(SlidingError::Invalid("every field resistor requires exactly one thermal destination"));
        }
        for route in input.resistor_cells {
            if !thermal.cells().contains_key(&route.cell)
                || destinations.insert(route.branch,route.cell).is_some()
                || !input.motion.circuit().resistors().iter().any(|r| r.branch == route.branch) {
                return Err(SlidingError::Invalid("duplicate, foreign or inactive field resistor destination"));
            }
        }
        let routes = contacts::admit_routes(thermal,input.motion.circuit(),input.contact_cells,checkpoint)
            .map_err(SlidingError::Source)?;
        let trials = Cell::new(0usize);
        let field_trials = Cell::new(0usize);
        let evaluate = |x: &[f64]| {
            checkpoint()?;
            if trials.get() >= controls.maximum_motion_trials { return Err(SlidingError::MotionBudget); }
            trials.set(trials.get()+1);
            let candidate = self.step_candidate_with_motion(thermal,target_s,
                SlidingDrive { normal_load_n: input.normal_load_n, relative_increment_m: [x[0],0.0] },
                compliance,[x[1],0.0],controls.sliding,energy_guess_j,|state,friction| {
                    if field_trials.get() >= controls.maximum_field_trials {
                        return Err(CircuitEnthalpyError::CircuitBudget);
                    }
                    let coefficients = resolve(state).map_err(CircuitEnthalpyError::Source)?;
                    contacts::admit_sources(input.motion.circuit(),&coefficients.contacts,
                        || (*cancel.borrow_mut())())?;
                    let bulk = coefficients.bulk;
                    let external_heat_j = finite(dt*state.external_step_power_w(
                        &bulk.sources_w,&bulk.conductances_w_k,checkpoint)?)?;
                    field_trials.set(field_trials.get()+1);
                    let field = input.motion.step_with_contacts_candidate(bulk.resistors,coefficients.contacts,
                        FieldMotionInput { electrical_ports: input.electrical_ports,
                            external_force_n: finite(input.external_force_n-friction.traction_n[0])?,
                            duration_s: dt },controls.circuit,|| (*cancel.borrow_mut())())?;
                    let mut powers = BTreeMap::new();
                    for heat in &field.step().resistors {
                        let cell = *destinations.get(&heat.material.branch)
                            .ok_or(CircuitEnthalpyError::Invalid("field resistor returned unrouted heat"))?;
                        let power = finite(heat.heat_j/dt)?;
                        if heat.heat_j != 0.0 && power == 0.0 {
                            return Err(CircuitEnthalpyError::Invalid("field Joule heat vanished in power conversion"));
                        }
                        let sum = powers.entry(cell).or_insert(0.0);
                        *sum = finite(*sum+power)?;
                    }
                    contacts::deposit(&field.step().contacts,&routes,dt,&mut powers,checkpoint)?;
                    let joule_j = finite(dt*powers.values().try_fold(0.0,|s,&p|finite(s+p))?)?;
                    if (joule_j-field.step().heat_j).abs() > controls.energy_tolerance_j {
                        return Err(CircuitEnthalpyError::Invalid("field Joule heat and receiving powers disagree"));
                    }
                    let mut sources_w = bulk.sources_w;
                    for (cell,power) in powers {
                        let sum = sources_w.entry(cell).or_insert(0.0);
                        *sum = finite(*sum+power)?;
                    }
                    Ok(StepCoefficients { sources_w,conductances_w_k:bulk.conductances_w_k,
                        payload: Trial { field,source:bulk.source,external_heat_j,joule_j } })
                },|| (*cancel.borrow_mut())()).map_err(|e| match e {
                    SlidingError::Source(CircuitEnthalpyError::Cancelled) => SlidingError::Cancelled,
                    other => other,
                })?;
            let report = candidate.report();
            let actual = number(report.source.field.model().gap_m()-input.motion.gap_m())?;
            let slip = report.friction.slip_increment_m[0];
            if report.friction.slip_increment_m[1] != 0.0 || report.friction.traction_n[1] != 0.0 {
                return Err(SlidingError::Invalid("field guide must use the original first contact tangent only"));
            }
            Ok((candidate,[number(actual-slip)?,number(x[1]-slip)?]))
        };
        let problem = FallibleFiniteDifference::new(2,controls.difference_step_m,
            |x: &[f64]| evaluate(x).map(|(_,r)|r.to_vec())).map_err(numerical)?;
        let construction = NewtonKrylovState::new(&problem,displacement_guess_m.to_vec(),controls.newton);
        if let Some(e) = problem.take_failure() { return Err(numerical(e)); }
        let mut solver = construction?;
        let result = solver.run_controlled(&problem,controls.maximum_newton_iterations,
            &mut || (*cancel.borrow_mut())() || problem.has_failed());
        if let Some(e) = problem.take_failure() { return Err(numerical(e)); }
        checkpoint()?;
        if result.progress != SolverRunProgress::Complete { return Err(SlidingError::Interrupted(result.progress)); }
        if !result.report.converged { return Err(SlidingError::NotConverged(result.report)); }
        let mut accepted = evaluate(&solver.x)?;
        // Sticking is the native exact zero-slip active set. A Newton tolerance
        // must not publish a fictitious speed/distance material source query.
        if accepted.0.report().friction.regime == FrictionRegime::Sticking && solver.x[1] != 0.0 {
            accepted = evaluate(&[solver.x[0],0.0])?;
            if accepted.0.report().friction.regime != FrictionRegime::Sticking {
                return Err(SlidingError::Invalid("zero-speed source does not admit the proposed native sticking state"));
            }
        }
        let (candidate,residual) = accepted;
        let material = &candidate.report().material;
        let queried_distance = number(dt*material.speed_mps)?;
        let solved_distance = candidate.report().friction.slip_increment_m[0].abs();
        let source_distance_residual = number(material.accumulated_distance_m
            -candidate.contact().accumulated_distance_m())?;
        if residual.iter().any(|r| r.abs() > controls.displacement_tolerance_m)
            || number(queried_distance-solved_distance)?.abs() > controls.displacement_tolerance_m
            || source_distance_residual.abs() > controls.displacement_tolerance_m {
            return Err(SlidingError::Invalid("field displacement, original solved slip and material query do not agree"));
        }
        let v0 = input.motion.velocity_m_s();
        let v1 = candidate.report().source.field.model().velocity_m_s();
        if v0 != 0.0 && ((v1 != 0.0 && v1.signum() != v0.signum())
            || candidate.report().friction.regime == FrictionRegime::Sticking) {
            return Err(SlidingError::Invalid("field contact step crosses an unresolved stop or reversal event"));
        }
        let (thermal,contact,sliding) = candidate.into_parts();
        let Trial { field,source,external_heat_j,joule_j } = sliding.source;
        let displacement = number(field.model().gap_m()-input.motion.gap_m())?;
        let mechanical_work = number(input.external_force_n*displacement)?;
        let transfer = number(sliding.drive_work_j-sliding.friction.traction_n[0]*displacement)?;
        let e = field.step();
        let total = number(e.electrical_energy_j[1]-e.electrical_energy_j[0]
            +e.kinetic_energy_j[1]-e.kinetic_energy_j[0]+e.elastic_energy_j[1]-e.elastic_energy_j[0]
            +sliding.thermal.stored_energy_change_j+sliding.friction.stored_change_j
            +sliding.partition.other_work_j()-e.electrical_supply_j-mechanical_work-external_heat_j)?;
        if transfer.abs() > controls.energy_tolerance_j || total.abs() > controls.energy_tolerance_j {
            return Err(SlidingError::Invalid("field, friction and finite thermal work do not close the same external supply"));
        }
        let sliding = SlidingReport { friction:sliding.friction,material:sliding.material,
            partition:sliding.partition,thermal:sliding.thermal,source,nonlinear:sliding.nonlinear,
            friction_trials:sliding.friction_trials,maximum_cell_residual_j:sliding.maximum_cell_residual_j,
            drive_work_j:sliding.drive_work_j,deposited_heat_j:sliding.deposited_heat_j,
            external_heat_j:sliding.external_heat_j,total_energy_residual_j:sliding.total_energy_residual_j };
        let (motion,field) = field.into_parts();
        checkpoint()?;
        Ok(Candidate { thermal,contact,motion,report:Report { sliding,field,nonlinear:result.report,
            motion_trials:trials.get(),field_trials:field_trials.get(),displacement_residual_m:residual,
            deposited_joule_j:joule_j,external_heat_j,external_mechanical_work_j:mechanical_work,
            contact_transfer_residual_j:transfer,total_energy_residual_j:total } })
    }
}
