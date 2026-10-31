//! Native Stokes and material response inside the original thermal residual.
//!
//! Unknowns are cell-energy increments [J]. Every residual evaluation closes
//! source rheology, native momentum, physical heating, directed enthalpy transport
//! and molecular conduction at the SAME proposed endpoint h/T. No trial advances
//! accepted history. BackwardEulerStep alone owns the storage update and its
//! local/global energy gates; the existing Newton--Krylov owner solves it.

use super::rheology::QuasisteadyRheologyReport;
use super::*;
use fs_conduction::enthalpy::backward_euler::{
    BackwardEulerStep, StepCoefficients, TransportStepCoefficients,
};
use fs_conduction::transient::deposition::EnergyAudit;
use fs_solver::finite_difference::{FallibleFiniteDifference, FiniteDifferenceError};
use fs_solver::{
    NewtonError, NewtonKrylovConfig, NewtonKrylovState, NewtonReport, SolverRunProgress,
};
use std::cell::{Cell, RefCell};

/// Numerical work and independent physical acceptance for one implicit interval.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ImplicitFlowControls {
    /// Includes thermal cells outside the moving-material map.
    pub max_unknowns: usize,
    pub max_newton_iterations: usize,
    /// Maximum-component perturbation in CELL energy, not specific enthalpy [J].
    pub difference_step_j: f64,
    pub cell_energy_tolerance_j: f64,
    pub total_energy_tolerance_j: f64,
    /// Residual coordinates are joules. Inner work remains explicitly bounded.
    pub newton: NewtonKrylovConfig,
}

#[derive(Debug)]
pub enum ImplicitFlowError {
    Flow(MaterialFlowError<DepositionError>),
    Invalid(&'static str),
    Newton(NewtonError),
    NotConverged(NewtonReport),
    Interrupted(SolverRunProgress),
}
impl From<MaterialFlowError<DepositionError>> for ImplicitFlowError {
    fn from(error: MaterialFlowError<DepositionError>) -> Self {
        Self::Flow(error)
    }
}
impl From<DepositionError> for ImplicitFlowError {
    fn from(error: DepositionError) -> Self {
        Self::Flow(error.into())
    }
}
impl From<NewtonError> for ImplicitFlowError {
    fn from(error: NewtonError) -> Self {
        Self::Newton(error)
    }
}
impl core::fmt::Display for ImplicitFlowError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Flow(e) => write!(f, "implicit material flow: {e}"),
            Self::Invalid(s) => write!(f, "implicit material-flow controls: {s}"),
            Self::Newton(e) => write!(f, "implicit material-flow Newton: {e}"),
            Self::NotConverged(r) => write!(
                f,
                "implicit material-flow residual {} J after {} attempts: {:?}",
                r.residual_norm, r.iterations, r.diagnosis
            ),
            Self::Interrupted(p) => write!(f, "implicit material-flow solve interrupted: {p:?}"),
        }
    }
}
impl std::error::Error for ImplicitFlowError {}

/// One accepted backward-Euler interval and its actual endpoint physical state.
/// These diagnostics do not become sources for the next interval.
#[derive(Debug)]
pub struct ImplicitFlowReport<F, T = FluidTransportState> {
    /// Includes open-stream enthalpy and physical heat exactly once.
    pub thermal: EnergyAudit,
    pub nonlinear: NewtonReport,
    pub maximum_cell_residual_j: f64,
    /// All native solves, including nonlinear/finite-difference trials and the
    /// final independent acceptance evaluation; bounded by the original owner.
    pub momentum_solves: usize,
    pub endpoint: MaterialFlowEpoch<F, T>,
    /// Observed rates are distinct from the queried coordinates in `endpoint`.
    pub solved_shear_rate_per_s: Option<Vec<f64>>,
    pub rheology: Option<QuasisteadyRheologyReport>,
    /// Same native wall/reservoir/body supply and physical heat, integrated at
    /// the accepted endpoint. Supply is not added again to thermal sources.
    pub mechanical: ViscousHeatingAudit,
    pub deposited_heat_j: f64,
    pub net_stream_enthalpy_j: f64,
    /// Independently evaluated authored heat and boundary supply, excluding
    /// stream enthalpy and mechanical-to-thermal transfer.
    pub other_thermal_input_j: f64,
    /// Delta h storage minus stream enthalpy, other heat and actuator work.
    pub combined_residual_j: f64,
}

struct Endpoint<F, T> {
    flow: MaterialFlowEpoch<F, T>,
    solved_rates: Option<Vec<f64>>,
    rheology: Option<QuasisteadyRheologyReport>,
    deposited_w: f64,
    stream_w: f64,
    other_w: f64,
}

fn numerical_error(
    error: FiniteDifferenceError<MaterialFlowError<DepositionError>>,
) -> ImplicitFlowError {
    match error {
        FiniteDifferenceError::Source(e) => e.into(),
        FiniteDifferenceError::Invalid(s) => ImplicitFlowError::Invalid(s),
    }
}

fn finite(value: f64) -> Result<f64, DepositionError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(DepositionError::Arithmetic)
    }
}

impl<F: Clone, R, T: FlowTransportSnapshot> MaterialFlowFeedback<'_, F, R, (), T>
where
    R: Clone + FnMut(usize, &DepositedCell, f64) -> Result<F, MaterialError>,
{
    /// Close ONE thermal interval against its own material response and flow.
    ///
    /// Bind the original viscous, body-force or reservoir-pressure heating mode
    /// first. Retained momentum, buoyancy and mutable wrapped observer histories
    /// are outside this stationary Stokes entry. Existing material, phase,
    /// porous and optional quasisteady-rheology sources remain authoritative.
    /// Callbacks must be repeatable and must not mutate shared physical history.
    /// Each trial clones this accepted material blueprint, including its source
    /// selections and rate seed. No nonlinear trial becomes a new starting state.
    ///
    /// Joule unknowns follow `thermal.cells()`'s stable CellId order. Their initial
    /// guess is not heat input. The SAME native advective graph supplies donor
    /// transfers and inlet/outlet enthalpy; they are never symmetric conductances.
    /// Existing full fluid-face conduction uses that endpoint's original source
    /// samples. Authored contacts, heat and reservoir boundaries retain meaning.
    /// No explicit CFL is imposed on this implicit residual. Algebraic convergence
    /// does not establish temporal accuracy, uniqueness or physical validity.
    ///
    /// All native work, including Jacobian probes, line-search trials and the
    /// final residual check, shares the original `max_solves` ceiling. Unsupported
    /// source stencils, budget exhaustion and cancellation refuse the whole
    /// interval. Topology/inventory are fixed; interior events refuse through the
    /// original thermal owner and target-time events remain pending.
    ///
    /// The return value contains an accepted thermal state and endpoint report,
    /// not a fabricated explicit observer history. Resume by passing that thermal
    /// state to this same blueprint. This remains reference-density creeping flow
    /// with no kinetic, moving-skeleton or free-surface store; refinement of time,
    /// space, constitutive closure and numerical Jacobian scale remains required.
    #[allow(clippy::too_many_arguments)]
    pub fn step_implicit<'t>(
        &self,
        thermal: &DepositionSimulation<'t>,
        target_s: f64,
        controls: ImplicitFlowControls,
        initial_energy_increments_j: &[f64],
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(DepositionSimulation<'t>, ImplicitFlowReport<F, T>), ImplicitFlowError> {
        let checkpoint = RefCell::new(checkpoint);
        let check = || check_gate(self.gate, &mut *checkpoint.borrow_mut());
        check()?;
        controls.newton.validate()?;
        if controls.max_unknowns == 0
            || thermal.cells().len() > controls.max_unknowns
            || initial_energy_increments_j.len() != thermal.cells().len()
            || controls.max_newton_iterations == 0
            || !controls.difference_step_j.is_finite()
            || controls.difference_step_j <= 0.0
            || !controls.cell_energy_tolerance_j.is_finite()
            || controls.cell_energy_tolerance_j < 0.0
            || !controls.total_energy_tolerance_j.is_finite()
            || controls.total_energy_tolerance_j < 0.0
        {
            return Err(ImplicitFlowError::Invalid(
                "invalid unknown/work budget, difference scale or joule tolerance",
            ));
        }
        if self.unsteady.is_some() || self.buoyancy.is_some() || self.viscous_heating.is_none() {
            return Err(ImplicitFlowError::Invalid(
                "implicit material flow needs the existing stationary Stokes heating owner",
            ));
        }
        self.check_map(thermal, &mut || check())?;
        let step = BackwardEulerStep::new(thermal, target_s, || check())?;
        let work = Cell::new(0usize);
        let evaluate = |x: &[f64]| {
            step.evaluate_with_transport(
                x,
                |state| {
                    let mut proposal = self.clone();
                    proposal.solves = work.get();
                    // Capture work even when a later source/rheology gate refuses.
                    let resolved = proposal.solve_rheology_epoch(state, &mut || check());
                    work.set(proposal.solves);
                    let (epoch, reference) = resolved?;
                    let active = proposal.momentum_domain(&epoch.samples, &mut || check())?;
                    let domain = active.as_ref().unwrap_or(proposal.domain);
                    let (flows, boundaries) = mass_graph(
                        domain,
                        &epoch.flow.field,
                        proposal.density_kg_m3,
                        proposal.input,
                        &mut || check(),
                    )?;
                    let transport = Advected::new((), flows, proposal.input.max_connections)?
                        .with_boundaries(boundaries, proposal.input.max_connections)?
                        .endpoint_transport(state, || check())?;
                    let mut sources_w = transport.heat.bulk.sources_w;
                    let stream_w = sources_w.values().try_fold(0.0, |s, &p| finite(s + p))?;
                    let power = epoch
                        .viscous_power
                        .as_ref()
                        .ok_or(invalid("implicit endpoint lacks physical Stokes power"))?;
                    let mut deposited_w = 0.0;
                    for sample in &epoch.samples {
                        check()?;
                        let watts = power.cell_heat_w[sample.voxel];
                        deposited_w = finite(deposited_w + watts)?;
                        let source = sources_w.entry(sample.cell).or_insert(0.0);
                        *source = finite(*source + watts)?;
                    }
                    proposal.reference = reference;
                    proposal.applied = Some(epoch);
                    let conductances_w_k = proposal.conductances(state, &mut || check())?;
                    let other_w =
                        state.external_step_power_w(&BTreeMap::new(), &conductances_w_k, || {
                            check()
                        })?;
                    let solved_rates = proposal.endpoint_shear_rate_per_s().map(<[f64]>::to_vec);
                    let rheology = proposal.quasisteady_rheology_report();
                    let endpoint = Endpoint {
                        flow: proposal.applied.take().expect("resolved endpoint retained"),
                        solved_rates,
                        rheology,
                        deposited_w,
                        stream_w,
                        other_w,
                    };
                    Ok::<_, MaterialFlowError<DepositionError>>(TransportStepCoefficients {
                        heat: StepCoefficients {
                            sources_w,
                            conductances_w_k,
                            payload: endpoint,
                        }
                        .into(),
                        transfers: transport.transfers,
                    })
                },
                || check(),
            )
        };
        let problem = FallibleFiniteDifference::new(
            step.cell_ids().len(),
            controls.difference_step_j,
            |x: &[f64]| evaluate(x).map(|trial| trial.residual_j().to_vec()),
        )
        .map_err(numerical_error)?;
        let construction = NewtonKrylovState::new(
            &problem,
            initial_energy_increments_j.to_vec(),
            controls.newton,
        );
        if let Some(error) = problem.take_failure() {
            return Err(numerical_error(error));
        }
        let mut solver = construction?;
        let result = solver.run_controlled(&problem, controls.max_newton_iterations, &mut || {
            self.gate.is_requested() || problem.has_failed()
        });
        if let Some(error) = problem.take_failure() {
            return Err(numerical_error(error));
        }
        check()?;
        if result.progress != SolverRunProgress::Complete {
            return Err(ImplicitFlowError::Interrupted(result.progress));
        }
        if !result.report.converged {
            return Err(ImplicitFlowError::NotConverged(result.report));
        }
        let trial = evaluate(&solver.x)?;
        let maximum_cell_residual_j = trial.maximum_cell_residual_j();
        let (next, endpoint, thermal_audit) = trial.accept(
            controls.cell_energy_tolerance_j,
            controls.total_energy_tolerance_j,
            || check(),
        )?;
        let dt = step.duration_s();
        let mut mechanical = self.viscous_heating.clone().expect("heating admitted");
        mechanical.begin();
        mechanical.accept(
            endpoint
                .flow
                .viscous_power
                .as_ref()
                .expect("endpoint power retained"),
            thermal.time_s(),
            dt,
            target_s,
        )?;
        let deposited_heat_j = finite(dt * endpoint.deposited_w)?;
        let net_stream_enthalpy_j = finite(dt * endpoint.stream_w)?;
        let other_thermal_input_j = finite(dt * endpoint.other_w)?;
        let work_j = finite(
            finite(
                mechanical.audit.supplied_wall_work_j + mechanical.audit.supplied_pressure_work_j,
            )? + mechanical.audit.supplied_body_work_j,
        )?;
        let combined_residual_j = finite(
            thermal_audit.stored_energy_change_j
                - finite(finite(net_stream_enthalpy_j + other_thermal_input_j)? + work_j)?,
        )?;
        let power_allowance = finite(dt * mechanical.controls.power_tolerance_w)?;
        let allowance = finite(controls.total_energy_tolerance_j + power_allowance)?;
        if (deposited_heat_j - mechanical.audit.generated_heat_j).abs() > power_allowance
            || combined_residual_j.abs() > allowance
        {
            return Err(invalid(
                "implicit physical heat, stream enthalpy and actuator work do not close",
            )
            .into());
        }
        check()?;
        Ok((
            next,
            ImplicitFlowReport {
                thermal: thermal_audit,
                nonlinear: result.report,
                maximum_cell_residual_j,
                momentum_solves: work.get(),
                endpoint: endpoint.flow,
                solved_shear_rate_per_s: endpoint.solved_rates,
                rheology: endpoint.rheology,
                mechanical: mechanical.audit,
                deposited_heat_j,
                net_stream_enthalpy_j,
                other_thermal_input_j,
                combined_residual_j,
            },
        ))
    }
}
