//! Coupled source-backed contact Newton on the original enthalpy balance.
//! No thermal evolution or material derivative is reimplemented here: each
//! physical query uses StepContext and the original contact-tangent owner.

use std::cell::Cell;

use fs_solver::Globalization;
use fs_sparse::Coo;

use super::*;
use crate::interface::temperature::enthalpy::ContactEnthalpyNewtonConfig;
use crate::interface::temperature::{TemperatureContactEvaluation, TemperatureContacts};

struct Physical {
    h: Vec<f64>,
    contacts: TemperatureContactEvaluation,
}

struct Tangent {
    h: Vec<f64>,
    matrix: Csr,
    inverse: Vec<f64>,
    slope: Vec<f64>,
    radiation: Option<PointwiseBoundaryAssembly>,
}

struct CoupledProblem<'a, 'm, 'c, 'cx> {
    base: StepContext<'a, 'm, 'c, 'cx>,
    contacts: &'a TemperatureContacts<'a, 'a>,
    config: ContactEnthalpyNewtonConfig,
    physical: RefCell<Option<Physical>>,
    tangent: RefCell<Option<Tangent>>,
    fatal: RefCell<Option<EnthalpyError>>,
    trial_refusal: RefCell<Option<EnthalpyError>>,
    evaluations: Cell<usize>,
    assemblies: Cell<usize>,
}

fn same_field(a: &[f64], b: &[f64]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(a, b)| a.to_bits() == b.to_bits())
}

impl<'m, 'c, 'cx> CoupledProblem<'_, 'm, 'c, 'cx> {
    fn context<'s>(&'s self, interfaces: &'s ThermalInterfaces) -> StepContext<'s, 'm, 'c, 'cx> {
        StepContext {
            interfaces: Some(interfaces),
            storage: self.base.storage,
            cx: self.base.cx,
            problem: self.base.problem,
            old: self.base.old,
            dt: self.base.dt,
            pointwise_radiation: self.base.pointwise_radiation,
        }
    }

    fn fail(&self, error: EnthalpyError) {
        let mut failure = self.fatal.borrow_mut();
        if failure.is_none() {
            *failure = Some(error);
        }
    }

    fn stopped(&self) -> bool {
        if self.fatal.borrow().is_some() {
            return true;
        }
        if let Err(error) = poll(self.base.cx, self.evaluations.get()) {
            self.fail(error);
            return true;
        }
        false
    }

    fn evaluate(&self, h: &[f64], output: &mut [f64]) -> Result<Physical, (EnthalpyError, bool)> {
        let recoverable = |error| {
            let domain = matches!(
                error,
                EnthalpyError::Phase {
                    source: PhaseStateError::OutsideEnthalpyDomain { .. },
                    ..
                } | EnthalpyError::Conduction(ConductionError::OutsideTemperatureSpan { .. })
            ) || matches!(
                &error,
                EnthalpyError::Conduction(ConductionError::MaterialQuery { property, .. })
                    if self.base.pointwise_radiation.is_some()
                        && property == crate::SURFACE_EMISSIVITY_PROPERTY
            );
            (error, domain)
        };
        let temperature = self.base.temperatures(h).map_err(recoverable)?;
        let mut source_domain = false;
        let contacts = self
            .contacts
            .evaluate_with_query_refusal(
                self.base.cx,
                self.base.problem.boundary,
                &temperature,
                &mut |error| {
                    source_domain = matches!(
                        error,
                        fs_matdb::MatDbError::NoClaimInDomain { .. }
                            | fs_matdb::MatDbError::PinnedClaimOutOfDomain { .. }
                            | fs_matdb::MatDbError::OutsideKnotSpan { .. }
                            | fs_matdb::MatDbError::OutsideSurfaceSupport { .. }
                    );
                },
            )
            .map_err(|error| (error.into(), source_domain))?;
        let context = self.context(contacts.frozen_interfaces());
        let scales = context.phase_scales(h).map_err(recoverable)?;
        let system = context
            .assemble(&temperature, scales.as_deref())
            .map_err(recoverable)?;
        context
            .assembled_residual(h, &temperature, &system, output)
            .map_err(recoverable)?;
        Ok(Physical {
            h: h.to_vec(),
            contacts,
        })
    }

    fn prepare_tangent(&self, h: &[f64]) -> Result<(), EnthalpyError> {
        if self
            .tangent
            .borrow()
            .as_ref()
            .is_some_and(|t| same_field(&t.h, h))
        {
            return Ok(());
        }
        let count = self
            .assemblies
            .get()
            .checked_add(1)
            .filter(|&n| n <= self.config.step.max_newton_iterations)
            .ok_or(EnthalpyError::InvalidInput(
                "contact enthalpy tangent allowance exhausted",
            ))?;
        self.assemblies.set(count);
        let physical = self.physical.borrow();
        let physical = physical.as_ref().filter(|p| same_field(&p.h, h)).ok_or(
            EnthalpyError::InvalidInput("contact tangent has no matching physical state"),
        )?;
        let context = self.context(physical.contacts.frozen_interfaces());
        let stage = context.stage(h)?;
        let temperature = context.temperatures(h)?;
        let contact = self.contacts.linearize_temperature_for_newton(
            self.base.cx,
            self.base.problem.boundary,
            &temperature,
            self.config.max_face_pairs,
        )?;
        let n = h.len();
        let mut full = Coo::new(n, n);
        for row in 0..n {
            if row % ASSEMBLY_TILE == 0 {
                poll(self.base.cx, row)?;
            }
            full.push(row, row, self.base.storage.masses[row]);
            if context.is_prescribed(row) {
                continue;
            }
            let (columns, values) = stage.tangent.row(row);
            for (&column, &value) in columns.iter().zip(values) {
                full.push(row, column, finite(self.base.dt * value)?);
            }
            // stage already includes K_contact*dT/dh. Add ONLY the source
            // feedback dK_contact/dT*T, then convert its column K to J/kg.
            // Multiplication by dt gives the original joule residual tangent.
            let (columns, values) = contact.feedback().row(row);
            for (&column, &value) in columns.iter().zip(values) {
                full.push(
                    row,
                    column,
                    finite(self.base.dt * finite(value * stage.slope[column])?)?,
                );
            }
        }
        *self.tangent.borrow_mut() = Some(Tangent {
            h: h.to_vec(),
            matrix: full.assemble(),
            inverse: stage.inverse,
            slope: stage.slope,
            radiation: stage.radiation,
        });
        Ok(())
    }
}

impl NonlinearProblem for CoupledProblem<'_, '_, '_, '_> {
    fn dimension(&self) -> usize {
        self.base.old.len()
    }

    fn residual(&self, h: &[f64], output: &mut [f64]) {
        output.fill(f64::NAN);
        if self.stopped() {
            return;
        }
        let Some(count) = self
            .evaluations
            .get()
            .checked_add(1)
            .filter(|&n| n <= self.config.max_residual_evaluations)
        else {
            self.fail(EnthalpyError::InvalidInput(
                "contact enthalpy residual allowance exhausted",
            ));
            return;
        };
        self.evaluations.set(count);
        match self.evaluate(h, output) {
            Ok(physical) => {
                *self.physical.borrow_mut() = Some(physical);
                self.trial_refusal.borrow_mut().take();
            }
            Err((error, domain)) => {
                if domain {
                    *self.trial_refusal.borrow_mut() = Some(error);
                } else {
                    self.fail(error);
                }
                output.fill(f64::NAN);
            }
        }
    }

    fn jacobian_apply(&self, h: &[f64], direction: &[f64], output: &mut [f64]) {
        output.fill(f64::NAN);
        if self.stopped() {
            return;
        }
        if let Err(error) = self.prepare_tangent(h) {
            self.fail(error);
            return;
        }
        let tangent = self.tangent.borrow();
        let tangent = tangent.as_ref().expect("prepared contact tangent");
        let matrix = &tangent.matrix;
        for (row, value) in output.iter_mut().enumerate() {
            if row % ASSEMBLY_TILE == 0 && self.stopped() {
                return;
            }
            let (columns, entries) = matrix.row(row);
            *value = columns
                .iter()
                .zip(entries)
                .map(|(&c, &entry)| entry * direction[c])
                .sum();
        }
        if let Some(radiation) = &tangent.radiation
            && let Err(error) = radiation.add_scaled_tangent(
                self.base.cx,
                direction,
                Some(&tangent.slope),
                self.base.dt,
                |row| !self.base.is_prescribed(row),
                output,
            )
        {
            self.fail(error.into());
            output.fill(f64::NAN);
        }
    }

    fn preconditioner_apply(
        &self,
        h: &[f64],
        _: usize,
        _: usize,
        residual: &[f64],
        output: &mut [f64],
    ) {
        output.fill(f64::NAN);
        if self.stopped() {
            return;
        }
        if let Err(error) = self.prepare_tangent(h) {
            self.fail(error);
            return;
        }
        let tangent = self.tangent.borrow();
        let inverse = &tangent.as_ref().expect("prepared contact tangent").inverse;
        for (row, value) in output.iter_mut().enumerate() {
            if row % ASSEMBLY_TILE == 0 && self.stopped() {
                return;
            }
            *value = residual[row] * inverse[row];
        }
    }
}

impl EnthalpyBackwardEuler<'_, '_> {
    #[allow(clippy::too_many_arguments)]
    pub(super) fn advance_with_temperature_contacts(
        &self,
        cx: &Cx<'_>,
        problem: ConductionProblem<'_>,
        contacts: &TemperatureContacts<'_, '_>,
        old_h: &[f64],
        dt_s: f64,
        config: ContactEnthalpyNewtonConfig,
        pointwise_radiation: Option<PointwiseBoundary<'_>>,
    ) -> Result<(EnthalpyStepSolution, TemperatureContactEvaluation), EnthalpyError> {
        self.admit_step(cx, problem, old_h, dt_s, config.step)?;
        contacts.check_mesh(problem.mesh)?;
        let radiation_matches = match pointwise_radiation {
            Some(PointwiseBoundary::Ambient(owner)) => {
                owner.matches(problem.mesh, problem.boundary)
            }
            Some(PointwiseBoundary::MaterialAmbient(owner)) => {
                owner.matches(problem.mesh, problem.boundary)
                    && owner.matches_curves(cx, |vertex| self.curve_for_vertex(vertex))?
            }
            Some(PointwiseBoundary::Enclosure(owner)) => {
                owner.matches(problem.mesh, problem.boundary)
            }
            None => true,
        };
        if !radiation_matches {
            return Err(EnthalpyError::InvalidInput(
                "pointwise radiation must retain the exact transport mesh and boundary",
            ));
        }
        if config.max_residual_evaluations == 0
            || config.max_face_pairs == 0
            || !config.residual_tolerance_j.is_finite()
            || config.residual_tolerance_j <= 0.0
            || config.step.newton.relative_tolerance != 0.0
            || config.step.newton.absolute_tolerance > config.residual_tolerance_j
            || !matches!(
                config.step.newton.globalization,
                Globalization::LineSearch(_)
            )
        {
            return Err(EnthalpyError::InvalidInput(
                "contact enthalpy Newton requires positive work/joule gates, line search, zero relative tolerance and a numerical target no looser than the physical target",
            ));
        }
        let base = StepContext {
            storage: self,
            cx,
            problem,
            interfaces: None,
            old: old_h,
            dt: dt_s,
            pointwise_radiation,
        };
        let initial = base.initial_values()?;
        let coupled = CoupledProblem {
            base,
            contacts,
            config,
            physical: RefCell::new(None),
            tangent: RefCell::new(None),
            fatal: RefCell::new(None),
            trial_refusal: RefCell::new(None),
            evaluations: Cell::new(0),
            assemblies: Cell::new(0),
        };
        let initialized = NewtonKrylovState::new(&coupled, initial, config.step.newton);
        if let Some(error) = coupled.fatal.borrow_mut().take() {
            return Err(error);
        }
        let mut newton = initialized.map_err(|error| {
            coupled
                .trial_refusal
                .borrow_mut()
                .take()
                .unwrap_or(EnthalpyError::Newton(error))
        })?;
        let run = newton.run_controlled(&coupled, config.step.max_newton_iterations, &mut || {
            coupled.stopped()
        });
        if let Some(error) = coupled.fatal.borrow_mut().take() {
            return Err(error);
        }
        match run.progress {
            SolverRunProgress::Complete => {}
            SolverRunProgress::Paused => {
                return Err(ConductionError::Cancelled {
                    stage: "contact-enthalpy-newton",
                    at: run.report.iterations,
                }
                .into());
            }
            SolverRunProgress::CallbackPanicked(callback) => {
                return Err(EnthalpyError::SolverCallbackPanicked(callback));
            }
            SolverRunProgress::DimensionMismatch { expected, actual } => {
                return Err(EnthalpyError::Newton(NewtonError::Dimension {
                    problem: actual,
                    state: expected,
                }));
            }
        }
        if !run.report.converged {
            return Err(coupled
                .trial_refusal
                .borrow_mut()
                .take()
                .unwrap_or(EnthalpyError::NotConverged(run.report)));
        }
        let physical = coupled
            .physical
            .borrow_mut()
            .take()
            .filter(|p| same_field(&p.h, &newton.x))
            .ok_or(EnthalpyError::InvalidInput(
                "accepted contact enthalpy lost its physical source state",
            ))?;
        let accepted = coupled
            .context(physical.contacts.frozen_interfaces())
            .finish(newton.x, run.report, config.step.energy_tolerance_j)?;
        poll(cx, accepted.newton.iterations)?;
        Ok((accepted, physical.contacts))
    }
}
