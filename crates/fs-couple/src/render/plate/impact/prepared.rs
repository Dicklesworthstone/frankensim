//! Prepared execution of the existing nonlinear impact model.
//!
//! Preparation changes scratch ownership, not the constitutive model or basis.
//! Both execution paths finish through the same physical acceptance gate below.
mod substeps;
pub use substeps::{ImpactSubstepConfig, ImpactSubstepReport, SubsteppedImpactSystem};
use std::ops::Deref;
use fs_exec::CancelGate;
use fs_material::fiber::{Uniaxial, WoolFeltState};
use fs_phs::{PreparedStepError, PreparedStepRecord, StepWorkspace};
use super::{ImpactError, ImpactFrame, ImpactSystem, invalid};
use super::audio::ImpactSource;

/// Persistent nonlinear mechanics with reusable solver and material-trial buffers.
///
/// Construction allocates; stepping does not allocate in this host or its solver.
/// The supplied storage callbacks must independently avoid allocation. This is
/// finite-difference Newton or analytic Newton with exact acoustic-pair
/// condensation. This is NOT a measured hard-real-time contract.
/// The model, basis, timestep, contact, Kelvin memory and felt conditioning are
/// unchanged. Only accepted steps update physical history and the sample counter.
pub struct PreparedImpactSystem {
    inner: ImpactSystem,
    workspace: StepWorkspace,
    candidate: Vec<f64>,
    output: Vec<f64>,
    histories: Vec<WoolFeltState>,
    analytic_newton: bool,
}

impl ImpactSystem {
    /// Prepare this exact system, preserving any already accepted motion/history.
    ///
    /// # Errors
    /// Returns the numerical owner's workspace preparation refusal.
    pub fn prepare(self) -> Result<PreparedImpactSystem, ImpactError> {
        let mut workspace = StepWorkspace::new(&self.system).map_err(ImpactError::PreparedSolve)?;
        if let Some(air)=&self.radiation {
            workspace.set_condensed_pairs(&air.state_pairs()).map_err(ImpactError::PreparedSolve)?;
        }
        let candidate = vec![0.0; self.x.len()];
        let output = vec![0.0; self.modes];
        let histories = self.histories.borrow().clone();
        Ok(PreparedImpactSystem { inner: self, workspace, candidate, output, histories, analytic_newton: false })
    }

    /// Prepare the same system using analytic storage tangents for Newton.
    /// Retains the complete Gonzalez energy correction, physical acceptance
    /// gates, and all current state/history. Independent acoustic pairs are
    /// eliminated from the factorization, not the physical equation or state.
    /// No difference probes, split acoustic updates, or new time integrator.
    ///
    /// # Errors
    /// Returns the existing workspace preparation refusal.
    pub fn prepare_analytic(self) -> Result<PreparedImpactSystem, ImpactError> {
        let mut prepared=self.prepare()?;
        prepared.analytic_newton=true;
        Ok(prepared)
    }

    // The reference and prepared images share this exact gate. Constitutive
    // state is immutable during Newton and staged until every check has passed.
    pub(super) fn accept_step(
        &mut self, state: &mut Vec<f64>, candidate: &mut [WoolFeltState],
        before: f64, record: PreparedStepRecord, gate: &CancelGate,
    ) -> Result<ImpactFrame, ImpactError> {
        for (_,offset,film) in &self.membranes { film.observe_interleaved(state,*offset)?; }
        for (_,offset,string) in &self.strings { string.observe_interleaved(state,*offset)?; }
        for (_,offset,support) in &self.supports {support.observe_interleaved(state,*offset)?;}
        let frozen=self.system.hamiltonian(state);let mut crush=0.0;
        candidate.clone_from_slice(&self.histories.borrow());
        for (pad,h) in self.pads.iter().zip(candidate.iter_mut()) {
            let strain=pad.strain(state,self.modes);
            if !strain.is_finite() || strain>pad.spec.law.eps_densify {return Err(invalid("felt trial exceeds densification validity"));}
            let new=pad.spec.law.update_state(strain,h);
            let loss=pad.spec.path_energy(strain,h)-pad.spec.recovered(strain,&new);
            if !loss.is_finite() || loss < -64.0*f64::EPSILON*frozen.abs() {return Err(invalid("felt history update creates energy"));}
            crush+=loss;*h=new;
        }
        let after=frozen-crush;let dissipated=record.dissipated+crush;
        let residual=after-before+dissipated-record.supplied;
        let tolerance=self.config.energy_absolute_tolerance_j+self.config.energy_relative_tolerance
            *(before.abs()+after.abs()+dissipated.abs()+record.supplied.abs());
        if ![after,dissipated,residual,record.supplied,record.solver_residual,tolerance].iter().all(|v|v.is_finite())
            || after<0.0 || after>self.config.maximum_energy_j || record.dissipated<0.0 {
            return Err(invalid("impact candidate exceeds finite energy limits"));
        }
        if residual.abs()>tolerance {return Err(ImpactError::Energy{residual_j:residual,tolerance_j:tolerance});}
        if gate.is_requested() {return Err(ImpactError::Cancelled);}
        std::mem::swap(&mut self.x,state);
        self.histories.borrow_mut().clone_from_slice(candidate);
        self.sample+=1;
        Ok(ImpactFrame{sample:self.sample,time_s:self.sample as f64*self.config.dt_s,stored_energy_j:after,
            dissipated_energy_j:dissipated,felt_crush_loss_j:crush,supplied_work_j:record.supplied,
            balance_residual_j:residual,solver_residual:record.solver_residual})
    }
}

// Read-only reference observations stay available without duplicating their
// implementation. No DerefMut: callers cannot advance behind the workspace.
impl Deref for PreparedImpactSystem {
    type Target = ImpactSystem;
    fn deref(&self) -> &ImpactSystem { &self.inner }
}

impl PreparedImpactSystem {
    /// Select the analytic or original finite-difference Newton Jacobian without
    /// resetting physical state, accepted time, or material history. Both solve
    /// the same equation and use the same physical acceptance gate.
    pub fn set_analytic_newton(&mut self, enabled: bool) { self.analytic_newton=enabled; }

    /// Enable/disable exact acoustic-pair elimination for analytic Newton.
    /// Enabled at preparation when radiation memory exists. Finite-difference
    /// execution ignores the plan. This is a numerical choice: state, accepted
    /// time, force ports, felt and material history remain untouched. Changing
    /// the plan allocates and belongs outside a real-time callback.
    ///
    /// A failed pair solve falls back to full dense LU on the same Jacobian;
    /// no coupling, acoustic pole or Gonzalez correction is dropped.
    pub fn set_radiation_condensation(&mut self, enabled: bool) -> Result<(), ImpactError> {
        let pairs=if enabled {self.inner.radiation.as_ref().map(|air|air.state_pairs()).unwrap_or_default()}
            else {Vec::new()};
        self.workspace.set_condensed_pairs(&pairs).map_err(ImpactError::PreparedSolve)
    }

    /// Planned dense border dimension when using analytic acoustic condensation.
    /// Full state storage and all acoustic histories still exist.
    pub fn condensed_newton_dimension(&self) -> Option<usize> {
        if self.analytic_newton {self.workspace.condensed_dimension()} else {None}
    }

    /// (Condensed, full dense) Newton solves in the most recent solver call.
    /// For substepped output this is the last internal attempt, not a tick sum.
    pub fn newton_linear_solve_counts(&self) -> (usize, usize) {
        self.workspace.linear_solve_counts()
    }

    /// Scalar J-R products in analytic Jacobian assembly for the most recent
    /// solver call. Excludes residuals, Hessian actions, and LU. For internal
    /// substeps this is the last attempt, not an aggregate for the output tick.
    pub fn newton_flow_product_count(&self) -> usize {
        self.workspace.jacobian_flow_product_count()
    }

    /// Resume reference execution without resetting motion, memory or time.
    #[must_use]
    pub fn into_reference(self) -> ImpactSystem { self.inner }

    /// Set the numerical owner's finite Newton-update ceiling.
    /// Zero permits only an already stationary step equation.
    ///
    /// # Errors
    /// Refuses a ceiling above the reference owner's maximum.
    pub fn set_iteration_limit(&mut self, limit: usize) -> Result<(), ImpactError> {
        self.workspace.set_iteration_limit(limit).map_err(ImpactError::PreparedSolve)
    }

    /// Advance the unchanged nonlinear model with held generalized forces.
    ///
    /// Cancellation is polled within Newton/Jacobian work and before physical
    /// publication. Failed trials never advance motion, felt history or time.
    ///
    /// # Errors
    /// The same input, budget, membrane-slope, felt-densification and energy
    /// refusals as the reference, plus typed prepared-solver/cancellation errors.
    pub fn step(&mut self, external: &[f64], gate: &CancelGate) -> Result<ImpactFrame, ImpactError> {
        let inner = &mut self.inner;
        if gate.is_requested() { return Err(ImpactError::Cancelled); }
        if inner.sample >= inner.config.max_steps { return Err(ImpactError::Budget); }
        if external.len() != inner.modes || external.iter().any(|v|
            !v.is_finite() || v.abs() > inner.config.maximum_generalized_force) {
            return Err(invalid("external generalized force shape or ceiling failed"));
        }
        let before = inner.stored_energy_j();
        let ledger = if inner.contact_loss {
            let action=|x:&[f64],d:&[f64],out:&mut[f64]|inner.hessian_vector(x,d,out);
            let tangent=|x:&[f64],e:&[f64],dx:&[f64],de:&[f64],out:&mut[f64]|
                inner.contact.dissipative_flow_tangent_into(x,e,dx,de,out);
            self.workspace.step_into_dissipative_controlled(
                &inner.system,&inner.x,external,inner.config.dt_s,&mut self.candidate,&mut self.output,
                &|x,e,out|inner.contact.dissipative_flow_into(x,e,out),
                if self.analytic_newton {Some(&action)}else{None},
                if self.analytic_newton {Some(&tangent)}else{None},||gate.is_requested(),
            )
        } else if self.analytic_newton {
            self.workspace.step_into_analytic_controlled(
                &inner.system, &inner.x, external, inner.config.dt_s,
                &mut self.candidate, &mut self.output,
                &|x,d,out|inner.hessian_vector(x,d,out), || gate.is_requested(),
            )
        } else {
            self.workspace.step_into_controlled(
                &inner.system, &inner.x, external, inner.config.dt_s,
                &mut self.candidate, &mut self.output, || gate.is_requested(),
            )
        }.map_err(|error| match error {
            PreparedStepError::Cancelled => ImpactError::Cancelled,
            PreparedStepError::Solver(error) => ImpactError::PreparedSolve(error),
        })?;
        inner.accept_step(&mut self.candidate, &mut self.histories, before, ledger, gate)
    }
}

impl ImpactSource for PreparedImpactSystem {
    type Frame = ImpactFrame;
    fn state(&self) -> &[f64] { self.inner.state() }
    fn mode_count(&self) -> usize { self.inner.modes }
    fn samples_rendered(&self) -> u64 { self.inner.sample }
    fn sample_period_s(&self) -> f64 { self.inner.config.dt_s }
    fn remaining_steps(&self) -> u64 { self.inner.config.max_steps - self.inner.sample }
    fn maximum_generalized_force(&self) -> f64 { self.inner.config.maximum_generalized_force }
    fn advance(&mut self, forces: &[f64], gate: &CancelGate) -> Result<Self::Frame, ImpactError> {
        self.step(forces, gate)
    }
}

#[cfg(test)]
#[path = "condensed_tests.rs"]
mod condensed_tests;
