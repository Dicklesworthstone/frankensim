//! Prepared execution of the existing nonlinear impact model.
//!
//! Preparation changes scratch ownership, not the constitutive model or basis.
//! Both execution paths finish through the same physical acceptance gate below.
mod substeps;
mod momentum;
pub use substeps::{ImpactSubstepConfig, ImpactSubstepReport, SubsteppedImpactSystem};
use std::ops::Deref;
use std::cell::RefCell;
use std::rc::Rc;
use fs_exec::CancelGate;
use fs_material::fiber::{Uniaxial, WoolFeltState};
use fs_phs::{PreparedStepError, PreparedStepRecord, StepWorkspace};
use super::{ImpactError, ImpactFrame, ImpactSystem, invalid};
use super::audio::ImpactSource;
use super::tangent::BodyTangents;

/// Persistent nonlinear mechanics with reusable solver and material-trial buffers.
///
/// Construction allocates; stepping does not allocate in this host or its solver.
/// The supplied storage callbacks must independently avoid allocation. This is
/// finite-difference Newton or analytic Newton with exact ordered momentum
/// and acoustic-pair condensation. This is NOT a measured hard-real-time contract.
/// The model, basis, timestep, contact, Kelvin memory and felt conditioning are
/// unchanged. Only accepted steps update physical history and the sample counter.
pub struct PreparedImpactSystem {
    inner: ImpactSystem,
    workspace: StepWorkspace,
    candidate: Vec<f64>,
    output: Vec<f64>,
    histories: Vec<WoolFeltState>,
    analytic_newton: bool,
    body_tangents: RefCell<BodyTangents>,
    batched_body_tangents: bool,
}

impl ImpactSystem {
    /// Prepare this exact system, preserving any already accepted motion/history.
    ///
    /// # Errors
    /// Returns the numerical owner's workspace preparation refusal.
    pub fn prepare(self) -> Result<PreparedImpactSystem, ImpactError> {
        let mut workspace = StepWorkspace::new(&self.system).map_err(ImpactError::PreparedSolve)?;
        momentum::prepare_plan(&self, &mut workspace)?;
        let candidate = vec![0.0; self.x.len()];
        let output = vec![0.0; self.modes];
        let histories = self.histories.borrow().clone();
        let body_tangents=RefCell::new(BodyTangents::new(Rc::clone(&self.mechanical)));
        Ok(PreparedImpactSystem { inner: self, workspace, candidate, output, histories,
            analytic_newton: false, body_tangents, batched_body_tangents: true })
    }

    /// Prepare the same system using analytic storage tangents for Newton.
    /// Retains the complete Gonzalez energy correction, physical acceptance
    /// gates, and all current state/history. Canonical momenta and acoustic pairs
    /// are eliminated from the factorization, not the physical equation or state.
    /// Ordered kinematic pivots retain cross-velocity damping, nonlinear contact
    /// loss and reciprocal acoustic coupling; displacements/material memories
    /// remain in the dense border. Simple unloaded models keep scalar leaves.
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
        self.validate_material_contacts(state,state)?;
        self.mechanical.validate_rigid_motion(state)?;
        for (_,offset,film) in &self.membranes { film.observe_interleaved(state,*offset)?; }
        for (_,offset,string) in &self.strings { string.observe_interleaved(state,*offset)?; }
        for (_,offset,support) in &self.supports {support.observe_interleaved(state,*offset)?;}
        self.validate_squeeze_configuration(state)?;
        if let Some(memory)=&self.relaxation {memory.admit_state(state)?;}
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

    /// Reuse each shell/film's exact tangent across analytic Newton columns.
    /// Enabled at preparation. Only identical body displacement bits reuse a
    /// matrix; no time extrapolation or approximate stiffness is introduced.
    /// Contact, felt, material-memory and acoustic terms are evaluated live.
    /// Disabling selects the original per-direction owner actions for comparison.
    /// This switch neither allocates nor changes physical state/history/time.
    /// Finite-difference stepping ignores it. Rounding/iteration counts may differ.
    pub fn set_batched_body_tangents(&mut self, enabled: bool) { self.batched_body_tangents=enabled; }

    /// (Body matrix assemblies, body matrix actions) in the most recent solver
    /// attempt. Excludes contact, felt, acoustic work, residuals and LU; these
    /// counts are not timing. Substepped output reports its last internal attempt.
    /// Counts saturate at usize::MAX; state/history acceptance is independent.
    pub fn body_tangent_work_counts(&self) -> (usize, usize) { self.body_tangents.borrow().counts() }

    /// Enable/disable joint kinematic/acoustic elimination for analytic Newton.
    /// Enabled at preparation when radiation memory exists. Finite-difference
    /// execution ignores the plan. This is a numerical choice: state, accepted
    /// time, force ports, felt and material history remain untouched. Changing
    /// the plan allocates and belongs outside a real-time callback.
    ///
    /// This replaces the current elimination plan; false selects full dense LU.
    /// On a model without radiation, use `set_momentum_condensation` to restore
    /// the unloaded momentum plan. Loaded models with at least two canonical
    /// modes eliminate momenta using displacement-equation pivots FIRST, then
    /// solve the acoustic pairs with their lower couplings retained. One-mode
    /// models retain the independent-acoustic-pair plan. The two physical sets
    /// are never treated as mutually independent. Failed elimination falls back
    /// to full dense LU on the same Jacobian;
    /// no coupling, acoustic pole or Gonzalez correction is dropped.
    pub fn set_radiation_condensation(&mut self, enabled: bool) -> Result<(), ImpactError> {
        if enabled && self.inner.radiation.is_some() {
            momentum::prepare_plan(&self.inner,&mut self.workspace)
        } else {
            self.workspace.set_condensed_pairs(&[]).map_err(ImpactError::PreparedSolve)
        }
    }

    /// Planned dense border dimension for the selected analytic elimination.
    /// Includes the scalar Gonzalez correction. Full physical state and all
    /// material/acoustic histories still exist; this is not a timing claim.
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
        self.step_checked(external, gate, |_,_,_,_,_|Ok(()))
    }

    // Stage dependent physical states before our common acceptance gate.
    // The callback receives the actual generalized port velocity AND complete
    // discrete effort used by fs-phs, including the Gonzalez correction. It
    // must not mutate accepted histories.
    pub(super) fn step_checked(&mut self, external: &[f64], gate: &CancelGate,
        check: impl FnOnce(&ImpactSystem, &[f64], &[f64], &[f64], PreparedStepRecord) -> Result<(),ImpactError>,
    ) -> Result<ImpactFrame, ImpactError> {
        let inner = &mut self.inner;
        if gate.is_requested() { return Err(ImpactError::Cancelled); }
        if inner.sample >= inner.config.max_steps { return Err(ImpactError::Budget); }
        if external.len() != inner.modes || external.iter().any(|v|
            !v.is_finite() || v.abs() > inner.config.maximum_generalized_force) {
            return Err(invalid("external generalized force shape or ceiling failed"));
        }
        let before = inner.stored_energy_j();
        let cache=&self.body_tangents;
        cache.borrow_mut().reset_counts();
        let batched=self.batched_body_tangents;
        let action=|x:&[f64],d:&[f64],out:&mut[f64]| {
            if batched {inner.hessian_vector_cached(x,d,out,&mut cache.borrow_mut())}
            else {inner.hessian_vector(x,d,out)}
        };
        let ledger = if inner.has_nonlinear_dissipation() {
            let tangent=|x:&[f64],e:&[f64],dx:&[f64],de:&[f64],out:&mut[f64]|
                inner.dissipative_flow_tangent_into(x,e,dx,de,out);
            self.workspace.step_into_dissipative_controlled(
                &inner.system,&inner.x,external,inner.config.dt_s,&mut self.candidate,&mut self.output,
                &|x,e,out|inner.dissipative_flow_into(x,e,out),
                if self.analytic_newton {Some(&action)}else{None},
                if self.analytic_newton {Some(&tangent)}else{None},||gate.is_requested(),
            )
        } else if self.analytic_newton {
            self.workspace.step_into_analytic_controlled(
                &inner.system, &inner.x, external, inner.config.dt_s,
                &mut self.candidate, &mut self.output,
                &action, || gate.is_requested(),
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
        let effort = self.workspace.last_step_discrete_effort()
            .ok_or_else(|| invalid("successful impact step is missing its discrete effort"))?;
        check(inner, &self.candidate, &self.output, effort, ledger)?;
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

impl PreparedImpactSystem {
    /// Replace unloaded source contacts without rebuilding or resetting prepared
    /// state, memory, caches or solver controls. See the reference owner's gate.
    pub fn refresh_unloaded_contacts(&mut self,
        updates:Vec<super::material::refresh::ContactSourceUpdate>,gate:&CancelGate,
    )->Result<usize,ImpactError> {
        self.inner.refresh_unloaded_contacts(updates,gate)
    }

    pub(crate) fn swap_contact_refresh(&mut self,next:super::material::refresh::ContactRefresh)
        ->super::material::refresh::ContactRefresh {
        self.inner.swap_contact_refresh(next)
    }

    // Only admitted existing islands may be replaced at an exclusive accepted
    // boundary. The caller stages every coefficient before this infallible
    // publication; no mechanical coordinate or prepared workspace is rebuilt.
    pub(crate) fn swap_squeeze_films(
        &mut self, updates: &mut [(usize, super::squeeze::ResistiveFilm)],
    ) {
        for (index, film) in updates {
            std::mem::swap(&mut self.inner.squeeze_films[*index], film);
        }
    }
}

