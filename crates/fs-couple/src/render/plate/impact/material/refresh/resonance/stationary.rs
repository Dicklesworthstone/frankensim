//! Force-loaded operating points of the original impact/film/material source.
//! The shared Newton--Krylov bridge owns the solve; this module only binds the
//! same source admission, Hessian and resisting laws used by local analysis.
use super::*;
use crate::stationary_equilibrium::{
    EquilibriumSolveError, EquilibriumSolveOptions, SolvedEquilibrium, solve_equilibrium,
};
use fs_phs::harmonic::equilibrium::stationary::{StationaryCoordinates, StationaryProblem};
use std::cell::RefCell;

impl ImpactSystem {
    // One gate for caller-resolved and solver-resolved stationary configurations.
    // Caller scratch avoids a fresh effort allocation in every line-search trial.
    pub(super) fn admit_resonance_state(
        &self, state: &[f64], input: &[f64], effort: &mut [f64],
    ) -> Result<(), ImpactError> {
        if state.len()!=self.x.len() || effort.len()!=state.len() || input.len()!=self.modes
            || state.iter().chain(input).any(|x|!x.is_finite())
            || input.iter().any(|f|f.abs()>self.config.maximum_generalized_force) {
            return Err(invalid("resonance point requires complete finite coordinates and admitted force"));
        }
        if !self.pads.is_empty() {
            return Err(invalid("irreversible felt histories need a smooth explicit-state linearization"));
        }
        self.mechanical.validate_rigid_motion(state)?;
        for (_,offset,film) in &self.membranes {film.observe_interleaved(state,*offset)?;}
        for (_,offset,string) in &self.strings {string.observe_interleaved(state,*offset)?;}
        for (_,offset,support) in &self.supports {support.observe_interleaved(state,*offset)?;}
        self.validate_squeeze_configuration(state)?;
        let energy=self.system.hamiltonian(state);
        if !energy.is_finite() || energy<0. || energy>self.config.maximum_energy_j {
            return Err(invalid("resonance point exceeds original physical energy admission"));
        }
        effort.fill(f64::NAN);
        self.system.storage().gradient(state,effort);
        if effort.iter().any(|x|!x.is_finite()) {
            return Err(invalid("original resonance-point effort is not representable"));
        }
        self.validate_material_contacts(state,effort)?;
        if !self.contact.has_two_sided_contact_tangent(state,effort) {
            return Err(invalid("contact point lacks a two-sided local storage/loss derivative"));
        }
        Ok(())
    }

    /// Resolve a held generalized force into the complete original stationary
    /// configuration and its local resonances. Coordinates use state(), including
    /// acoustic and material histories; scales are explicit physical caller data.
    ///
    /// Seed selection chooses a local branch, not a global minimum or unique
    /// root. All original admissibility limits apply to each proposed point;
    /// source-domain trials can backtrack without clipping penetration or gaps.
    /// Every material memory is solved in place in the equations, not reset to
    /// an assumed history. A nonstationary free load or unresolved mechanism can
    /// refuse. No diagonal stiffness or damping is substituted to make it settle.
    ///
    /// This operation is read-only: the result is NOT installed into the accepted
    /// mechanical trajectory. Samples, thermodynamic/material histories and the
    /// time integrator remain untouched. The existing resonance no-claims for
    /// felt, derivative corners, finite amplitude and contact switching remain.
    pub fn solve_equilibrium(
        &self, initial_state: &[f64], coordinates: StationaryCoordinates<'_>,
        budget: LinearizationBudget, max_loss_tangent_actions: usize,
        options: EquilibriumSolveOptions, mut cancelled: impl FnMut()->bool,
    ) -> Result<SolvedEquilibrium, ImpactError> {
        if cancelled() {return Err(ImpactError::Cancelled);}
        let n=self.x.len();
        if n>fs_phs::harmonic::MAX_HARMONIC_DIMENSION || n>budget.harmonic.max_dimension
            || n>budget.max_hessian_actions || initial_state.len()!=n {
            return Err(invalid("complete mechanical resonance state exceeds its declared capacity"));
        }
        let has_loss=self.has_nonlinear_dissipation();
        if has_loss && n>max_loss_tangent_actions {
            return Err(invalid("complete nonlinear-loss tangent exceeds its action capacity"));
        }
        let mut scratch=vec![f64::NAN;n];
        self.admit_resonance_state(initial_state,coordinates.input,&mut scratch)?;
        if cancelled() {return Err(ImpactError::Cancelled);}
        let scratch=RefCell::new(scratch);
        let domain=|x:&[f64]|self.admit_resonance_state(x,coordinates.input,&mut scratch.borrow_mut()).is_ok();
        let hessian=|x:&[f64],d:&[f64],out:&mut[f64]|self.hessian_vector(x,d,out);
        let flow=|x:&[f64],e:&[f64],out:&mut[f64]|self.dissipative_flow_into(x,e,out);
        let tangent=|x:&[f64],e:&[f64],dx:&[f64],de:&[f64],out:&mut[f64]|
            self.dissipative_flow_tangent_into(x,e,dx,de,out);
        let loss=has_loss.then_some(NonlinearDissipation {
            flow:&flow,tangent:&tangent,max_tangent_actions:max_loss_tangent_actions,
        });
        let problem=StationaryProblem::from_ode(&self.system,coordinates,budget,
            &hessian,loss,&domain).map_err(|e|ImpactError::Owner(e.to_string()))?;
        solve_equilibrium(&problem,initial_state,options,&mut cancelled).map_err(|error|match error {
            EquilibriumSolveError::Cancelled=>ImpactError::Cancelled,
            other=>ImpactError::Owner(other.to_string()),
        })
    }
}
