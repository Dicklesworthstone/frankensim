//! Finite isothermal reservoirs and sorbed film storage in one implicit solve.
//!
//! Each reservoir explicitly supplies its molar pressure capacity C [mol/Pa],
//! so n = C p. No equation of state is inferred here. The same tridiagonal
//! finite-volume matrix gains two endpoint storage rows: C_L dp_L/dt = -A J_L
//! and C_R dp_R/dt = A J_R. Every membrane face therefore has one transfer,
//! while the entire reservoir/film/reservoir system retains its total amount.
//! At connected equilibrium p = n_total/(C_L+C_R+A sum(S_i L_i)).
//!
//! Temperature and capacities stay fixed. Maintaining isothermal conditions
//! requires an external thermostat; conservation here is of amount, not energy.

use super::*;
use fs_qty::Dims;

pub const MOLAR_PRESSURE_CAPACITY_DIMS:Dims=Dims([1,-1,2,0,0,1]);

/// Whole-reservoir dn/dp [mol/Pa] at its declared fixed external conditions.
/// A pure ideal gas in a rigid isothermal volume has C=V/(R T); that EOS and
/// its source qualification belong to the caller, not this numerical owner.
#[derive(Clone,Copy,Debug,PartialEq)]
pub struct MolarPressureCapacity(f64);
impl MolarPressureCapacity {
    pub fn try_new(value_mol_pa:f64)->Result<Self,MembraneDiffusionError> {
        if !positive(value_mol_pa) {return Err(Error::Invalid("finite reservoirs require positive representable molar pressure capacities"));}
        Ok(Self(value_mol_pa))
    }
    pub const fn value_si(self)->f64 {self.0}
}

#[derive(Clone,Debug,PartialEq)]
pub struct FiniteMembraneDiffusionState {
    membrane:MembraneDiffusionState,
    capacities:[MolarPressureCapacity;2],
    reservoir_pressures_pa:[f64;2],
}

/// Retained factors for the original grid plus two finite reservoir rows.
/// Factor reuse binds the exact membrane model, both capacities and duration.
#[derive(Debug)]
pub struct FiniteMembraneDiffusionWorkspace {
    inner:MembraneDiffusionWorkspace,
}
impl FiniteMembraneDiffusionWorkspace {
    /// `cells` counts membrane cells; the two reservoir rows are added here.
    pub fn new(cells:usize)->Result<Self,MembraneDiffusionError> {
        if cells==0||cells>MAX_TRANSIENT_MEMBRANE_CELLS {
            return Err(Error::Invalid("bounded nonempty finite-reservoir membrane workspace required"));
        }
        Ok(Self {inner:MembraneDiffusionWorkspace::with_rows(cells+2)?})
    }
}

#[derive(Clone,Debug,PartialEq)]
pub struct FiniteMembraneDiffusionStep {
    /// The original face-flux, sorbed-storage and balance report. Its reservoir
    /// pressures are the jointly solved END pressures used by backward Euler.
    pub membrane:MembraneDiffusionStep,
    pub reservoir_pressures_before_pa:[f64;2],
    pub reservoir_pressures_after_pa:[f64;2],
    pub reservoir_amounts_before_mol:[f64;2],
    pub reservoir_amounts_after_mol:[f64;2],
    /// C*(p_new-p_old), avoiding subtraction of two large reservoir inventories.
    pub reservoir_amount_changes_mol:[f64;2],
    pub reservoir_balance_residuals_mol:[f64;2],
    pub total_before_mol:f64,
    pub total_after_mol:f64,
    /// Sum of both reservoir changes and the membrane storage change.
    pub amount_residual_mol:f64,
}

#[derive(Clone,Debug,PartialEq)]
pub struct FiniteMembraneDiffusionCandidate {
    state:FiniteMembraneDiffusionState,
    step:FiniteMembraneDiffusionStep,
}
impl FiniteMembraneDiffusionCandidate {
    pub const fn state(&self)->&FiniteMembraneDiffusionState {&self.state}
    pub const fn step(&self)->&FiniteMembraneDiffusionStep {&self.step}
    pub fn into_parts(self)->(FiniteMembraneDiffusionState,FiniteMembraneDiffusionStep) {(self.state,self.step)}
}

impl FiniteMembraneDiffusionState {
    /// Attach two finite inventories to an existing explicitly initialized
    /// sorbed state. Zero pressure is a valid empty numerical reservoir.
    pub fn new(membrane:MembraneDiffusionState,capacities:[MolarPressureCapacity;2],
        reservoir_pressures_pa:[f64;2],mut cancelled:impl FnMut()->bool)->Result<Self,MembraneDiffusionError>
    {
        poll(&mut cancelled)?;
        if !reservoir_pressures_pa.into_iter().all(nonnegative) {
            return Err(Error::Invalid("finite reservoir pressures must be nonnegative and representable"));
        }
        let state=Self {membrane,capacities,reservoir_pressures_pa};
        state.total_amount_mol()?;
        poll(&mut cancelled)?;
        Ok(state)
    }
    pub const fn membrane(&self)->&MembraneDiffusionState {&self.membrane}
    /// Move the sorbed state into its source-bound owner without cloning cells.
    pub fn into_membrane(self)->MembraneDiffusionState {self.membrane}
    pub const fn capacities(&self)->[MolarPressureCapacity;2] {self.capacities}
    pub const fn reservoir_pressures_pa(&self)->[f64;2] {self.reservoir_pressures_pa}
    pub fn reservoir_amounts_mol(&self)->Result<[f64;2],MembraneDiffusionError> {
        let mut amounts=[0.0;2];
        for i in 0..2 {
            amounts[i]=finite(self.capacities[i].0*self.reservoir_pressures_pa[i])?;
            if self.reservoir_pressures_pa[i]>0.0&&amounts[i]==0.0 {
                return Err(Error::Invalid("finite reservoir amount underflowed"));
            }
        }
        Ok(amounts)
    }
    pub fn total_amount_mol(&self)->Result<f64,MembraneDiffusionError> {
        let amounts=self.reservoir_amounts_mol()?;
        let mut total=Sum::default();
        total.add(amounts[0])?;total.add(self.membrane.stored_amount_mol()?)?;total.add(amounts[1])?;
        Ok(total.value)
    }

    /// Solve both reservoirs and every sorbed cell together. There is no
    /// prescribed-pressure predictor followed by a separate inventory debit.
    /// All original cell, reservoir and whole-system balances must pass before
    /// an immutable candidate is returned. Equal-size steps retain O(n) solves
    /// with cached factors; stability does not remove time/grid accuracy error.
    pub fn step_candidate(&self,duration_s:f64,controls:MembraneDiffusionSolve,
        workspace:&mut FiniteMembraneDiffusionWorkspace,mut cancelled:impl FnMut()->bool)
        ->Result<FiniteMembraneDiffusionCandidate,MembraneDiffusionError>
    {
        poll(&mut cancelled)?;controls.validate()?;
        let elapsed=self.membrane.elapsed_s+duration_s;
        if !positive(duration_s)||!elapsed.is_finite()||elapsed<=self.membrane.elapsed_s {
            return Err(Error::Invalid("positive representable finite-reservoir duration required"));
        }
        let work=&mut workspace.inner;
        let n=self.membrane.partial_pressures_pa.len();
        work.factor_with_reservoirs(&self.membrane.model,duration_s,
            Some(self.capacities.map(MolarPressureCapacity::value_si)),&mut cancelled)?;
        for i in 0..n+2 {
            if i%128==0 {poll(&mut cancelled)?;}
            let pressure=if i==0 {self.reservoir_pressures_pa[0]}else if i==n+1 {self.reservoir_pressures_pa[1]}
                else {self.membrane.partial_pressures_pa[i-1]};
            work.rhs[i]=finite(work.mass[i]*pressure)?;
        }
        work.linear.solve_factored_into(&work.rhs,&mut work.solution,&mut cancelled)?;
        let pressures=[work.solution[0],work.solution[n+1]];
        let candidate=self.membrane.finish_step(pressures,duration_s,&work.solution[1..n+1],controls,&mut cancelled)?;
        let (membrane,report)=candidate.into_parts();
        let state=Self::new(membrane,self.capacities,pressures,&mut cancelled)?;
        let before=self.reservoir_amounts_mol()?;let after=state.reservoir_amounts_mol()?;
        let mut changes=[0.0;2];let mut residuals=[0.0;2];
        let outward=[report.left_reservoir_to_membrane_mol,-report.membrane_to_right_reservoir_mol];
        let mut total_change=Sum::default();
        for i in 0..2 {
            poll(&mut cancelled)?;
            changes[i]=finite(self.capacities[i].0*(pressures[i]-self.reservoir_pressures_pa[i]))?;
            residuals[i]=finite(changes[i]+outward[i])?;
            let tolerance=controls.tolerance(changes[i].abs()+outward[i].abs())?;
            if residuals[i].abs()>tolerance {
                return Err(Error::AmountBalance {residual_mol:residuals[i],tolerance_mol:tolerance});
            }
            total_change.add(changes[i])?;
        }
        total_change.add(report.stored_change_mol)?;
        let tolerance=controls.tolerance(changes[0].abs()+changes[1].abs()+report.stored_change_mol.abs())?;
        if total_change.value.abs()>tolerance {
            return Err(Error::AmountBalance {residual_mol:total_change.value,tolerance_mol:tolerance});
        }
        let step=FiniteMembraneDiffusionStep {membrane:report,
            reservoir_pressures_before_pa:self.reservoir_pressures_pa,reservoir_pressures_after_pa:pressures,
            reservoir_amounts_before_mol:before,reservoir_amounts_after_mol:after,
            reservoir_amount_changes_mol:changes,reservoir_balance_residuals_mol:residuals,
            total_before_mol:self.total_amount_mol()?,total_after_mol:state.total_amount_mol()?,
            amount_residual_mol:total_change.value};
        poll(&mut cancelled)?;
        Ok(FiniteMembraneDiffusionCandidate {state,step})
    }
}

#[cfg(test)]
#[path="finite/tests.rs"]
mod tests;
