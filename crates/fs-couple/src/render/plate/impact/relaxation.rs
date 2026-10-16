//! Material memory composed through fs-phs, not an endpoint damping correction.
//! The equilibrium solid remains in the original mechanical storage. Additional
//! Maxwell arms use the existing owner's internal-strain storage and resistance.
use super::{ImpactError, ImpactSystem, invalid};
use fs_math::det;
use fs_phs::RelaxationBranch;

#[path = "shell_relaxation.rs"]
mod shell;
pub use shell::ShellBendingSpectrum;

/// Initial viscous strain is a physical initial condition, not a solver choice.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum InitialMemory {
    /// Each added arm is equilibrated at the declared initial displacement.
    Relaxed,
    /// Zero prior viscous strain; initial displacement also loads every arm.
    Unrelaxed,
}

/// Read-only endpoint quantities. This energy is ALREADY in total storage.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct RelaxationObservation {
    pub states: usize,
    pub stored_energy_j: f64,
    /// Instantaneous power, not the accepted interval's dissipated energy.
    pub dissipated_power_w: f64,
}

// Only derivative/observation metadata. fs-phs alone owns the actual storage,
// resistance and time equation; retained rows use its exact sqrt(k) scaling.
pub(super) struct Memory {
    pub(super) base_dim: usize,
    rows: Vec<Vec<f64>>,
    rates: Vec<f64>,
}
impl Memory {
    pub(super) fn add_hessian(&self, direction: &[f64], output: &mut [f64]) {
        output[self.base_dim..].fill(0.0);
        for (i, row) in self.rows.iter().enumerate() {
            let change = row.iter().zip(direction).map(|(b,d)| b*d).sum::<f64>()
                - direction[self.base_dim+i];
            for (out,b) in output[..self.base_dim].iter_mut().zip(row) { *out += b*change; }
            output[self.base_dim+i] = -change;
        }
    }
    fn observe(&self, state: &[f64]) -> RelaxationObservation {
        let mut result = RelaxationObservation { states:self.rows.len(), ..Default::default() };
        for (i,(row,rate)) in self.rows.iter().zip(&self.rates).enumerate() {
            let effort = row.iter().zip(state).map(|(b,x)| b*x).sum::<f64>() - state[self.base_dim+i];
            result.stored_energy_j += 0.5*effort*effort;
            result.dissipated_power_w += rate*effort*effort;
        }
        result
    }
}
impl ImpactSystem {
    /// Append declared hereditary material arms BEFORE starting the instrument.
    ///
    /// Reuses `fs_phs::PortHamiltonian::with_relaxation_branches` verbatim:
    /// `H_arm=(sqrt(k)*projection.x-z)^2/2`, `R_zz=1/tau`. No new integrator,
    /// split force, or independent energy debit. Existing equilibrium stiffness
    /// must NOT already include the relaxing branch stiffness. Projection uses
    /// the COMPLETE current state layout but may load only mechanical q entries;
    /// momenta, felt/Kelvin histories, and external force ports are unchanged.
    ///
    /// The caller declares a branch ceiling (at most 256); exhausted capacity
    /// refuses rather than dropping material modes. One scalar memory follows
    /// the original state per arm. Preparation, analytic Newton, contact loss,
    /// cancellation and substep rollback all retain that full state. Initial
    /// memory is explicit, and its energy is subject to the original ceiling.
    ///
    /// # Errors
    /// Refuses invalid projections/coefficients/budgets, repeated attachment,
    /// or attachment after motion has started. This is not an in-flight material
    /// replacement with work accounting. Empty arms preserve the original image.
    pub fn with_relaxation_branches(mut self, branches: Vec<RelaxationBranch>,
        initial: InitialMemory, maximum_branches: usize) -> Result<Self, ImpactError>
    {
        if branches.is_empty() { return Ok(self); }
        if maximum_branches>256 || branches.len()>maximum_branches || self.sample!=0
            || self.relaxation.is_some() {
            return Err(invalid("material relaxation requires one cold attachment within its declared 256-arm ceiling"));
        }
        let base_dim=self.x.len();
        let mut rows=Vec::with_capacity(branches.len());
        let mut rates=Vec::with_capacity(branches.len());
        let mut memories=Vec::with_capacity(branches.len());
        for branch in &branches {
            if branch.projection.len()!=base_dim || branch.projection.iter().enumerate().any(|(i,b)|
                !b.is_finite() || *b!=0.0 && (i>=2*self.modes || i%2!=0))
                || !branch.stiffness.is_finite() || branch.stiffness<=0.0
                || !branch.relaxation_time_s.is_finite() || branch.relaxation_time_s<=0.0 {
                return Err(invalid("relaxation needs finite displacement-only projection, positive stiffness and positive time"));
            }
            let scale=det::sqrt(branch.stiffness);
            let row:Vec<_>=branch.projection.iter().map(|b| b*scale).collect();
            let rate=1.0/branch.relaxation_time_s;
            if row.iter().any(|b|!b.is_finite()) || !row.iter().any(|b|*b!=0.0)
                || !rate.is_finite() || rate<=0.0 {
                return Err(invalid("unrepresentable relaxation arm"));
            }
            let z:f64=match initial {
                InitialMemory::Relaxed=>row.iter().zip(&self.x).map(|(b,x)|b*x).sum(),
                InitialMemory::Unrelaxed=>0.0,
            };
            if !z.is_finite() {return Err(invalid("initial material memory overflow"));}
            rows.push(row);rates.push(rate);memories.push(z);
        }
        self.system=self.system.with_relaxation_branches(branches).map_err(ImpactError::PreparedSolve)?;
        self.x.extend(memories);
        let energy=self.stored_energy_j();
        if !energy.is_finite() || energy<0.0 || energy>self.config.maximum_energy_j {
            return Err(invalid("initial material memory exceeds the mechanical energy ceiling"));
        }
        self.relaxation=Some(Memory {base_dim,rows,rates});
        Ok(self)
    }

    /// Observe added material memory without advancing mechanics or its history.
    pub fn relaxation_observation(&self) -> Option<RelaxationObservation> {
        self.relaxation.as_ref().map(|memory|memory.observe(&self.x))
    }
}

#[cfg(test)]
#[path = "relaxation_tests.rs"]
mod tests;
