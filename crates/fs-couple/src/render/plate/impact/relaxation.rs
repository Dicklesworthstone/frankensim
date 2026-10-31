//! Material memory composed through fs-phs, not an endpoint damping correction.
//! The equilibrium solid remains in the original mechanical storage. Additional
//! Maxwell arms use the existing owner's internal-strain storage and resistance.
use super::{ImpactError, ImpactSystem, invalid};
use super::radiation::MAX_JOINT_STATES;
use fs_math::det;
use fs_phs::{RelaxationBranch,MappedRelaxationBranch};

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
    mapped: Vec<MappedRelaxationBranch>,
    supports: Vec<(usize,super::supported::SupportedStrings,usize,usize)>,
    films: Vec<(usize,super::membrane::MembranePotential,usize,usize)>,
    joint_omega_bound: f64,
}
impl Memory {
    pub(crate) fn branch_count(&self) -> usize { self.rates.len() }

    /// Original diagonal Maxwell resistance at the time owner's FULL discrete
    /// effort. Appended acoustic coordinates are allowed; missing material
    /// coordinates and unrepresentable nonzero power are refused.
    pub(crate) fn branch_dissipated_power(&self, discrete_effort: &[f64], branch: usize)
        -> Result<f64, ImpactError>
    {
        if discrete_effort.len() < self.base_dim + self.rates.len() {
            return Err(invalid("material heat requires the complete discrete effort"));
        }
        let rate = *self.rates.get(branch)
            .ok_or_else(|| invalid("material heat names an absent relaxation arm"))?;
        let effort = discrete_effort[self.base_dim + branch];
        // Keep the original e^T R e multiplication order, not an inferred
        // endpoint power or a mechanical energy deficit.
        let power = effort * (rate * effort);
        if !effort.is_finite() || !power.is_finite() || power < 0.0
            || (effort != 0.0 && power == 0.0) {
            return Err(invalid("material relaxation heat is unrepresentable"));
        }
        Ok(power)
    }

    pub(super) fn add_hessian(&self, state:&[f64], direction: &[f64], output: &mut [f64]) {
        output[self.base_dim..].fill(0.0);
        for (i, row) in self.rows.iter().enumerate() {
            let change = row.iter().zip(direction).map(|(b,d)| b*d).sum::<f64>()
                - direction[self.base_dim+i];
            for (out,b) in output[..self.base_dim].iter_mut().zip(row) { *out += b*change; }
            output[self.base_dim+i] = -change;
        }
        for (i,branch) in self.mapped.iter().enumerate() {
            let scale=det::sqrt(branch.stiffness);
            let x=&state[..self.base_dim];let d=&direction[..self.base_dim];
            let effort=scale*branch.strain.strain(x)-state[self.base_dim+i];
            let mut g=vec![0.0;self.base_dim];let mut h=vec![0.0;self.base_dim];
            branch.strain.gradient(x,&mut g);branch.strain.hessian_vector(x,d,&mut h);
            for value in &mut g {*value*=scale;}
            let change=g.iter().zip(d).map(|(g,d)|g*d).sum::<f64>()-direction[self.base_dim+i];
            for ((out,g),h) in output[..self.base_dim].iter_mut().zip(g).zip(h) {
                *out+=g*change+effort*(scale*h);
            }
            output[self.base_dim+i]=-change;
        }
    }
    pub(super) fn admit_state(&self,state:&[f64])->Result<(),ImpactError> {
        for (offset,support,start,count) in &self.supports {
            let memory=state.get(self.base_dim+start..self.base_dim+start+count)
                .ok_or_else(||invalid("supported material history prefix is incomplete"))?;
            support.admit_memory(*offset,state,memory)?;
        }
        for (offset,film,start,count) in &self.films {
            let memory=state.get(self.base_dim+start..self.base_dim+start+count)
                .ok_or_else(||invalid("film material history prefix is incomplete"))?;
            film.admit_memory(*offset,state,memory)?;
        }
        Ok(())
    }
    pub(super) fn admit_radiation(&self,omega:f64,coupling_bound:f64,dt:f64)
        -> Result<Option<f64>,ImpactError>
    {
        if self.supports.is_empty() && self.films.is_empty() {return Ok(None);}
        // In energy-scaled canonical coordinates, the old and acoustic
        // oscillator blocks have norm max(Omega_old,Omega_air). The reciprocal
        // velocity-port skew block adds at most ||L||_2, bounded by the original
        // radiation owner's max absolute row/column sum. Memory has no J row.
        let bound=self.joint_omega_bound.max(omega)+coupling_bound;
        if !bound.is_finite() {return Err(invalid("supported hereditary radiation bound overflows"));}
        fs_phs::linear_exact::admit_midpoint_band(bound,dt).map_err(ImpactError::TimeResolution)?;
        for (_,support,_,_) in &self.supports {
            if let Some(config)=support.memory_configuration() {config.admit_domain(bound,dt)?;}
        }
        for (_,film,_,_) in &self.films {
            if let Some(config)=film.memory_configuration() {config.admit_domain(bound,dt)?;}
        }
        Ok(Some(bound))
    }
    pub(super) fn commit_radiation_bound(&mut self,bound:f64) {self.joint_omega_bound=bound;}
    fn observe(&self, state: &[f64]) -> RelaxationObservation {
        let mut result = RelaxationObservation { states:self.rates.len(), ..Default::default() };
        for (i,(row,rate)) in self.rows.iter().zip(&self.rates).enumerate() {
            let effort = row.iter().zip(state).map(|(b,x)| b*x).sum::<f64>() - state[self.base_dim+i];
            result.stored_energy_j += 0.5*effort*effort;
            result.dissipated_power_w += rate*effort*effort;
        }
        for (i,(branch,rate)) in self.mapped.iter().zip(&self.rates).enumerate() {
            let effort=det::sqrt(branch.stiffness)*branch.strain.strain(&state[..self.base_dim])-state[self.base_dim+i];
            result.stored_energy_j+=0.5*effort*effort;
            result.dissipated_power_w+=rate*effort*effort;
        }
        result
    }
}
impl ImpactSystem {
    /// A later load needs a new complete strain/history and source-band bound.
    /// Preserve the old linear-memory composition path; the mapped supported
    /// path must not silently reuse authority from its previous operator.
    pub(super) fn require_recomposition_without_mapped_memory(&self)->Result<(),ImpactError> {
        if self.relaxation.as_ref().is_some_and(|memory|!memory.supports.is_empty() || !memory.films.is_empty()) {
            return Err(invalid("nonlinear gap loading needs complete strain/history domain recomposition"));
        }
        Ok(())
    }

    /// Original attached Maxwell arms, in their unchanged insertion order.
    /// Each arm has one physical loss channel; zero means no memory attachment.
    pub fn relaxation_branch_count(&self) -> usize {
        self.relaxation.as_ref().map_or(0, Memory::branch_count)
    }

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
    /// attachment after motion has started, or memory beyond the joint
    /// `MAX_JOINT_STATES` envelope (states already attached included). This is
    /// not an in-flight material replacement with work accounting. Empty arms
    /// preserve the original image.
    pub fn with_relaxation_branches(mut self, branches: Vec<RelaxationBranch>,
        initial: InitialMemory, maximum_branches: usize) -> Result<Self, ImpactError>
    {
        if branches.is_empty() { return Ok(self); }
        if maximum_branches>256 || branches.len()>maximum_branches || self.sample!=0
            || self.relaxation.is_some() {
            return Err(invalid("material relaxation requires one cold attachment within its declared 256-arm ceiling"));
        }
        let base_dim=self.x.len();
        if base_dim+branches.len()>MAX_JOINT_STATES {
            return Err(invalid("material relaxation exceeds the joint 1024-state work envelope"));
        }
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
        self.relaxation=Some(Memory {base_dim,rows,rates,mapped:vec![],supports:vec![],films:vec![],joint_omega_bound:0.0});
        Ok(self)
    }

    /// Called by the original constructor after all body/felt state is known.
    /// Immutable supported-string and film data cannot be dropped when composed
    /// with an impact exciter. Memory remains in the one existing pHS state and ledger.
    pub(super) fn with_body_relaxation(mut self,joint_omega_bound:f64)->Result<Self,ImpactError> {
        let base_dim=self.x.len();let mut branches=Vec::new();let mut supports=Vec::new();
        let mut memories=Vec::new();let mut rates=Vec::new();let mut films=Vec::new();
        for (_,offset,support) in &self.supports {
            let Some(config)=support.memory_configuration() else {continue;};
            config.admit_domain(joint_omega_bound,self.config.dt_s)?;
            let local=support.memory_branches(*offset,base_dim)?;
            let start=branches.len();
            if start+local.len()>256 || base_dim+start+local.len()>MAX_JOINT_STATES {
                return Err(invalid("complete supported material histories exceed the joint state ceiling"));
            }
            for branch in &local {
                let rate=1.0/branch.relaxation_time_s;
                if !rate.is_finite() || rate<=0.0 || self.config.dt_s*rate>2.0 {
                    return Err(invalid("original impact clock does not resolve supported material relaxation"));
                }
                let z=match config.initial {
                    InitialMemory::Relaxed=>det::sqrt(branch.stiffness)*branch.strain.strain(&self.x),
                    InitialMemory::Unrelaxed=>0.0,
                };
                if !z.is_finite() {return Err(invalid("supported material initial history overflows"));}
                memories.push(z);rates.push(rate);
            }
            supports.push((*offset,support.clone(),start,local.len()));branches.extend(local);
        }
        for (_,offset,film) in &self.membranes {
            let Some(config)=film.memory_configuration() else {continue;};
            config.admit_domain(joint_omega_bound,self.config.dt_s)?;
            let local=film.memory_branches(*offset,base_dim)?;
            let start=branches.len();
            if start+local.len()>256 || base_dim+start+local.len()>MAX_JOINT_STATES {
                return Err(invalid("complete film and companion material histories exceed the joint state ceiling"));
            }
            for branch in &local {
                let rate=1.0/branch.relaxation_time_s;
                let z=match config.initial {
                    InitialMemory::Relaxed=>det::sqrt(branch.stiffness)*branch.strain.strain(&self.x),
                    InitialMemory::Unrelaxed=>0.0,
                };
                if !rate.is_finite() || rate<=0.0 || !z.is_finite() {
                    return Err(invalid("film material initial history or rate is unrepresentable"));
                }
                memories.push(z);rates.push(rate);
            }
            films.push((*offset,film.clone(),start,local.len()));branches.extend(local);
        }
        if branches.is_empty() && films.is_empty() && supports.iter().all(|(_,s,_,_)|
            s.memory_configuration().is_none_or(|m|!m.has_source_domain())) {return Ok(self);}
        self.system=self.system.with_mapped_relaxation_branches(branches.clone()).map_err(ImpactError::PreparedSolve)?;
        self.x.extend(memories);
        let memory=Memory {base_dim,rows:vec![],rates,mapped:branches,supports,films,joint_omega_bound};
        memory.admit_state(&self.x)?;
        let energy=self.stored_energy_j();
        if !energy.is_finite() || energy<0.0 || energy>self.config.maximum_energy_j {
            return Err(invalid("supported material initial history exceeds the energy ceiling"));
        }
        self.relaxation=Some(memory);Ok(self)
    }

    /// Observe added material memory without advancing mechanics or its history.
    pub fn relaxation_observation(&self) -> Option<RelaxationObservation> {
        self.relaxation.as_ref().map(|memory|memory.observe(&self.x))
    }
}

#[cfg(test)]
#[path = "relaxation_tests.rs"]
mod tests;
