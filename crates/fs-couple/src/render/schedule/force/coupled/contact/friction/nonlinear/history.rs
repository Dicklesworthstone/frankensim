//! Hereditary constitutive state in the original joint mechanical transaction.
//! Every nonlinear/contact/thermal trial starts from accepted history. Only the
//! existing mechanical publication boundary can commit the proposed history.
use super::*;

/// A passive material with internal coordinates and original mechanical ports.
/// The constitutive owner supplies its discrete update, stored energy and force.
/// The host independently checks Delta(storage)+loss+force.dot(Delta(q))=0.
pub trait ModalMaterialHistory: Send + Sync {
    fn coordinate_count(&self) -> usize;
    fn memory_count(&self) -> usize;
    /// Finite nonnegative upper bound [1/s^2] on the largest eigenvalue of
    /// additional instantaneous mechanical stiffness at fixed material memory,
    /// in these ORIGINAL mass-normalized coordinates. It must cover every
    /// admitted future q/history reachable from this initial memory, including
    /// geometric stress terms. The host retains and aggregates this bound for
    /// its shared material-frequency and Nyquist admission; no zero is inferred.
    fn stiffness_upper_bound_per_s2(&self,initial_memory:&[f64]) -> Result<f64,ModalCouplingError>;
    /// Cold relaxed memory at the actual initial displacement.
    fn relaxed_memory(&self, q: &[f64]) -> Result<Vec<f64>, ModalCouplingError>;
    fn energy(&self, q: &[f64], memory: &[f64]) -> Result<f64, ModalCouplingError>;
    /// Instantaneous force on the mechanical coordinates, for initial guesses
    /// and endpoint acceleration observation. This must not advance history.
    fn reaction(&self, q: &[f64], memory: &[f64], out: &mut [f64]) -> Result<(), ModalCouplingError>;
    /// Propose a complete memory and held work-conjugate force for this physical
    /// interval. Return nonnegative physical dissipation [J]. No accepted state
    /// may be changed by this call, including failed or repeated trials.
    fn stage(&self, before: &[f64], after: &[f64], memory: &[f64], dt: f64,
        candidate: &mut [f64], reaction: &mut [f64]) -> Result<f64, ModalCouplingError>;
}

/// Additional material storage over original coordinates. It may share those
/// coordinates with elastic storage or independent parallel physical mechanisms.
pub struct MappedModalHistoryBinding {
    pub coordinates: Vec<ModalCoordinate>,
    pub law: Box<dyn ModalMaterialHistory>,
    /// None initializes the material relaxed, never at guessed zero strain.
    pub initial_memory: Option<Vec<f64>>,
}

pub(super) struct BoundHistory {
    coordinates: Vec<usize>,
    law: Box<dyn ModalMaterialHistory>,
    accepted: Vec<f64>,
    candidate: Vec<f64>,
    stiffness_bound_per_s2: f64,
}
impl BoundHistory {
    fn select(&self, q: &[f64]) -> Vec<f64> {
        self.coordinates.iter().map(|&i|q[i]).collect()
    }
}

impl FrictionModalSystem {
    /// Add immutable constitutive laws and staged memory to the SAME nonlinear
    /// owner before playback. No new propagator or split history update exists.
    /// Static equilibrium must be initialized before this binding; otherwise a
    /// later static solve could silently reset genuine stored material stress.
    pub fn with_mapped_material_histories(mut self, bindings: Vec<MappedModalHistoryBinding>,
        maximum_states: usize, gate: &CancelGate) -> Result<Self, ModalFrictionError>
    {
        poll(Some(gate))?;
        if self.samples_rendered()!=0 || bindings.is_empty() || self.nonlinear.is_none()
            || !(1..=256).contains(&maximum_states) {
            return Err(invalid("material histories need a cold nonlinear owner and 1..=256 memory allowance").into());
        }
        let q=coordinates(self.components());
        let mut added=Vec::with_capacity(bindings.len());
        let mut count=self.nonlinear.as_ref().expect("checked owner").histories.iter()
            .map(|binding|binding.accepted.len()).sum::<usize>();
        for binding in bindings {
            poll(Some(gate))?;
            let n=binding.law.coordinate_count();let memory=binding.law.memory_count();
            count=count.checked_add(memory).ok_or_else(||invalid("material history count overflow"))?;
            if n==0 || binding.coordinates.len()!=n || memory==0 || count>maximum_states {
                return Err(invalid("complete material coordinate/history budget exceeded").into());
            }
            let mut indices=Vec::with_capacity(n);
            for coordinate in binding.coordinates {
                if self.components().get(coordinate.component)
                    .is_none_or(|c|coordinate.mode>=c.modes().len()) {
                    return Err(invalid("material history names an absent original coordinate").into());
                }
                let index=self.network.offsets[coordinate.component]+coordinate.mode;
                if indices.contains(&index) {return Err(invalid("material history repeats a coordinate").into());}
                indices.push(index);
            }
            let selected:Vec<_>=indices.iter().map(|&i|q[i]).collect();
            let accepted=match binding.initial_memory {
                Some(memory)=>memory,
                None=>binding.law.relaxed_memory(&selected)?,
            };
            if accepted.len()!=memory || accepted.iter().any(|v|!v.is_finite()) {
                return Err(invalid("material history needs complete finite initial memory").into());
            }
            let energy=finite(binding.law.energy(&selected,&accepted)?)?;
            if energy<0.0 {return Err(invalid("material history storage is negative").into());}
            let stiffness_bound_per_s2=finite(binding.law.stiffness_upper_bound_per_s2(&accepted)?)?;
            if stiffness_bound_per_s2<0.0 {
                return Err(invalid("material history requires a nonnegative admitted stiffness bound").into());
            }
            added.push(BoundHistory {coordinates:indices,law:binding.law,
                candidate:accepted.clone(),accepted,stiffness_bound_per_s2});
        }
        self.nonlinear.as_mut().expect("checked owner").histories.extend(added);
        self.validate_material_history_band()?;
        limit("complete initial hereditary material energy",self.total_energy_j()?,
            self.network.config.maximum_total_energy_j)?;
        poll(Some(gate))?;
        Ok(self)
    }

    /// Accepted hereditary histories in binding order; these are separate from
    /// the original linear Maxwell coordinates and retain their declared units.
    pub fn nonlinear_material_memory(&self) -> Vec<&[f64]> {
        self.nonlinear.as_ref().map_or_else(Vec::new,|s|s.histories.iter()
            .map(|binding|binding.accepted.as_slice()).collect())
    }

    pub(in crate::render::schedule::force::coupled::contact::friction) fn validate_material_history_band(&self)
        -> Result<(),ModalCouplingError> {
        if self.nonlinear.as_ref().is_none_or(|state|state.histories.is_empty()) {return Ok(());}
        let high=finite(self.linear_frequency_upper_bound_rad_s())?;
        let guard=core::f64::consts::PI/self.sample_period_s()*self.network.config.nyquist_guard_fraction;
        if high>guard {return Err(invalid("complete hereditary stiffness bound exceeds the original Nyquist guard"));}
        Ok(())
    }
}

impl NonlinearFrictionState {
    pub(in crate::render::schedule::force::coupled::contact::friction) fn history_stiffness_bound_per_s2(&self)->f64 {
        self.histories.iter().map(|binding|binding.stiffness_bound_per_s2).sum()
    }
    pub(in crate::render::schedule::force::coupled::contact::friction) fn publish_histories(&mut self) {
        for binding in &mut self.histories {
            std::mem::swap(&mut binding.accepted,&mut binding.candidate);
        }
    }

    pub(super) fn history_energy(&self,q:&[f64]) -> Result<f64,ModalCouplingError> {
        let mut energy=0.0;
        for binding in &self.histories {
            let value=finite(binding.law.energy(&binding.select(q),&binding.accepted)?)?;
            if value<0.0 {return Err(invalid("material history storage is negative"));}
            energy=finite(energy+value)?;
        }
        Ok(energy)
    }

    pub(super) fn history_reactions(&self,q:&[f64],candidate:bool,out:&mut[f64],gate:Option<&CancelGate>)
        -> Result<(),ModalFrictionError>
    {
        out.fill(0.0);
        for binding in &self.histories {
            poll(gate)?;
            let memory=if candidate {&binding.candidate} else {&binding.accepted};
            let mut force=vec![0.0;binding.coordinates.len()];
            binding.law.reaction(&binding.select(q),memory,&mut force)?;
            for (&i,value) in binding.coordinates.iter().zip(force) {out[i]=finite(out[i]+value)?;}
        }
        Ok(())
    }

    pub(super) fn stage_histories(&mut self,before:&[f64],after:&[f64],dt:f64,out:&mut[f64],
        absolute:f64,relative:f64,gate:Option<&CancelGate>) -> Result<(f64,f64),ModalFrictionError>
    {
        out.fill(0.0);let (mut energy,mut loss)=(0.0,0.0);
        for binding in &mut self.histories {
            poll(gate)?;
            let q0=binding.select(before);let q1=binding.select(after);
            let mut force=vec![0.0;q0.len()];
            let dissipated=finite(binding.law.stage(&q0,&q1,&binding.accepted,dt,
                &mut binding.candidate,&mut force)?)?;
            let initial=finite(binding.law.energy(&q0,&binding.accepted)?)?;
            let stored=finite(binding.law.energy(&q1,&binding.candidate)?)?;
            let work=finite(force.iter().zip(q0.iter().zip(&q1)).map(|(f,(a,b))|f*(b-a)).sum())?;
            let residual=finite(stored-initial+dissipated+work)?;
            let tolerance=finite(absolute+relative*initial.abs().max(stored.abs()).max(dissipated).max(work.abs()))?;
            if initial<0.0 || stored<0.0 || dissipated<0.0 || binding.candidate.iter().any(|v|!v.is_finite()) {
                return Err(invalid("hereditary material must retain finite storage and positive dissipation").into());
            }
            if residual.abs()>tolerance {
                return Err(ModalCouplingError::EnergyBalance {residual_j:residual,tolerance_j:tolerance}.into());
            }
            for (&i,value) in binding.coordinates.iter().zip(force) {out[i]=finite(out[i]+value)?;}
            energy=finite(energy+stored)?;loss=finite(loss+dissipated)?;
        }
        Ok((energy,loss))
    }
}
