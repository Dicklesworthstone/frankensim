//! Constitutive memory in the original coupled held-force owner.
//! q/v remain mechanical coordinates; each arm appends one first-order z.
use super::*;
use fs_math::det;
use fs_phs::RelaxationBranch;

/// A full source tensor tied to original mechanical and material coordinates.
/// This observes the accepted candidate; it introduces no additional storage,
/// forces or integration. `strain_scale` maps an energy quadrature point to
/// the physical surface whose complete source strain must remain admissible.
#[derive(Clone, Debug)]
pub struct ModalPlaneStressPoint {
    pub modes: core::ops::Range<usize>,
    pub memories: core::ops::Range<usize>,
    pub strain_scale: f64,
    pub image: crate::tensor_material::PlaneStressRelaxation,
}

impl ModalPlaneStressPoint {
    fn validate(&self, displacement: &[f64], memory: &[f64]) -> Result<(),ModalCouplingError> {
        let q=displacement.get(self.modes.clone())
            .ok_or_else(||invalid("complete tensor material displacement layout changed"))?;
        let z=memory.get(self.memories.clone())
            .ok_or_else(||invalid("complete tensor material history layout changed"))?;
        let q:Vec<_>=q.iter().map(|q|self.strain_scale*q).collect();
        let z:Vec<_>=z.iter().map(|z|self.strain_scale*z).collect();
        self.image.strain(&q,&z).map_err(|_|invalid("complete tensor material surface strain exceeds its original source domain"))?;
        // Recompute through the ORIGINAL complete law, including normal history;
        // never substitute zero by definition or admit nonfinite normal traction.
        self.image.normal_stress_pa(&q,&z)
            .map_err(|_|invalid("complete tensor material normal traction is unrepresentable"))?;
        Ok(())
    }
}

impl CoupledModalSystem {
    /// Attach source-domain checks after the complete material channels exist.
    /// The fixed geometry/history bindings survive every contact extension.
    /// Both cold state and each staged physical endpoint are checked before
    /// publication, including force-loaded equilibrium. No free predictor is
    /// mistaken for the actual contact-constrained displacement.
    pub fn with_plane_stress_material_points(mut self, points:Vec<ModalPlaneStressPoint>,
        gate:&CancelGate) -> Result<Self,ModalCouplingError> {
        poll(Some(gate))?;
        if self.samples_rendered()!=0 || !self.plane_stress_material_points.is_empty()
            || points.len()>4096 {
            return Err(invalid("tensor material point binding requires one bounded cold attachment"));
        }
        for point in &points {
            poll(Some(gate))?;
            if point.modes.start>point.modes.end || point.modes.end>self.mode_count()
                || point.modes.len()!=point.image.point().coordinates()
                || point.memories.start>point.memories.end || point.memories.end>self.material_memory.len()
                || point.memories.len()!=point.image.point().memory_coordinates()
                || !point.strain_scale.is_finite() || point.strain_scale<=0.0 {
                return Err(invalid("tensor material point requires complete fixed original coordinate ranges"));
            }
        }
        let q:Vec<_>=self.models.iter().flat_map(|m|m.states())
            .map(|s|s.displacement_m_sqrt_kg).collect();
        for point in &points { poll(Some(gate))?; point.validate(&q,&self.material_memory)?; }
        self.plane_stress_material_points=points;
        Ok(self)
    }

    pub(super) fn validate_staged_material_points(&self,gate:Option<&CancelGate>) -> Result<(),ModalCouplingError> {
        if self.plane_stress_material_points.is_empty() {return Ok(());}
        let q:Vec<_>=self.candidates.iter().flat_map(|m|m.states())
            .map(|s|s.displacement_m_sqrt_kg).collect();
        for point in &self.plane_stress_material_points {
            poll(gate)?; point.validate(&q,&self.candidate_material_memory)?;
        }
        Ok(())
    }

    /// Bind Maxwell arms to the SAME complete held-force propagator, once before
    /// playback. Projections use the original interleaved `[q0,v0,q1,v1,...]`
    /// coordinates; velocity coefficients must be zero. The original modal and
    /// connection stiffnesses describe the equilibrium material. Each arm adds
    /// `H = (sqrt(k) projection.q - z)^2 / 2`, with genuine internal strain
    /// `z' = (sqrt(k) projection.q-z)/tau`, never a fictitious mechanical mass.
    ///
    /// `None` explicitly initializes every arm fully relaxed at the current q.
    /// `Some(z)` supplies energy-normalized material history in branch order.
    /// In particular, zero history on a displaced body is an unrelaxed pluck.
    /// Source validity, geometry and modal projection belong to the producer.
    ///
    /// The conservative instantaneous-stiffness Nyquist bound, original total
    /// energy/setup limits and augmented-state capacity all remain binding.
    /// The prepared backend uses its existing budgeted exponential and separate
    /// material-loss Gramian; neither attachment order discards history.
    pub fn with_material_relaxation(
        self,
        branches: Vec<RelaxationBranch>,
        initial_memory: Option<Vec<f64>>,
        gate: &CancelGate,
    ) -> Result<Self, ModalCouplingError> {
        poll(Some(gate))?;
        if self.samples_rendered() != 0 || !self.relaxation_branches.is_empty() {
            return Err(invalid("material relaxation requires one cold binding"));
        }
        self.with_additional_material_relaxation(branches, initial_memory, gate)
    }

    /// Compose another explicitly selected material spectrum before playback.
    /// Existing branches and their accepted internal strains remain an unchanged
    /// prefix. `initial_memory` addresses ONLY the additional branches; `None`
    /// relaxes those new branches at the existing mechanical displacement, never
    /// the already bound material. All branches join one exact propagator and
    /// share its original total-energy, work, memory and bandwidth limits.
    ///
    /// This is a cold constitutive addition, not a hot material replacement.
    /// Callers retain source identities and must not add the same physical loss
    /// twice. Empty additions are a no-op only before the first accepted sample.
    pub fn with_additional_material_relaxation(
        mut self,
        branches: Vec<RelaxationBranch>,
        initial_memory: Option<Vec<f64>>,
        gate: &CancelGate,
    ) -> Result<Self, ModalCouplingError> {
        poll(Some(gate))?;
        if self.samples_rendered() != 0 {
            return Err(invalid("additional material relaxation requires cold construction"));
        }
        if branches.is_empty() {
            if initial_memory.as_ref().is_some_and(|z| !z.is_empty()) {
                return Err(invalid("material memory must match every relaxation branch"));
            }
            return Ok(self);
        }
        let base = self.mode_count().checked_mul(2)
            .ok_or_else(|| invalid("material base-state extent overflow"))?;
        let count = self.relaxation_branches.len().checked_add(branches.len())
            .ok_or_else(|| invalid("material branch count overflow"))?;
        // Bound assembly independently of the per-island hard exponential
        // capacity. No caller's old setup budget silently grows with memory.
        self.admit_augmented_setup(count, self.radiation_state.len())?;
        let mut upper_squared = self.radiation.as_ref().map_or(
            self.linear_frequency_upper_bound_rad_s, |r| r.mechanical_bound).powi(2);
        for branch in &branches {
            poll(Some(gate))?;
            if branch.projection.len() != base
                || !branch.stiffness.is_finite() || branch.stiffness <= 0.0
                || !branch.relaxation_time_s.is_finite() || branch.relaxation_time_s <= 0.0
                || !(1.0 / branch.relaxation_time_s).is_finite()
                || branch.projection.iter().any(|v| !v.is_finite())
                || branch.projection.iter().skip(1).step_by(2).any(|v| *v != 0.0)
            {
                return Err(invalid("material branches require finite displacement-only projections and positive stiffness/time"));
            }
            let root = det::sqrt(branch.stiffness);
            let mut norm = 0.0;
            for &b in branch.projection.iter().step_by(2) {
                let coupling = finite(root * b)?;
                norm = finite(norm + coupling * coupling)?;
            }
            if norm == 0.0 && branch.projection.iter().any(|v| *v!=0.0) {
                return Err(invalid("material stiffness projection is not representable"));
            }
            upper_squared = finite(upper_squared + norm)?;
        }
        let mechanical_bound = det::sqrt(upper_squared);
        let bound = self.radiation.as_ref().map_or(Ok(mechanical_bound), |r|
            r.frequency_bound(mechanical_bound))?;
        if bound > core::f64::consts::PI / self.dt * self.config.nyquist_guard_fraction {
            return Err(invalid("material instantaneous stiffness exceeds the Nyquist guard"));
        }
        let memory = match initial_memory {
            Some(memory) => {
                if memory.len() != branches.len() || memory.iter().any(|z| !z.is_finite()) {
                    return Err(invalid("material memory must be finite and match every relaxation branch"));
                }
                memory
            }
            None => relaxed_memory(&branches, &self.models)?,
        };
        let material_energy = stored_energy(&branches, &self.models, &memory)?;
        limit("total initial material energy", finite(self.total_energy_j()? + material_energy)?,
            self.config.maximum_total_energy_j)?;
        // Existing history is neither reinitialized nor inferred from q. The
        // frequency bound above already includes its instantaneous stiffness;
        // only the newly admitted increments were added to that bound.
        let mut all_branches = self.relaxation_branches.clone();
        all_branches.extend(branches);
        let mut all_memory = self.material_memory.clone();
        all_memory.extend(memory);
        let exact = if let Some(budget) = self.prepared_budget {
            exact::ExactNetwork::new_prepared_with_relaxation(
                &self.models, &self.columns, &self.connections, self.mode_count(), self.dt,
                budget, &all_branches, self.radiation.as_ref(), gate,
            )?
        } else {
            exact::ExactNetwork::new_with_extensions(
                &self.models, &self.columns, &self.connections, self.mode_count(), &all_branches,
                self.radiation.as_ref(), gate,
            )?
        };
        poll(Some(gate))?;
        self.exact = Some(exact);
        self.relaxation_branches = all_branches;
        self.candidate_material_memory = all_memory.clone();
        self.material_memory = all_memory;
        self.linear_frequency_upper_bound_rad_s = bound;
        if let Some(r) = &mut self.radiation { r.mechanical_bound = mechanical_bound; }
        Ok(self)
    }

    /// Original admitted projections and constitutive parameters, in memory order.
    #[must_use]
    pub fn material_relaxation_branches(&self) -> &[RelaxationBranch] {
        &self.relaxation_branches
    }

    /// Accepted energy-normalized internal strains; rejected trials do not change these.
    #[must_use]
    pub fn material_memory(&self) -> &[f64] { &self.material_memory }

    /// Current Maxwell storage, independent of modal and bridge storage [J].
    pub fn material_energy_j(&self) -> Result<f64, ModalCouplingError> {
        stored_energy(&self.relaxation_branches, &self.models, &self.material_memory)
    }

    // Stage material history and its actual endpoint restoring load together.
    // `forces` then serves the existing acceleration-pressure observation.
    pub(super) fn stage_material_state(&mut self, endpoint: &[f64]) -> Result<f64, ModalCouplingError> {
        let base = 2 * self.mode_count();
        let end = base + self.relaxation_branches.len();
        if endpoint.len() != end + self.radiation_state.len() {
            return Err(invalid("complete material endpoint dimension"));
        }
        self.candidate_material_memory.copy_from_slice(&endpoint[base..end]);
        let mut energy = 0.0;
        for (branch, &z) in self.relaxation_branches.iter().zip(&self.candidate_material_memory) {
            let root = det::sqrt(branch.stiffness);
            let strain = branch.projection.iter().zip(&endpoint[..base])
                .try_fold(0.0, |sum, (b, x)| finite(sum + root * b * x))?;
            let mismatch = finite(strain - z)?;
            energy = finite(energy + 0.5 * mismatch * mismatch)?;
            for (force, &b) in self.forces.iter_mut().zip(branch.projection.iter().step_by(2)) {
                *force = finite(*force - root * b * mismatch)?;
            }
        }
        Ok(energy)
    }
}

pub(super) fn relaxed_memory(branches: &[RelaxationBranch], models: &[ModalAcousticTimeModel])
    -> Result<Vec<f64>, ModalCouplingError>
{
    branches.iter().map(|branch| {
        let root = det::sqrt(branch.stiffness);
        branch.projection.iter().step_by(2).zip(models.iter().flat_map(|m| m.states()))
            .try_fold(0.0, |sum, (b, state)| finite(sum + root * b * state.displacement_m_sqrt_kg))
    }).collect()
}

fn stored_energy(branches: &[RelaxationBranch], models: &[ModalAcousticTimeModel], memory: &[f64])
    -> Result<f64, ModalCouplingError>
{
    let relaxed = relaxed_memory(branches, models)?;
    relaxed.iter().zip(memory).try_fold(0.0, |sum, (strain, z)| {
        let mismatch = finite(strain - z)?;
        finite(sum + 0.5 * mismatch * mismatch)
    })
}

#[cfg(test)]
mod tests;
