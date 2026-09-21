//! Compile the SAME distributed-air Hamiltonian into the existing prepared
//! modal/contact solver. A spring may involve both heads, acoustic inertia and
//! openings. Concatenating unchanged diagonal coordinates makes that complete
//! signed column one attachment, not several pairwise springs (which would
//! introduce the wrong cross terms). This is a direct sum, not an eigensolve,
//! homogenization, shared wire state, or an additional numerical integrator.
//!
//! All original state addresses survive. Mechanical contact columns gain zeros
//! for air coordinates. Pressure is still observed through CavityCoupling, and
//! exterior sources still project only actual solid motion. Every reaction is
//! solved in the same step as contact; no delayed pressure-forcing loop exists.
use super::CavityCoupling;
use super::super::{BodyPotential, ImpactBody, ImpactError, invalid};
use super::super::linear::{LinearImpactConfig, LinearImpactSystem, VolumeConnection};
use fs_dcontact::Obstacle;
use fs_exec::CancelGate;

impl CavityCoupling {
    /// Prepare linear structural bodies, distributed air and nonlinear contact.
    /// Cavity springs replace, never supplement, the old compact gas spring.
    /// `reference_area_m2` only scales the spring coordinate; it is not a new
    /// physical aperture. The aggregate diagonal component uses `config.component`
    /// as its state/energy budget. The original modal, port and contact limits
    /// are also enforced. Initial body states and contact loss are unchanged.
    ///
    /// This image deliberately refuses nonlinear body potentials, felt pads
    /// (not accepted by this signature), and nonzero acoustic/neck momentum drag.
    /// The existing free-coordinate ZOH owner has no drag law; none is dropped
    /// or replaced by an arbitrary small resonance. Use the nonlinear reference
    /// image for those cases. A finite-step port solve needs time refinement;
    /// it is not the exact full coupled exponential or a real-time certificate.
    pub fn build_linear(self, bodies: Vec<ImpactBody>, contacts: Vec<Obstacle>,
        reference_area_m2: f64, config: LinearImpactConfig, gate: &CancelGate)
        -> Result<(LinearImpactSystem, Self), ImpactError>
    {
        if gate.is_requested() { return Err(ImpactError::Cancelled); }
        if config.sample_rate_hz == 0 || self.total > config.coupling.max_modes
            || !reference_area_m2.is_finite() || reference_area_m2 <= 0.0 {
            return Err(invalid("prepared cavity needs a positive clock/area and sufficient modal budget"));
        }
        if self.damping.iter().any(|d| *d != 0.0)
            || self.necks.iter().any(|n| n.drag_per_s != 0.0) {
            return Err(invalid("prepared cavity cannot discard acoustic or neck momentum drag"));
        }
        let dt_s = f64::from(config.sample_rate_hz).recip();
        let (parts, contacts, _) = self.extend_parts(bodies, contacts, Vec::new(), dt_s, gate)?;
        let mut frequencies = Vec::with_capacity(self.total);
        let mut initial = Vec::with_capacity(self.total);
        let mut damping = Vec::with_capacity(self.total);
        for body in parts {
            if gate.is_requested() { return Err(ImpactError::Cancelled); }
            let BodyPotential::Linear(omega) = body.potential else {
                return Err(invalid("prepared cavity requires linear bodies; nonlinear storage is not discarded"));
            };
            if omega.is_empty() || body.initial.len() != omega.len()
                || body.damping_per_s.len() != omega.len() {
                return Err(invalid("prepared cavity body has an inconsistent diagonal state layout"));
            }
            frequencies.extend(omega);
            initial.extend(body.initial);
            damping.extend(body.damping_per_s);
        }
        let body = ImpactBody { potential: BodyPotential::Linear(frequencies),
            initial, damping_per_s: damping };
        let volumes = self.springs.iter().cloned().map(|spring|
            VolumeConnection { spring, reference_area_m2 }).collect();
        let system = LinearImpactSystem::new(vec![body], contacts, volumes, config, gate)?;
        if gate.is_requested() { return Err(ImpactError::Cancelled); }
        Ok((system, self))
    }
}

#[cfg(test)]
#[path = "prepared_tests.rs"]
mod tests;
