//! Compile the SAME distributed-air Hamiltonian into the existing prepared
//! modal/contact solver. A spring may involve both heads, acoustic inertia and
//! openings. Complete signed columns retain every cross term while keeping
//! each original solid body and the appended acoustic inertia as independently
//! budgeted components. No bodies are regrouped, no wire is homogenized, and
//! no additional eigenbasis or numerical integrator is introduced.
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
    /// physical aperture. Each original body keeps its `config.component`
    /// state/energy budget; appended acoustic/neck inertia has its own component
    /// with the same budget. The whole-system modal, port, contact and energy
    /// limits still apply. Initial body states and contact loss are unchanged.
    ///
    /// Acoustic/neck momentum drag is retained through the linear impact
    /// owner's simultaneous grounded viscous links, consuming one connection
    /// per nonzero drag in addition to the cavity springs and any solid ports.
    /// No loss is substituted onto the heads, counted twice, or omitted.
    /// Nonlinear body potentials and felt pads (absent from this signature)
    /// still require the reference image. Finite-step ports need refinement;
    /// it is not the exact full coupled exponential or a real-time certificate.
    pub fn build_linear(self, bodies: Vec<ImpactBody>, contacts: Vec<Obstacle>,
        reference_area_m2: f64, config: LinearImpactConfig, gate: &CancelGate)
        -> Result<(LinearImpactSystem, Self), ImpactError>
    {
        self.build_linear_with_dampers(bodies, contacts, Vec::new(), reference_area_m2, config, gate)
    }

    /// Preserve solid viscous ports while appending acoustic inertia. Original
    /// mode rows gain exact gas zeros before the shared bilateral/contact solve.
    pub fn build_linear_with_dampers(self, bodies: Vec<ImpactBody>, contacts: Vec<Obstacle>,
        dampers: Vec<super::super::damping::ViscousDamper>, reference_area_m2: f64,
        config: LinearImpactConfig, gate: &CancelGate) -> Result<(LinearImpactSystem, Self), ImpactError>
    {
        if gate.is_requested() { return Err(ImpactError::Cancelled); }
        if config.sample_rate_hz == 0 || self.total > config.coupling.max_modes
            || !reference_area_m2.is_finite() || reference_area_m2 <= 0.0 {
            return Err(invalid("prepared cavity needs a positive clock/area and sufficient modal budget"));
        }
        let dt_s = f64::from(config.sample_rate_hz).recip();
        let (parts, contacts, _) = self.extend_parts(bodies, contacts, Vec::new(), dt_s, gate)?;
        for body in &parts {
            if gate.is_requested() { return Err(ImpactError::Cancelled); }
            let BodyPotential::Linear(omega) = &body.potential else {
                return Err(invalid("prepared cavity requires linear bodies; nonlinear storage is not discarded"));
            };
            if omega.is_empty() || body.initial.len() != omega.len()
                || body.damping_per_s.len() != omega.len() {
                return Err(invalid("prepared cavity body has an inconsistent diagonal state layout"));
            }
        }
        let volumes = self.springs.iter().cloned().map(|spring|
            VolumeConnection { spring, reference_area_m2 }).collect();
        let dampers = super::super::damping::extend(dampers, self.structural, self.total)?;
        let system = LinearImpactSystem::new_with_dampers(parts, contacts, volumes, dampers, config, gate)?;
        if gate.is_requested() { return Err(ImpactError::Cancelled); }
        Ok((system, self))
    }
}

#[cfg(test)]
#[path = "prepared_tests.rs"]
mod tests;
