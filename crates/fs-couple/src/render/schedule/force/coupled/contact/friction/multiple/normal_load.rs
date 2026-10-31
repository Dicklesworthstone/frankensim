//! Explicit normal/tangential contact composition on the existing joint owner.
//!
//! This is a ONE-tangent friction law with normal load equal to the actual held
//! compressive reaction in the same discrete step. It is not a vector friction
//! cone, adhesive/thermal history model or a continuous-time impact solution.
//! Contact enrichment and time refinement remain necessary at engagement.

use super::*;

/// Bind one physical tangential contact to its corresponding normal reaction.
/// The producer owns co-location, normals/tangents and material correspondence;
/// modal columns alone cannot establish that two contacts occupy the same point.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ContactNormalLoad {
    /// Original tangential-port index, including port zero.
    pub friction_port: usize,
    /// Existing compliant-obstacle index in its original source order.
    pub obstacle: usize,
}

impl MultiFrictionModalSystem {
    /// Use simultaneously solved normal reactions instead of prescribed loads.
    ///
    /// Supply the complete binding once at sample zero, after obstacles exist.
    /// One obstacle may supply only ONE tangential law: independent scalar cones
    /// on the same normal force would overstate the admissible friction force.
    /// Other ports keep their original prescribed-load semantics. All accepted
    /// motion, preload, interfaces, observers and source laws remain unchanged.
    ///
    /// For each bound port, `FrictionDrive.normal_force_n` MUST be zero: it means
    /// "no additional prescribed load", not "disable friction". A nonzero value
    /// refuses rather than being silently ignored or added. Speed remains the
    /// signed prescribed tangential speed. To release, separate the normal
    /// contact using its existing mechanics/support trajectory. Zero reaction
    /// automatically yields zero traction, regardless of sliding speed.
    ///
    /// Normal and tangential forces are iterated together with both mixed body
    /// responses. The existing normal root, friction root and integration/work
    /// owners are unchanged. This can fail to converge for strong/nonmonotone
    /// coupling; the original bounded sweeps and atomic refusal remain in force.
    /// Nominal pressure/applicability checks use this step's held normal reaction
    /// and initial, average and endpoint tangential slip, not a previous report.
    pub fn with_contact_normal_loads(
        mut self, bindings: Vec<ContactNormalLoad>, gate: &CancelGate,
    ) -> Result<Self, ModalFrictionError> {
        poll(Some(gate))?;
        if bindings.is_empty() { return Ok(self); }
        if self.samples_rendered() != 0 || self.contact_loads.iter().any(Option::is_some)
            || bindings.len() > self.port_count() {
            return Err(invalid("contact-normal loads require one complete binding before playback").into());
        }
        let mut mapping = vec![None; self.port_count()];
        for binding in bindings {
            poll(Some(gate))?;
            if binding.friction_port >= mapping.len() || self.obstacle(binding.obstacle).is_none()
                || mapping[binding.friction_port].is_some()
                || self.mindlin_contacts[binding.friction_port].is_some()
                || mapping.contains(&Some(binding.obstacle)) {
                return Err(invalid("one actual normal contact is required per distinct tangential contact").into());
            }
            mapping[binding.friction_port] = Some(binding.obstacle);
        }
        poll(Some(gate))?;
        self.contact_loads = mapping;
        Ok(self)
    }

    /// Original normal-contact source for this tangential port, or None for a
    /// prescribed-load port/unknown index. This cannot mutate admitted geometry.
    #[must_use]
    pub fn contact_normal_load(&self, friction_port: usize) -> Option<usize> {
        self.contact_loads.get(friction_port).copied().flatten()
    }

    pub(super) fn validate_contact_normal_drives(&self, drives: &[FrictionDrive])
        -> Result<(), ModalFrictionError>
    {
        for (source, drive) in self.contact_loads.iter().zip(drives) {
            if source.is_some() && drive.normal_force_n != 0.0 {
                return Err(invalid("contact-loaded friction cannot also prescribe a normal load").into());
            }
        }
        for (source, drive) in self.mindlin_contacts.iter().zip(drives) {
            if source.is_some() && drive.normal_force_n != 0.0 {
                return Err(invalid("fixed Mindlin preload requires zero additional prescribed normal load; geometry and materials own its retained Hertz load").into());
            }
        }
        Ok(())
    }

    pub(super) fn contact_loaded_drive(&self, port: usize, mut drive: FrictionDrive)
        -> Result<FrictionDrive, ModalFrictionError>
    {
        if let Some(source) = &self.mindlin_contacts[port] {
            drive.normal_force_n = source.normal_force_n();
        } else if let Some(index) = self.contact_loads[port] {
            drive.normal_force_n = self.obstacles.as_ref().and_then(|o| o.normal_reaction(index))
                .ok_or_else(|| invalid("contact-normal binding lost its original obstacle"))?;
            if !drive.normal_force_n.is_finite() || drive.normal_force_n < 0.0 {
                return Err(invalid("resolved normal load must be finite and nonnegative").into());
            }
        }
        Ok(drive)
    }
}

#[cfg(test)]
#[path = "obstacles/normal_load_tests.rs"]
mod tests;
