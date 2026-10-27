//! Instantaneous configuration updates of the same closed mass inventory.
use super::*;

/// Numerical and physical-work limits for one configuration update.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConfigurationTransferControls {
    /// Relative spin-momentum round-trip residual, in (0, 1e-8]. This does not
    /// bound the intervening dynamics integrator's momentum error.
    pub relative_momentum_tolerance: f64,
    /// Maximum absolute signed change of kinetic energy [J]. An expanding
    /// rotor generally returns negative work. This is parameter work, not heat
    /// or numerical damping, and must be retained by the composing owner.
    pub maximum_parameter_work_j: f64,
}

impl ConfigurationTransferControls {
    pub(crate) fn validate(self) -> Result<(), PrincipalFrameError> {
        if !self.relative_momentum_tolerance.is_finite()
            || self.relative_momentum_tolerance <= 0.0
            || self.relative_momentum_tolerance > 1e-8
            || !self.maximum_parameter_work_j.is_finite()
            || self.maximum_parameter_work_j < 0.0
        {
            return Err(PrincipalFrameError::Invalid(
                "finite transfer work bound and momentum tolerance in (0,1e-8] required",
            ));
        }
        Ok(())
    }
}

/// One immutable proposed state with explicit configuration work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConfigurationTransfer {
    pub state_before: RigidBodyState,
    pub state_after: RigidBodyState,
    pub kinetic_before_j: f64,
    pub kinetic_after_j: f64,
    /// Change due to inertia/configuration alone at conserved world momenta.
    /// No reciprocal thermal energy transfer is inferred here.
    pub parameter_work_j: f64,
    /// Actual round-trip spin-momentum error divided by its original norm.
    /// Zero original spin has exactly zero residual.
    pub relative_momentum_residual: f64,
}

impl PrincipalFrame {
    /// Continue a free body's state through a changed centroid and full inertia.
    ///
    /// Keep the world centre-of-mass position, world linear momentum, world
    /// spin angular momentum, and authored body orientation. The geometry-frame
    /// origin moves as necessary to keep its new centroid at the unchanged
    /// world centre of mass. This also preserves angular momentum about any
    /// fixed world origin, including its orbital contribution. Eigenvector
    /// reorderings and signs therefore cannot rotate the physical geometry or
    /// reset the body's angular velocity.
    ///
    /// The mass must be bit-identical. The caller must establish unchanged
    /// material inventory and body-frame meaning; equal mass alone does not
    /// prove that an arbitrary replacement is this body's configuration.
    /// This quasi-static-shape update omits internal deformation velocity and
    /// its kinetic energy. It is not a contact impulse, an imposed world-space
    /// support, or a finite-strain constitutive/thermal-work closure.
    pub fn transfer_state_from(
        &self,
        previous: &Self,
        state: RigidBodyState,
        controls: ConfigurationTransferControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ConfigurationTransfer, PrincipalFrameError> {
        poll(&mut cancelled)?;
        controls.validate()?;
        if self.properties.mass().to_bits() != previous.properties.mass().to_bits() {
            return Err(PrincipalFrameError::Invalid(
                "configuration transfer requires the identical closed mass inventory",
            ));
        }
        let kinetic_before_j = state.kinetic_energy(previous.properties)?;
        // No source change must retain exact state bits, including at repeated
        // eigenvalues; repeated frame conversions must not accumulate drift.
        if self == previous {
            poll(&mut cancelled)?;
            return Ok(ConfigurationTransfer {
                state_before: state,
                state_after: state,
                kinetic_before_j,
                kinetic_after_j: kinetic_before_j,
                parameter_work_j: 0.0,
                relative_momentum_residual: 0.0,
            });
        }
        let world_from_body = previous.world_from_body(state)?;
        let world_from_principal = world_from_body
            .rotation()
            .compose(self.body_from_principal.rotation())?;
        let q = *world_from_principal.as_quat();
        let orientation = UnitQuaternion::new(q.w, q.x, q.y, q.z)?;
        let spin_world = state
            .pose()
            .orientation()
            .rotate_body_to_world(state.angular_momentum_body());
        let momentum = orientation.rotate_world_to_body(spin_world);
        let state_after = RigidBodyState::new(
            Pose::new(state.pose().position_world(), orientation)?,
            state.linear_momentum_world(),
            momentum,
        )?;
        let actual = orientation.rotate_body_to_world(momentum);
        let residual = actual.sub(spin_world);
        let norm = spin_world.x.hypot(spin_world.y).hypot(spin_world.z);
        let relative_momentum_residual = if norm == 0.0 {
            0.0
        } else {
            residual.x.hypot(residual.y).hypot(residual.z) / norm
        };
        if !relative_momentum_residual.is_finite()
            || relative_momentum_residual > controls.relative_momentum_tolerance
        {
            return Err(PrincipalFrameError::Invalid(
                "configuration transfer exceeded the spin-momentum residual bound",
            ));
        }
        let kinetic_after_j = state_after.kinetic_energy(self.properties)?;
        // The translational term is unchanged. Difference the rotational
        // energies before adding a possibly enormous common translation.
        let rotational_before = 0.5
            * state.angular_momentum_body().dot(
                previous
                    .properties
                    .angular_velocity_body_checked(state.angular_momentum_body())?,
            );
        let rotational_after =
            0.5 * momentum.dot(self.properties.angular_velocity_body_checked(momentum)?);
        let parameter_work_j = rotational_after - rotational_before;
        if !parameter_work_j.is_finite()
            || parameter_work_j.abs() > controls.maximum_parameter_work_j
        {
            return Err(PrincipalFrameError::Invalid(
                "configuration transfer exceeded the signed kinetic-work bound",
            ));
        }
        poll(&mut cancelled)?;
        Ok(ConfigurationTransfer {
            state_before: state,
            state_after,
            kinetic_before_j,
            kinetic_after_j,
            parameter_work_j,
            relative_momentum_residual,
        })
    }
}
