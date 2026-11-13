//! Ideal stationary supports for the original finite planar body/flow owner.

use super::{PlanarRigidError, PlanarRigidMotion, invalid, rates};
use fs_material::fluid_transport::FluidTransport;

/// Homogeneous stationary constraints in fixed world x/y and COM-z rotation.
///
/// A locked translation coordinate keeps the body's current COM coordinate.
/// A locked rotation keeps its current orientation. Unlocked coordinates keep
/// their original finite momentum equation. These are ideal world-aligned
/// supports, not off-COM hinges, moving supports, or arbitrary Pfaffian rows.
/// No impulse is inferred when changing this selection.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct PlanarRigidConstraints {
    pub translation_locked: [bool; 2],
    pub rotation_locked: bool,
}
impl PlanarRigidConstraints {
    pub const FREE: Self = Self {
        translation_locked: [false; 2],
        rotation_locked: false,
    };
    /// Fixed COM with free principal-z rotation: an ideal shaft through the COM.
    pub const SHAFT: Self = Self {
        translation_locked: [true; 2],
        rotation_locked: false,
    };
    /// No body motion; fluid motion and source resolution remain native.
    pub const FIXED: Self = Self {
        translation_locked: [true; 2],
        rotation_locked: true,
    };

    pub(super) const fn free_axes(self) -> [bool; 3] {
        [
            !self.translation_locked[0],
            !self.translation_locked[1],
            !self.rotation_locked,
        ]
    }
}

/// Reaction exerted by the ideal stationary support ON the body.
///
/// The opposite reaction acts on the external support, whose own dynamics are
/// not advanced here. Reactions are not discarded fluid loads: together with
/// fluid, authored external, spring and damping loads they close each actual
/// body momentum equation. Angular quantities are about the body's COM.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct PlanarRigidSupport {
    pub force_n: [f64; 2],
    pub torque_n_m: f64,
    pub impulse_n_s: [f64; 2],
    pub angular_impulse_n_m_s: f64,
    /// Reaction paired with the actual represented coordinate drift [J].
    /// Exact stationarity is required, so an admitted support does zero work.
    pub work_j: f64,
}

impl<F: FluidTransport + ?Sized> PlanarRigidMotion<'_, '_, '_, F> {
    /// Select stationary coordinates without changing accepted momentum,
    /// geometry, fluid history or clock.
    ///
    /// Every newly or previously locked native rate must be EXACTLY zero.
    /// A moving body must first be stopped by an explicit impulse/actuator
    /// owner; this API never silently projects it or discards its kinetic
    /// energy. The same constraints apply to ordinary and rheological steps.
    pub fn set_constraints(
        &mut self,
        constraints: PlanarRigidConstraints,
    ) -> Result<(), PlanarRigidError> {
        let velocity = rates(self.body, self.properties)?;
        let free = constraints.free_axes();
        if (0..3).any(|i| !free[i] && velocity[i] != 0.0) {
            return Err(invalid(
                "stationary supports require exactly zero native rates on locked coordinates",
            ));
        }
        self.constraints = constraints;
        Ok(())
    }

    #[must_use]
    pub const fn constraints(&self) -> PlanarRigidConstraints {
        self.constraints
    }
}

impl PlanarRigidSupport {
    pub(super) fn resolve(
        free: [bool; 3],
        residual: [f64; 3],
        drift: [f64; 3],
        timestep_s: f64,
    ) -> Result<Self, PlanarRigidError> {
        if (0..3).any(|i| !free[i] && drift[i] != 0.0) {
            return Err(invalid(
                "stationary support acquired nonzero represented drift",
            ));
        }
        let reaction: [f64; 3] = core::array::from_fn(|i| if free[i] { 0.0 } else { residual[i] });
        let impulse = reaction.map(|x| timestep_s * x);
        // Check impulses separately: finite forces alone do not make a finite
        // long-interval momentum transfer to the external support.
        if !reaction.into_iter().chain(impulse).all(f64::is_finite) {
            return Err(invalid(
                "unrepresentable stationary support reaction or impulse",
            ));
        }
        Ok(Self {
            force_n: [reaction[0], reaction[1]],
            torque_n_m: reaction[2],
            impulse_n_s: [impulse[0], impulse[1]],
            angular_impulse_n_m_s: impulse[2],
            work_j: 0.0,
        })
    }
}
