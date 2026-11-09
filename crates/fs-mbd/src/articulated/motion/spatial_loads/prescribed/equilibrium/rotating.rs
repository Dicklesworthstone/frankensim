//! Steady shape in a uniformly rotating supported frame.
//!
//! The SAME original zero-joint-acceleration RNEA residual includes centrifugal
//! loading through real body velocities and inertia. No centrifugal body force,
//! rotating potential, mass matrix or second equilibrium solver is assembled.
//! This is relative equilibrium: bodies have real kinetic energy and the support
//! can supply actual power. Neutral or unstable stationary shapes can be roots.

pub mod frame;
pub use frame::{UniformRotationError, UniformRotationFrame};

use super::super::super::MAX_SPATIAL_MOTION_LOADS;
use super::*;
use crate::articulated::BaseState;
use crate::articulated::motion::sum_products;

/// Native shape balance and the original uniformly rotating boundary.
#[derive(Clone, Debug, PartialEq)]
pub struct RotatingArticulatedEquilibrium {
    /// Original physical joint residuals, material storage and support reaction.
    pub equilibrium: ArticulatedEquilibrium,
    /// Anchored at the seed's clock; use elapsed time from that same instant.
    pub frame: UniformRotationFrame,
    /// Actual support wrench dotted with its body twist, in W into the tree.
    /// A co-rotating applied drag may require nonzero drive power at equilibrium.
    /// This is neither relaxation work nor a heat destination.
    pub support_power_w: f64,
}

impl RotatingArticulatedEquilibrium {
    /// Analytic boundary samples for native prescribed-motion continuation.
    /// Midpoint and endpoint samples must use elapsed times from the root's
    /// original clock. Native geometry admission and roundoff still apply;
    /// the original motion owner controls its timestep.
    pub fn base_at_elapsed(&self, elapsed_s: f64) -> Result<BaseState, UniformRotationError> {
        Ok(BaseState::prescribed(
            self.frame.pose_at_elapsed(elapsed_s)?,
            Twist::new(self.frame.angular_velocity_body(), Twist::zero().linear),
            Twist::zero(),
        ))
    }
}

/// Solve a steady shape under a constant supported angular velocity.
///
/// The seed has a fixed support origin, zero linear base velocity, zero spatial
/// base acceleration, and zero joint rates. Its nonzero angular velocity is
/// retained throughout the same original force-balance solve. Body-frame held
/// loads, original internal axial/bushing laws and constant generalized efforts
/// co-rotate. Every world gravity/point-force/couple vector must be invariant
/// about the rotation axis to the frame's machine-roundoff admission.
/// Arbitrary time/pose-dependent actuator callbacks are intentionally outside
/// this symmetry-admitted lane. Zero spin delegates to the stationary owner.
///
/// This prescribes the spin and reports its support reaction/power; it does not
/// solve free-body angular-momentum equilibrium, select a stable branch, or add
/// centrifugal loads a second time. Its output initializes the existing native
/// prescribed motion owner using `base_at_elapsed` for the boundary samples.
pub fn equilibrate_uniform_rotation(
    model: &ArticulatedModel,
    seed: &PrescribedBaseMotionState,
    loads: SpatialMotionLoads<'_>,
    bushings: EquilibriumBushings<'_>,
    config: EquilibriumConfig<'_>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<RotatingArticulatedEquilibrium, EquilibriumError> {
    let frame = admit_frame(seed, loads, &mut cancelled)?;
    let mode = if seed.base.twist_body.angular == Twist::zero().angular {
        FrameMode::Stationary
    } else {
        FrameMode::UniformRotation
    };
    let equilibrium = dispatch(
        model,
        seed,
        loads,
        bushings,
        config,
        None,
        mode,
        &mut cancelled,
    )?;
    let support_power_w = sum_products(
        equilibrium.support_wrench_body.to_array(),
        equilibrium.state.base.twist_body.to_array(),
    )
    .to_f64();
    if !support_power_w.is_finite() {
        return Err(EquilibriumError::Invalid(
            "support power is not representable",
        ));
    }
    poll(&mut cancelled)?;
    Ok(RotatingArticulatedEquilibrium {
        equilibrium,
        frame,
        support_power_w,
    })
}

// Shared symmetry admission for the root and its original local response.
pub(super) fn admit_frame(
    seed: &PrescribedBaseMotionState,
    loads: SpatialMotionLoads<'_>,
    mut cancelled: impl FnMut() -> bool,
) -> Result<UniformRotationFrame, EquilibriumError> {
    poll(&mut cancelled)?;
    let frame = UniformRotationFrame::new(seed.base.world_from_base, seed.base.twist_body.angular)
        .map_err(frame_error)?;
    if loads.world_points.len() > MAX_SPATIAL_MOTION_LOADS {
        return Err(EquilibriumError::Invalid(
            "bounded world point inventory required",
        ));
    }
    if !frame.admits_world_vector(loads.held.gravity_world) {
        return Err(EquilibriumError::Invalid(
            "world gravity is not invariant under the supported rotation",
        ));
    }
    for (i, load) in loads.world_points.iter().enumerate() {
        if i % 64 == 0 {
            poll(&mut cancelled)?;
        }
        if !frame.admits_world_vector(load.force_world)
            || !frame.admits_world_vector(load.couple_world)
        {
            return Err(EquilibriumError::Invalid(
                "world point load is not invariant under the supported rotation",
            ));
        }
    }
    poll(&mut cancelled)?;
    Ok(frame)
}

fn frame_error(error: UniformRotationError) -> EquilibriumError {
    match error {
        UniformRotationError::Invalid(reason) => EquilibriumError::Invalid(reason),
        UniformRotationError::Geometry(error) => PrescribedBaseMotionError::from(error).into(),
    }
}
