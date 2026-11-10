//! The native rigid-body step as a composition-kernel Component.
//!
//! [`RigidBody`] wraps [`RigidBodyIntegrator::step`] behind
//! `fs_ports::component::Component` (contract v0), so a free rigid body can be
//! joined to springs, dampers, joints and heat routes by the kernel in
//! `fs-system` without code written for the combination.
//!
//! Ports, both receiving effort:
//!
//! - `force`: world-frame force on the centre of mass [N], three components.
//!   It returns the mean centre-of-mass velocity of the step, `(x1 - x0) / dt`
//!   [m/s], which is the displacement the integrator pairs with a held force.
//! - `torque`: body-frame torque [N m], three components. It returns the mean
//!   body angular velocity `(omega0 + omega1) / 2` [rad/s], which is the
//!   angular increment the midpoint step pairs with a held torque.
//!
//! With these returns the work the kernel books from the wire values equals
//! the integrator's own applied work to roundoff. Stored energy is kinetic
//! plus uniform-gravity potential (zero at the world origin). The receipt
//! books the integrator's work-energy residual as numerical loss: it is
//! splitting and rounding error, not friction or heat.
//!
//! No claim: no contact, constraint or collision is handled here; the port
//! map has no analytic tangent, so a kernel loop differences it; the split
//! (kick-drift, Verlet) variants of the integrator are not wrapped.
//!
//! [`network::ViscousNetwork`] is the second owner here: several bodies and
//! their viscous connectors as one Component, with a heat port per connector.

/// A damped rigid-body network behind the Component contract.
pub mod network;

use core::num::NonZeroUsize;

use fs_exec::Cx;
use fs_ports::component::{
    Class, Component, ContractViolation, PortRole, PortSpec, Proposal, Receipt,
};
use fs_ports::{PortKind, StableId};

use crate::{DynamicsError, MassProperties, RigidBodyIntegrator, RigidBodyState, Vec3, Wrench};

/// Why a [`RigidBody`] refused a call.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RigidBodyRefusal {
    /// The native integrator refused; its error is unchanged.
    Dynamics(DynamicsError),
    /// Cancellation was requested before the step began.
    Cancelled,
}

impl core::fmt::Display for RigidBodyRefusal {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Dynamics(error) => write!(f, "rigid-body step refused: {error}"),
            Self::Cancelled => write!(f, "cancelled before the rigid-body step"),
        }
    }
}

impl std::error::Error for RigidBodyRefusal {}

impl From<DynamicsError> for RigidBodyRefusal {
    fn from(error: DynamicsError) -> Self {
        Self::Dynamics(error)
    }
}

/// The fixed port list could not be built. Not reachable for the identifiers
/// and roles this module uses; kept as a value so that nothing here panics.
#[derive(Debug, Clone, PartialEq)]
pub enum PortBuildError {
    /// A port identifier was refused.
    Identifier(fs_ports::CoupleError),
    /// A port specification was refused.
    Port(ContractViolation),
}

impl core::fmt::Display for PortBuildError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Identifier(error) => write!(f, "port identifier refused: {error:?}"),
            Self::Port(violation) => write!(f, "port refused: {violation}"),
        }
    }
}

impl std::error::Error for PortBuildError {}

/// One free rigid body, stepped by the energy-momentum midpoint rule.
#[derive(Debug, Clone)]
pub struct RigidBody {
    integrator: RigidBodyIntegrator,
    properties: MassProperties,
    ports: [PortSpec; 2],
}

impl RigidBody {
    /// A body with the given mass properties under the integrator's gravity.
    ///
    /// # Errors
    /// [`PortBuildError`], which the fixed identifiers and roles used here do
    /// not produce.
    pub fn new(
        integrator: RigidBodyIntegrator,
        properties: MassProperties,
    ) -> Result<Self, PortBuildError> {
        let three = NonZeroUsize::MIN.saturating_add(2);
        let port = |name: &str, kind: PortKind| -> Result<PortSpec, PortBuildError> {
            let id = StableId::new(name).map_err(PortBuildError::Identifier)?;
            PortSpec::power(id, kind, PortRole::EffortIn, three).map_err(PortBuildError::Port)
        };
        Ok(Self {
            integrator,
            properties,
            ports: [
                port("force", PortKind::MechanicalForceVelocity)?,
                port("torque", PortKind::RotationalTorqueAngularVelocity)?,
            ],
        })
    }

    /// Mass properties of the body.
    #[must_use]
    pub const fn properties(&self) -> MassProperties {
        self.properties
    }
}

fn components(vector: Vec3) -> [f64; 3] {
    [vector.x, vector.y, vector.z]
}

impl Component for RigidBody {
    type State = RigidBodyState;
    type Error = RigidBodyRefusal;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn stored_energy_j(&self, state: &RigidBodyState) -> Result<f64, RigidBodyRefusal> {
        Ok(self
            .integrator
            .diagnostics(*state, self.properties)?
            .mechanical_energy)
    }

    fn propose(
        &self,
        old: &RigidBodyState,
        dt_s: f64,
        inputs: &[f64],
        new: &mut RigidBodyState,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, RigidBodyRefusal> {
        cx.checkpoint().map_err(|_| RigidBodyRefusal::Cancelled)?;
        let wrench = Wrench {
            force_world: Vec3::new(inputs[0], inputs[1], inputs[2]),
            torque_body: Vec3::new(inputs[3], inputs[4], inputs[5]),
        };
        let receipt = self.integrator.step(*old, self.properties, wrench, dt_s)?;
        let after = receipt.state_after;
        *new = after;

        let start = components(old.pose().position_world());
        let end = components(after.pose().position_world());
        let omega_before = components(
            self.properties
                .angular_velocity_body(old.angular_momentum_body()),
        );
        let omega_after = components(
            self.properties
                .angular_velocity_body(after.angular_momentum_body()),
        );
        let (mut force_sq, mut velocity_sq, mut torque_sq, mut omega_sq) = (0.0_f64, 0.0, 0.0, 0.0);
        for axis in 0..3 {
            let velocity = (end[axis] - start[axis]) / dt_s;
            let omega = 0.5 * (omega_before[axis] + omega_after[axis]);
            outputs[axis] = velocity;
            outputs[3 + axis] = omega;
            force_sq += inputs[axis] * inputs[axis];
            velocity_sq += velocity * velocity;
            torque_sq += inputs[3 + axis] * inputs[3 + axis];
            omega_sq += omega * omega;
        }

        // What the owner's tolerance has to cover, term by term:
        // - the kernel's f64 port products. The angular velocity comes out of
        //   a frame rotation, so each component carries roundoff of the whole
        //   vector's size: products of norms, not of components.
        // - the roundoff of `stored_energy_j` at both ends. Kinetic energy
        //   and gravity potential cancel in free fall, so the measured total
        //   is only as exact as its terms, not as itself.
        let work = receipt.work;
        let (before, after_energy) = (receipt.diagnostics_before, receipt.diagnostics_after);
        let paired_j = dt_s * ((force_sq * velocity_sq).sqrt() + (torque_sq * omega_sq).sqrt());
        let stored_terms_j = before.translational_kinetic_energy
            + before.rotational_kinetic_energy
            + before.gravitational_potential_energy.abs()
            + after_energy.translational_kinetic_energy
            + after_energy.rotational_kinetic_energy
            + after_energy.gravitational_potential_energy.abs();
        let scale_j = paired_j
            + stored_terms_j
            + work.kinetic_energy_change_j.abs()
            + work.gravitational_potential_energy_change_j.abs()
            + work.work_energy_residual_j.abs();
        Ok(Proposal::Step(Receipt {
            stored_change_j: work.mechanical_energy_change_j,
            dissipated_j: 0.0,
            // The residual is `delta(K + V) - applied work`; a positive value
            // is energy the integrator created.
            numerical_loss_j: -work.work_energy_residual_j,
            entropy_generated_j_per_k: 0.0,
            energy_tolerance_j: 64.0 * f64::EPSILON * scale_j + f64::MIN_POSITIVE,
            residual: 0.0,
            tolerance: 0.0,
        }))
    }
}
