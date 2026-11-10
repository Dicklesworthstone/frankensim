//! A damped rigid-body network as one composition-kernel Component.
//!
//! [`ViscousNetwork`] wraps
//! [`step_viscous_network_with_loads`]: any number of rigid bodies joined by
//! passive viscous connectors, advanced by the owner's own splitting (half a
//! connector source at fixed geometry, the native rigid motion of every body,
//! half a source again). Nothing is integrated here.
//!
//! Ports, in this order:
//!
//! - `force<i>`, one per body, receiving effort: the world-frame force on the
//!   centre of mass [N], three components. It returns the mean
//!   centre-of-mass velocity of the step, `(x1 - x0) / dt` [m/s]. The
//!   connector sources change momenta at fixed geometry, so a body moves only
//!   in its native rigid motion and this is the displacement the owner pairs
//!   with a held force.
//! - `heat<e>`, one per connector, a heat port receiving effort: the
//!   temperature at which that connector's heat is released [K]. It returns
//!   minus the connector's mean dissipation rate, `-D_e / dt` [W]: heat
//!   leaves the network there, to wherever the system routes it.
//!
//! Receipt: the stored change is kinetic plus uniform-gravity potential
//! energy as the owner reduces them; the dissipation is the owner's total
//! connector heat; the owner's work-heat residual is booked as numerical
//! loss; the entropy generated is `sum(D_e / T_e)`.
//!
//! No claim: there is no torque port, because the owner reports couple work
//! as one total and not per body, so a couple cannot be attached to a body
//! through the kernel yet; the port map has no analytic tangent, so a kernel
//! loop differences it; one propose allocates the load and state banks; the
//! owner's accuracy statements (second order for smooth motion, nothing at
//! unresolved moving geometry) are unchanged by the wrapping.

use core::num::NonZeroUsize;

use fs_exec::Cx;
use fs_ports::component::{
    Class, Component, ContractViolation, PortRole, PortSpec, Proposal, Receipt,
};
use fs_ports::{PortKind, StableId};

use crate::spatial_loads::internal_impulse::viscous::network::ViscousNetworkLink;
use crate::spatial_loads::internal_motion::HeldExternalLoads;
use crate::spatial_loads::internal_motion::viscous::network::{
    ViscousNetworkMotionControls, ViscousNetworkMotionError, step_viscous_network_with_loads,
};
use crate::{
    DynamicsError, Gravity, MassProperties, RigidBodyIntegrator, RigidBodyState, Vec3, Wrench,
};

/// Why a [`ViscousNetwork`] could not be built.
#[derive(Debug, Clone, PartialEq)]
pub enum ViscousNetworkBuildError {
    /// A network needs at least one body.
    NoBodies,
    /// A connector names a body that does not exist, or the same body twice.
    Link {
        /// Index of the connector in the supplied order.
        connector: usize,
    },
    /// A port identifier was refused.
    Identifier(fs_ports::CoupleError),
    /// A port specification was refused.
    Port(ContractViolation),
}

impl core::fmt::Display for ViscousNetworkBuildError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NoBodies => write!(f, "a viscous network needs at least one body"),
            Self::Link { connector } => write!(
                f,
                "connector {connector} must join two different bodies of the network"
            ),
            Self::Identifier(error) => write!(f, "port identifier refused: {error:?}"),
            Self::Port(violation) => write!(f, "port refused: {violation}"),
        }
    }
}

impl std::error::Error for ViscousNetworkBuildError {}

/// Why a [`ViscousNetwork`] refused a call.
#[derive(Debug, Clone, PartialEq)]
pub enum ViscousNetworkRefusal {
    /// The owner refused the step; its error is unchanged.
    Motion(ViscousNetworkMotionError),
    /// The energy diagnostics refused a state.
    Dynamics(DynamicsError),
    /// The state does not hold one entry per body.
    Bodies {
        /// Bodies the network was built with.
        expected: usize,
        /// Entries in the state.
        found: usize,
    },
    /// A heat port received a temperature that is not positive and finite.
    Temperature {
        /// Index of the connector in the supplied order.
        connector: usize,
        /// The value received [K].
        kelvin: f64,
    },
    /// Cancellation was requested.
    Cancelled,
}

impl core::fmt::Display for ViscousNetworkRefusal {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Motion(error) => write!(f, "{error}"),
            Self::Dynamics(error) => write!(f, "network energy refused: {error}"),
            Self::Bodies { expected, found } => write!(
                f,
                "the network has {expected} bodies and the state holds {found}"
            ),
            Self::Temperature { connector, kelvin } => write!(
                f,
                "connector {connector} releases heat at {kelvin} K, which is not a temperature"
            ),
            Self::Cancelled => write!(f, "cancelled during the network step"),
        }
    }
}

impl std::error::Error for ViscousNetworkRefusal {}

impl From<ViscousNetworkMotionError> for ViscousNetworkRefusal {
    fn from(error: ViscousNetworkMotionError) -> Self {
        match error {
            ViscousNetworkMotionError::Cancelled => Self::Cancelled,
            error => Self::Motion(error),
        }
    }
}

impl From<DynamicsError> for ViscousNetworkRefusal {
    fn from(error: DynamicsError) -> Self {
        Self::Dynamics(error)
    }
}

/// Rigid bodies joined by passive viscous connectors, stepped by their owner.
#[derive(Debug, Clone)]
pub struct ViscousNetwork {
    properties: Vec<MassProperties>,
    links: Vec<ViscousNetworkLink>,
    gravity: Gravity,
    controls: ViscousNetworkMotionControls,
    ports: Vec<PortSpec>,
}

impl ViscousNetwork {
    /// A network of `properties.len()` bodies and `links.len()` connectors
    /// under uniform `gravity`, stepped with the owner's `controls`.
    ///
    /// # Errors
    /// [`ViscousNetworkBuildError`] for no bodies or a connector that does
    /// not join two different bodies of the network.
    pub fn new(
        properties: Vec<MassProperties>,
        links: Vec<ViscousNetworkLink>,
        gravity: Gravity,
        controls: ViscousNetworkMotionControls,
    ) -> Result<Self, ViscousNetworkBuildError> {
        let bodies = properties.len();
        if bodies == 0 {
            return Err(ViscousNetworkBuildError::NoBodies);
        }
        for (connector, link) in links.iter().enumerate() {
            let [first, second] = link.bodies;
            if first >= bodies || second >= bodies || first == second {
                return Err(ViscousNetworkBuildError::Link { connector });
            }
        }
        let three = NonZeroUsize::MIN.saturating_add(2);
        let mut ports = Vec::with_capacity(bodies + links.len());
        for body in 0..bodies {
            let id = StableId::new(format!("force{body}"))
                .map_err(ViscousNetworkBuildError::Identifier)?;
            ports.push(
                PortSpec::power(
                    id,
                    PortKind::MechanicalForceVelocity,
                    PortRole::EffortIn,
                    three,
                )
                .map_err(ViscousNetworkBuildError::Port)?,
            );
        }
        for connector in 0..links.len() {
            let id = StableId::new(format!("heat{connector}"))
                .map_err(ViscousNetworkBuildError::Identifier)?;
            ports.push(
                PortSpec::heat(id, PortRole::EffortIn).map_err(ViscousNetworkBuildError::Port)?,
            );
        }
        Ok(Self {
            properties,
            links,
            gravity,
            controls,
            ports,
        })
    }

    /// Number of bodies; ports `0..bodies()` are their force ports.
    #[must_use]
    pub fn bodies(&self) -> usize {
        self.properties.len()
    }

    /// Number of connectors; ports `bodies()..bodies() + connectors()` are
    /// their heat ports.
    #[must_use]
    pub fn connectors(&self) -> usize {
        self.links.len()
    }

    /// Kinetic terms and gravity potential of one bank, as sizes for a
    /// roundoff allowance and as the total.
    fn energy(&self, state: &[RigidBodyState]) -> Result<(f64, f64), ViscousNetworkRefusal> {
        if state.len() != self.properties.len() {
            return Err(ViscousNetworkRefusal::Bodies {
                expected: self.properties.len(),
                found: state.len(),
            });
        }
        let integrator = RigidBodyIntegrator::new(self.gravity);
        let (mut total, mut terms) = (0.0_f64, 0.0_f64);
        for (body, properties) in state.iter().zip(&self.properties) {
            let diagnostics = integrator.diagnostics(*body, *properties)?;
            total += diagnostics.mechanical_energy;
            terms += diagnostics.translational_kinetic_energy
                + diagnostics.rotational_kinetic_energy
                + super::potential_terms_j(self.gravity, *properties, *body);
        }
        Ok((total, terms))
    }
}

impl Component for ViscousNetwork {
    type State = Vec<RigidBodyState>;
    type Error = ViscousNetworkRefusal;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn stored_energy_j(&self, state: &Vec<RigidBodyState>) -> Result<f64, ViscousNetworkRefusal> {
        Ok(self.energy(state)?.0)
    }

    fn propose(
        &self,
        old: &Vec<RigidBodyState>,
        dt_s: f64,
        inputs: &[f64],
        new: &mut Vec<RigidBodyState>,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ViscousNetworkRefusal> {
        cx.checkpoint()
            .map_err(|_| ViscousNetworkRefusal::Cancelled)?;
        let bodies = self.properties.len();
        let (_, terms_before_j) = self.energy(old)?;
        let heat_inputs = &inputs[3 * bodies..];
        for (connector, &kelvin) in heat_inputs.iter().enumerate() {
            if !(kelvin.is_finite() && kelvin > 0.0) {
                return Err(ViscousNetworkRefusal::Temperature { connector, kelvin });
            }
        }
        let wrenches: Vec<Wrench> = (0..bodies)
            .map(|body| Wrench {
                force_world: Vec3::new(
                    inputs[3 * body],
                    inputs[3 * body + 1],
                    inputs[3 * body + 2],
                ),
                torque_body: Vec3::ZERO,
            })
            .collect();
        let step = step_viscous_network_with_loads(
            old,
            &self.properties,
            &self.links,
            HeldExternalLoads {
                wrenches: &wrenches,
                gravity: self.gravity,
            },
            dt_s,
            self.controls,
            || cx.is_cancel_requested(),
        )?;
        let after = step.states_after();

        // Work the kernel will book on the force ports, as a size.
        let mut paired_j = 0.0_f64;
        for body in 0..bodies {
            let start = old[body].pose().position_world();
            let end = after[body].pose().position_world();
            let velocity = [
                (end.x - start.x) / dt_s,
                (end.y - start.y) / dt_s,
                (end.z - start.z) / dt_s,
            ];
            let (mut force_sq, mut velocity_sq) = (0.0_f64, 0.0_f64);
            for axis in 0..3 {
                outputs[3 * body + axis] = velocity[axis];
                force_sq += inputs[3 * body + axis] * inputs[3 * body + axis];
                velocity_sq += velocity[axis] * velocity[axis];
            }
            paired_j += dt_s * (force_sq * velocity_sq).sqrt();
        }
        // A connector's heat is a sum of two source integrals and is never
        // negative; a rounding residue below zero is not entropy destroyed.
        let (mut heat_terms_j, mut entropy_j_per_k) = (0.0_f64, 0.0_f64);
        for (connector, &heat_j) in step.connector_dissipation_j.iter().enumerate() {
            outputs[3 * bodies + connector] = -heat_j / dt_s;
            heat_terms_j += heat_j.abs();
            entropy_j_per_k += heat_j.max(0.0) / heat_inputs[connector];
        }
        let stored_change_j = step.mechanical_energy_change_j;
        let residual_j = step.work_heat_residual_j;
        let dissipated_j = step.dissipated_energy_j.max(0.0);

        new.clear();
        new.extend_from_slice(after);
        let (_, terms_after_j) = self.energy(new)?;

        // The tolerance covers the kernel's f64 port products, the roundoff
        // of `stored_energy_j` at both ends (kinetic energy and gravity
        // potential cancel in free fall), and the owner's own rounded sums.
        let scale_j = paired_j
            + heat_terms_j
            + terms_before_j
            + terms_after_j
            + stored_change_j.abs()
            + residual_j.abs();
        Ok(Proposal::Step(Receipt {
            stored_change_j,
            dissipated_j,
            // The residual is `delta(K + V) + heat - applied work`; a
            // positive value is energy the splitting created.
            numerical_loss_j: -residual_j,
            entropy_generated_j_per_k: entropy_j_per_k,
            energy_tolerance_j: 64.0 * f64::EPSILON * scale_j + f64::MIN_POSITIVE,
            residual: 0.0,
            tolerance: 0.0,
        }))
    }
}
