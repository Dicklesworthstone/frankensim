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
//! - `temperature<e>`, one per connector and only for a network built by
//!   [`ViscousNetwork::temperature_dependent`], a signal input: the
//!   temperature [K] at which that connector's [`DampingLaw`] gives its
//!   passive form for the step, from the form the network was built with.
//!   The form is held over the step. The signal is separate from the heat
//!   port's temperature: where a connector's heat goes and what sets its
//!   damping are authored independently. A kernel solves the signal with
//!   the rest of its block: it needs a first guess for it inside the law's
//!   range, it hands the law the temperature to its solve tolerance and not
//!   to the bit, and one that differences the port map also asks the law
//!   for a form about 1.5e-8 of the temperature above it, so the law has to
//!   be defined that far past any temperature the system reaches.
//! - `torque<i>`, only for a body given one by [`ViscousNetwork::torqued`],
//!   receiving effort: the torque held on that body in its own frame [N m],
//!   three components. It returns the mean of the body's angular velocities
//!   at the two ends of its drift [rad/s], which is the rate the native
//!   step pairs with a held body torque. These ports follow all others, in
//!   the order they were asked for.
//!
//! Receipt: the stored change is kinetic plus uniform-gravity potential
//! energy as the owner reduces them; the dissipation is the owner's total
//! connector heat; the owner's work-heat residual is booked as numerical
//! loss; the entropy generated is `sum(D_e / T_e)`.
//!
//! No claim: a torque is held in the body's frame over the whole step (the
//! owner's world-held couples and point loads are not wrapped); a torque
//! port repeats that body's drift to find the rate, because the owner
//! reports couple work as one total and keeps no state at the end of the
//! drift; the port map has no analytic tangent, so a kernel
//! loop differences it; one propose allocates the load and state banks; the
//! owner's accuracy statements (second order for smooth motion, nothing at
//! unresolved moving geometry) are unchanged by the wrapping.

use core::num::NonZeroUsize;
use std::rc::Rc;

use fs_exec::Cx;
use fs_ports::component::{
    Class, Component, ContractViolation, PortRole, PortSpec, Proposal, Receipt,
};
use fs_ports::{PortKind, StableId};

use crate::spatial_loads::internal_impulse::viscous::PassiveBushingMatrix;
use crate::spatial_loads::internal_impulse::viscous::network::ViscousNetworkLink;
use crate::spatial_loads::internal_motion::HeldExternalLoads;
use crate::spatial_loads::internal_motion::viscous::network::{
    ViscousNetworkMotionControls, ViscousNetworkMotionError, step_viscous_network_with_loads,
};
use crate::{DynamicsError, Gravity, MassProperties, RigidBodyIntegrator, RigidBodyState};

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
    /// A torque port was asked for on a body that does not exist, or on one
    /// that already has it.
    Torque {
        /// Index of the body in the supplied order.
        body: usize,
    },
    /// A port identifier was refused.
    Identifier(fs_ports::CoupleError),
    /// A port specification was refused.
    Port(ContractViolation),
}

impl From<super::PortBuildError> for ViscousNetworkBuildError {
    fn from(error: super::PortBuildError) -> Self {
        match error {
            super::PortBuildError::Identifier(error) => Self::Identifier(error),
            super::PortBuildError::Port(violation) => Self::Port(violation),
        }
    }
}

impl core::fmt::Display for ViscousNetworkBuildError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NoBodies => write!(f, "a viscous network needs at least one body"),
            Self::Link { connector } => write!(
                f,
                "connector {connector} must join two different bodies of the network"
            ),
            Self::Torque { body } => {
                write!(f, "body {body} does not exist or already has a torque port")
            }
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
    /// A temperature signal that is not positive and finite.
    DampingTemperature {
        /// Index of the connector in the supplied order.
        connector: usize,
        /// The value received [K].
        kelvin: f64,
    },
    /// The damping law has no form for a connector at its temperature.
    Law {
        /// Index of the connector in the supplied order.
        connector: usize,
        /// The law's own reason.
        reason: String,
    },
    /// Cancellation was requested.
    Cancelled,
}

/// How a connector's passive form depends on temperature: a memoryless
/// viscous law. Rate- or history-dependent materials, elastic storage and
/// thermal expansion need their own owners.
///
/// The law must be a function: the same arguments, the same form, with no
/// history of its own, because a kernel proposes one step several times.
/// The network declares itself `Smooth`, which holds only as far as the law
/// is differentiable in temperature.
pub trait DampingLaw: core::fmt::Debug {
    /// The passive form of `connector` at `temperature_k`, from `reference`,
    /// the form the network was built with (never an already shifted one).
    ///
    /// # Errors
    /// The law's reason, when it has no form at that temperature.
    fn at_temperature(
        &self,
        connector: usize,
        reference: &PassiveBushingMatrix,
        temperature_k: f64,
    ) -> Result<PassiveBushingMatrix, String>;
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
            Self::DampingTemperature { connector, kelvin } => write!(
                f,
                "connector {connector}'s damping is asked for at {kelvin} K, which is not a \
                 temperature"
            ),
            Self::Law { connector, reason } => {
                write!(f, "connector {connector}'s damping law refused: {reason}")
            }
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
    law: Option<Rc<dyn DampingLaw>>,
    /// The body torque ports, in the order of their ports.
    torques: Vec<super::TorquePort>,
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
            law: None,
            torques: Vec::new(),
            ports,
        })
    }

    /// Gives `body` a torque port `torque<body>`, appended after the ports
    /// the network already has.
    ///
    /// # Errors
    /// [`ViscousNetworkBuildError::Torque`] for a body that does not exist
    /// or already has the port.
    pub fn torqued(mut self, body: usize) -> Result<Self, ViscousNetworkBuildError> {
        if body >= self.properties.len() || self.torques.iter().any(|port| port.body == body) {
            return Err(ViscousNetworkBuildError::Torque { body });
        }
        let port = super::TorquePort::declare(body, &mut self.ports)?;
        self.torques.push(port);
        Ok(self)
    }

    /// The same network with damping that follows temperature: each
    /// connector's form for a step is `law` at the temperature its
    /// `temperature<e>` signal input receives, from the form in `links`.
    ///
    /// # Errors
    /// As [`ViscousNetwork::new`].
    pub fn temperature_dependent(
        properties: Vec<MassProperties>,
        links: Vec<ViscousNetworkLink>,
        gravity: Gravity,
        controls: ViscousNetworkMotionControls,
        law: Rc<dyn DampingLaw>,
    ) -> Result<Self, ViscousNetworkBuildError> {
        let mut network = Self::new(properties, links, gravity, controls)?;
        for connector in 0..network.links.len() {
            let id = StableId::new(format!("temperature{connector}"))
                .map_err(ViscousNetworkBuildError::Identifier)?;
            network
                .ports
                .push(PortSpec::signal_in(id, NonZeroUsize::MIN));
        }
        network.law = Some(law);
        Ok(network)
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

    /// The connectors with the law's forms at the temperatures the signals
    /// carry; nothing for a network without a law.
    fn forms(&self, inputs: &[f64]) -> Result<Vec<ViscousNetworkLink>, ViscousNetworkRefusal> {
        let mut shifted = Vec::new();
        if let Some(law) = &self.law {
            let first = 3 * self.properties.len() + self.links.len();
            let signals = &inputs[first..first + self.links.len()];
            for (connector, link) in self.links.iter().enumerate() {
                let kelvin = signals[connector];
                if !(kelvin.is_finite() && kelvin > 0.0) {
                    return Err(ViscousNetworkRefusal::DampingTemperature { connector, kelvin });
                }
                let damping = law
                    .at_temperature(connector, &link.damping, kelvin)
                    .map_err(|reason| ViscousNetworkRefusal::Law { connector, reason })?;
                shifted.push(ViscousNetworkLink { damping, ..*link });
            }
        }
        Ok(shifted)
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
        let connectors = self.links.len();
        let heat_inputs = &inputs[3 * bodies..3 * bodies + connectors];
        for (connector, &kelvin) in heat_inputs.iter().enumerate() {
            if !(kelvin.is_finite() && kelvin > 0.0) {
                return Err(ViscousNetworkRefusal::Temperature { connector, kelvin });
            }
        }
        // The forms of the step: the built ones, or the law's at the
        // temperatures the signals carry.
        let shifted = self.forms(inputs)?;
        let links = if self.law.is_some() {
            &shifted
        } else {
            &self.links
        };
        let wrenches = super::held_wrenches(bodies, &self.torques, inputs);
        let step = step_viscous_network_with_loads(
            old,
            &self.properties,
            links,
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
        let torques_j = super::torque_rates(
            &self.torques,
            &super::Drift {
                gravity: self.gravity,
                properties: &self.properties,
                starts: &step.half_sources[0].states_after,
                wrenches: &wrenches,
                dt_s,
            },
            inputs,
            outputs,
        )?;
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
            + residual_j.abs()
            + torques_j;
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
