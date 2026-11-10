//! A rigid-body network with Maxwell memory as one composition-kernel
//! Component.
//!
//! [`RelaxationNetwork`] wraps [`step_relaxation_network_with_loads`] (and,
//! for kinetics that follow temperature, its held-kinetics twin): rigid
//! bodies joined by connectors that have, beside their Kelvin viscosity,
//! Maxwell (or lossless equilibrium) branches, each with its own elastic
//! memory. The owner solves bodies and branches as one system and admits
//! what fits its dense solve (at this writing, bodies plus branches up to
//! 21); a larger network is built here and refused by the owner at its
//! first step. The owner advances bodies and every memory coordinate
//! together (an exact simultaneous source at fixed geometry, the native
//! rigid motion, a source again). Nothing is integrated here.
//!
//! State: [`RelaxationState`], the bank of bodies and the complete memory, a
//! six-component energy amplitude `a` [sqrt(J)] per branch with storage
//! `|a|^2 / 2`. Zero is relaxed memory; the current pose does not determine
//! a preloaded history, so the memory is part of the state and never
//! inferred.
//!
//! Ports, as for [`super::network::ViscousNetwork`]: the force ports in body
//! order, then the heat ports in connector order, then the two ports of each
//! actuator, the temperature signals and the body torque ports, in the
//! order they were asked for.
//!
//! - `force<i>`, one per body, receiving effort: the world-frame force on
//!   the centre of mass [N], three components. It returns the mean
//!   centre-of-mass velocity of the step, `(x1 - x0) / dt` [m/s].
//! - `heat<e>`, one per connector, a heat port receiving effort: the
//!   temperature at which that connector's heat is released [K]. It returns
//!   minus the connector's mean dissipation rate, `-D_e / dt` [W], Kelvin
//!   and Maxwell loss together.
//! - for each connector given an actuator by [`RelaxationNetwork::driven`],
//!   two more ports receiving effort: `drive<e>`, the force the actuator
//!   holds on the connector's second attachment [N], and `couple<e>`, the
//!   couple it holds there [N m], three components each, in the axes of the
//!   port's [`DriveFrame`]. The first body takes the opposite wrench,
//!   carried to the second attachment. The ports return the rates conjugate
//!   to those efforts in the same axes: the velocity of the second
//!   attachment against the point of the first body that is there [m/s],
//!   and the relative angular velocity [rad/s], each as the mean of the
//!   owner's two half-step averages (the mean over the step to second
//!   order in its length). Their products with the efforts, times the step,
//!   are the owner's actuator work on that connector, which may be negative
//!   (recovery).
//! - `temperature<e>`, one per connector and only for a network that
//!   follows temperature, a signal input [K], appended when
//!   [`RelaxationNetwork::temperature_dependent`] or
//!   [`RelaxationNetwork::shifted`] is called. With the first, the
//!   connector's Kelvin form for the step is its [`DampingLaw`] at that
//!   temperature, from the form the network was built with, and the Maxwell
//!   branches do not change. With the second, the connector's [`ShiftLaw`]
//!   gives one factor at that temperature, by which its Kelvin viscosity
//!   and every relaxation time of its branches are multiplied; elastic
//!   factors, memory and storage do not change. Either is held over the
//!   step. What is said of this signal for the viscous network holds here.
//! - `torque<i>`, only for a body given one by
//!   [`RelaxationNetwork::torqued`], receiving effort: the torque held on
//!   that body in its own frame [N m]. It returns the mean of the body's
//!   angular velocities at the two ends of its drift [rad/s], the rate the
//!   viscous network's torque port returns.
//!
//! Receipt: the stored change is kinetic plus uniform-gravity potential plus
//! elastic memory, as the owner reduces them; the dissipation is the owner's
//! total connector heat; the owner's work-heat residual is booked as
//! numerical loss; the entropy generated is `sum(D_e / T_e)`. Elastic
//! recovery can speed a body up while its connector still heats: only the
//! independently integrated loss is heat.
//!
//! No claim: an effort is held over the whole step (the owner's source is
//! exact for a held effort at fixed geometry; a varying command is sampled
//! by whatever feeds the port); a connector without drive ports has no
//! actuator; a body torque is held in the body's frame, and its port
//! repeats that body's drift, as for the viscous network; a
//! network follows temperature by a damping law or by a shift factor, not
//! both, and a shift factor is one number per connector (the owner's
//! kinetics take one per Kelvin form and one per branch); the port map has
//! no analytic tangent;
//! cancellation raised while the owner is inside a step is passed to it but
//! not exercised by a test; the owner's
//! accuracy statements (exact sources at fixed geometry, second order for
//! the moving geometry) are unchanged by the wrapping.

use core::num::NonZeroUsize;
use std::rc::Rc;

use fs_exec::Cx;
use fs_ports::component::{
    Class, Component, ContractViolation, PortRole, PortSpec, Proposal, Receipt,
};
use fs_ports::{PortKind, StableId};

use crate::spatial_loads::internal_impulse::viscous::network::relaxation::{
    RelaxationNetworkLink, RelaxationNetworkMemory,
};
use crate::spatial_loads::internal_motion::HeldExternalLoads;
use crate::spatial_loads::internal_motion::viscous::network::relaxation::{
    HeldInternalEffort, HeldRelaxationKinetics, RelaxationNetworkMotionControls,
    RelaxationNetworkMotionError, RelaxationNetworkMotionStep, step_relaxation_network_with_loads,
    step_relaxation_network_with_loads_and_kinetics,
};
use crate::{
    DynamicsError, Gravity, MassProperties, RigidBodyIntegrator, RigidBodyState, Vec3, Wrench,
};

pub use super::network::DampingLaw;

/// How a connector's kinetics follow temperature: one factor per connector,
/// `a_T = tau(T) / tau_ref`, by which its Kelvin viscosity and every
/// relaxation time of its Maxwell branches are multiplied (a
/// thermorheologically simple connector). Elastic factors, the memory and
/// its storage do not change with it, and a lossless equilibrium branch
/// stays lossless.
///
/// A factor is a ratio of times: it must be positive and finite, and a
/// network refuses a step for which the law returns anything else. The law
/// must be a function, as a [`DampingLaw`] must, and the network's `Smooth`
/// class holds only as far as the factor is differentiable in temperature.
pub trait ShiftLaw: core::fmt::Debug {
    /// The factor of `connector` at `temperature_k`.
    ///
    /// # Errors
    /// The law's reason, when it has no factor at that temperature.
    fn factor(&self, connector: usize, temperature_k: f64) -> Result<f64, String>;
}

/// What the `temperature<e>` signals set.
#[derive(Debug, Clone)]
enum Thermal {
    /// Nothing: the network has no such signals.
    Fixed,
    /// Each connector's Kelvin form.
    Damping(Rc<dyn DampingLaw>),
    /// Each connector's Kelvin viscosity and relaxation times.
    Kinetics(Rc<dyn ShiftLaw>),
}

/// One actuator: its connector, its axes, and where its two ports'
/// coefficients start in the input and output vectors.
#[derive(Debug, Clone, Copy)]
struct Drive {
    connector: usize,
    frame: DriveFrame,
    input: usize,
    output: usize,
}

/// The axes in which a drive port's effort is held during a step, and in
/// which its conjugate rate is returned.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DriveFrame {
    /// The material axes of the connector's first body: the held components
    /// turn with that body.
    MaterialFirst,
    /// World axes: the held components stay fixed while the bodies turn.
    World,
}

/// Why a [`RelaxationNetwork`] could not be built.
#[derive(Debug, Clone, PartialEq)]
pub enum RelaxationNetworkBuildError {
    /// A network needs at least one body.
    NoBodies,
    /// A connector names a body that does not exist, or the same body twice.
    Link {
        /// Index of the connector in the supplied order.
        connector: usize,
    },
    /// An actuator was asked for on a connector that does not exist, or on
    /// one that already has it.
    Drive {
        /// Index of the connector in the supplied order.
        connector: usize,
    },
    /// The network already follows temperature.
    Thermal,
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

impl core::fmt::Display for RelaxationNetworkBuildError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::NoBodies => write!(f, "a relaxation network needs at least one body"),
            Self::Link { connector } => write!(
                f,
                "connector {connector} must join two different bodies of the network"
            ),
            Self::Drive { connector } => write!(
                f,
                "connector {connector} does not exist or already has an actuator"
            ),
            Self::Thermal => write!(f, "the network already follows temperature"),
            Self::Torque { body } => {
                write!(f, "body {body} does not exist or already has a torque port")
            }
            Self::Identifier(error) => write!(f, "port identifier refused: {error:?}"),
            Self::Port(violation) => write!(f, "port refused: {violation}"),
        }
    }
}

impl std::error::Error for RelaxationNetworkBuildError {}

impl From<super::PortBuildError> for RelaxationNetworkBuildError {
    fn from(error: super::PortBuildError) -> Self {
        match error {
            super::PortBuildError::Identifier(error) => Self::Identifier(error),
            super::PortBuildError::Port(violation) => Self::Port(violation),
        }
    }
}

/// Why a [`RelaxationNetwork`] refused a call.
#[derive(Debug, Clone, PartialEq)]
pub enum RelaxationNetworkRefusal {
    /// The owner refused the step; its error is unchanged.
    Motion(RelaxationNetworkMotionError),
    /// The energy diagnostics refused a state.
    Dynamics(DynamicsError),
    /// The state does not hold one entry per body.
    Bodies {
        /// Bodies the network was built with.
        expected: usize,
        /// Entries in the state.
        found: usize,
    },
    /// The state's memory does not hold one amplitude per branch of every
    /// connector, or holds a value that is not finite.
    Memory {
        /// Index of the first connector whose memory does not fit; the
        /// number of connectors when the memory has entries for more.
        connector: usize,
    },
    /// A heat port received a temperature that is not positive and finite.
    Temperature {
        /// Index of the connector in the supplied order.
        connector: usize,
        /// The value received [K].
        kelvin: f64,
    },
    /// A drive port received an effort that is not finite.
    Effort {
        /// Index of the connector in the supplied order.
        connector: usize,
    },
    /// A temperature signal that is not positive and finite.
    DampingTemperature {
        /// Index of the connector in the supplied order.
        connector: usize,
        /// The value received [K].
        kelvin: f64,
    },
    /// The damping law or the shift law has nothing for a connector at its
    /// temperature.
    Law {
        /// Index of the connector in the supplied order.
        connector: usize,
        /// The law's own reason.
        reason: String,
    },
    /// Cancellation was requested.
    Cancelled,
}

impl core::fmt::Display for RelaxationNetworkRefusal {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Motion(error) => write!(f, "{error}"),
            Self::Dynamics(error) => write!(f, "network energy refused: {error}"),
            Self::Bodies { expected, found } => write!(
                f,
                "the network has {expected} bodies and the state holds {found}"
            ),
            Self::Memory { connector } => write!(
                f,
                "the memory of connector {connector} is not one finite amplitude per branch"
            ),
            Self::Temperature { connector, kelvin } => write!(
                f,
                "connector {connector} releases heat at {kelvin} K, which is not a temperature"
            ),
            Self::Effort { connector } => write!(
                f,
                "the actuator of connector {connector} received an effort that is not finite"
            ),
            Self::DampingTemperature { connector, kelvin } => write!(
                f,
                "connector {connector}'s law is asked for at {kelvin} K, which is not a \
                 temperature"
            ),
            Self::Law { connector, reason } => {
                write!(f, "connector {connector}'s temperature law refused: {reason}")
            }
            Self::Cancelled => write!(f, "cancelled during the network step"),
        }
    }
}

impl std::error::Error for RelaxationNetworkRefusal {}

impl From<RelaxationNetworkMotionError> for RelaxationNetworkRefusal {
    fn from(error: RelaxationNetworkMotionError) -> Self {
        match error {
            RelaxationNetworkMotionError::Cancelled => Self::Cancelled,
            error => Self::Motion(error),
        }
    }
}

impl From<DynamicsError> for RelaxationNetworkRefusal {
    fn from(error: DynamicsError) -> Self {
        Self::Dynamics(error)
    }
}

/// The bodies and the complete elastic memory of a [`RelaxationNetwork`].
#[derive(Debug, Clone, PartialEq)]
pub struct RelaxationState {
    /// One state per body, in the network's order.
    pub bodies: Vec<RigidBodyState>,
    /// One amplitude per branch of every connector, in the network's order
    /// [sqrt(J)].
    pub memory: RelaxationNetworkMemory,
}

/// Rigid bodies joined by connectors with Kelvin viscosity and Maxwell
/// memory, stepped by their owner.
#[derive(Debug, Clone)]
pub struct RelaxationNetwork {
    properties: Vec<MassProperties>,
    links: Vec<RelaxationNetworkLink>,
    gravity: Gravity,
    controls: RelaxationNetworkMotionControls,
    /// The actuators, in the order of their ports.
    drives: Vec<Drive>,
    /// What the temperature signals set, and where the first of them is in
    /// the input vector.
    thermal: Thermal,
    temperature_input: usize,
    /// The body torque ports, in the order of their ports.
    torques: Vec<super::TorquePort>,
    ports: Vec<PortSpec>,
}

impl RelaxationNetwork {
    /// A network of `properties.len()` bodies and `links.len()` connectors
    /// under uniform `gravity`, stepped with the owner's `controls`.
    ///
    /// # Errors
    /// [`RelaxationNetworkBuildError`] for no bodies or a connector that
    /// does not join two different bodies of the network.
    pub fn new(
        properties: Vec<MassProperties>,
        links: Vec<RelaxationNetworkLink>,
        gravity: Gravity,
        controls: RelaxationNetworkMotionControls,
    ) -> Result<Self, RelaxationNetworkBuildError> {
        let bodies = properties.len();
        if bodies == 0 {
            return Err(RelaxationNetworkBuildError::NoBodies);
        }
        for (connector, link) in links.iter().enumerate() {
            let [first, second] = link.viscous.bodies;
            if first >= bodies || second >= bodies || first == second {
                return Err(RelaxationNetworkBuildError::Link { connector });
            }
        }
        let three = NonZeroUsize::MIN.saturating_add(2);
        let mut ports = Vec::with_capacity(bodies + links.len());
        for body in 0..bodies {
            let id = StableId::new(format!("force{body}"))
                .map_err(RelaxationNetworkBuildError::Identifier)?;
            ports.push(
                PortSpec::power(
                    id,
                    PortKind::MechanicalForceVelocity,
                    PortRole::EffortIn,
                    three,
                )
                .map_err(RelaxationNetworkBuildError::Port)?,
            );
        }
        for connector in 0..links.len() {
            let id = StableId::new(format!("heat{connector}"))
                .map_err(RelaxationNetworkBuildError::Identifier)?;
            ports.push(
                PortSpec::heat(id, PortRole::EffortIn)
                    .map_err(RelaxationNetworkBuildError::Port)?,
            );
        }
        Ok(Self {
            properties,
            links,
            gravity,
            controls,
            drives: Vec::new(),
            thermal: Thermal::Fixed,
            temperature_input: 0,
            torques: Vec::new(),
            ports,
        })
    }

    /// Coefficients the ports declared so far take in the input vector and
    /// in the output vector.
    fn widths(&self) -> (usize, usize) {
        super::widths(&self.ports)
    }

    /// Gives `body` a torque port `torque<body>`, appended after the ports
    /// the network already has.
    ///
    /// # Errors
    /// [`RelaxationNetworkBuildError::Torque`] for a body that does not
    /// exist or already has the port.
    pub fn torqued(mut self, body: usize) -> Result<Self, RelaxationNetworkBuildError> {
        if body >= self.properties.len() || self.torques.iter().any(|port| port.body == body) {
            return Err(RelaxationNetworkBuildError::Torque { body });
        }
        let port = super::TorquePort::declare(body, &mut self.ports)?;
        self.torques.push(port);
        Ok(self)
    }

    /// Appends one `temperature<e>` signal input per connector, read by
    /// `thermal`.
    fn following(mut self, thermal: Thermal) -> Result<Self, RelaxationNetworkBuildError> {
        if !matches!(self.thermal, Thermal::Fixed) {
            return Err(RelaxationNetworkBuildError::Thermal);
        }
        self.temperature_input = self.widths().0;
        for connector in 0..self.links.len() {
            let id = StableId::new(format!("temperature{connector}"))
                .map_err(RelaxationNetworkBuildError::Identifier)?;
            self.ports.push(PortSpec::signal_in(id, NonZeroUsize::MIN));
        }
        self.thermal = thermal;
        Ok(self)
    }

    /// The same network with Kelvin damping that follows temperature: each
    /// connector's Kelvin form for a step is `law` at the temperature its
    /// `temperature<e>` signal input receives, from the form the network
    /// was built with. The Maxwell branches are not changed.
    ///
    /// # Errors
    /// [`RelaxationNetworkBuildError::Thermal`] for a network that already
    /// follows temperature.
    pub fn temperature_dependent(
        self,
        law: Rc<dyn DampingLaw>,
    ) -> Result<Self, RelaxationNetworkBuildError> {
        self.following(Thermal::Damping(law))
    }

    /// The same network with kinetics that follow temperature: each
    /// connector's Kelvin viscosity and every relaxation time of its
    /// branches are, for a step, the built ones times `law`'s factor at the
    /// temperature its `temperature<e>` signal input receives.
    ///
    /// # Errors
    /// [`RelaxationNetworkBuildError::Thermal`] for a network that already
    /// follows temperature.
    pub fn shifted(self, law: Rc<dyn ShiftLaw>) -> Result<Self, RelaxationNetworkBuildError> {
        self.following(Thermal::Kinetics(law))
    }

    /// Gives `connector` an actuator: the ports `drive<connector>` (held
    /// force) and `couple<connector>` (held couple), appended after the
    /// ports the network already has, with their efforts held in `frame`.
    ///
    /// # Errors
    /// [`RelaxationNetworkBuildError::Drive`] for a connector that does not
    /// exist or already has an actuator.
    pub fn driven(
        mut self,
        connector: usize,
        frame: DriveFrame,
    ) -> Result<Self, RelaxationNetworkBuildError> {
        if connector >= self.links.len()
            || self.drives.iter().any(|drive| drive.connector == connector)
        {
            return Err(RelaxationNetworkBuildError::Drive { connector });
        }
        let (input, output) = self.widths();
        let three = NonZeroUsize::MIN.saturating_add(2);
        for (name, kind) in [
            ("drive", PortKind::MechanicalForceVelocity),
            ("couple", PortKind::RotationalTorqueAngularVelocity),
        ] {
            let id = StableId::new(format!("{name}{connector}"))
                .map_err(RelaxationNetworkBuildError::Identifier)?;
            self.ports.push(
                PortSpec::power(id, kind, PortRole::EffortIn, three)
                    .map_err(RelaxationNetworkBuildError::Port)?,
            );
        }
        self.drives.push(Drive {
            connector,
            frame,
            input,
            output,
        });
        Ok(self)
    }

    /// The state for `bodies` with every memory relaxed (zero).
    #[must_use]
    pub fn relaxed(&self, bodies: Vec<RigidBodyState>) -> RelaxationState {
        RelaxationState {
            bodies,
            memory: self
                .links
                .iter()
                .map(|link| vec![[0.0; 6]; link.branches.len()])
                .collect(),
        }
    }

    /// Number of bodies; ports `0..bodies()` are their force ports.
    #[must_use]
    pub fn bodies(&self) -> usize {
        self.properties.len()
    }

    /// Number of connectors; ports `bodies()..bodies() + connectors()` are
    /// their heat ports. The two ports of each actuator, the temperature
    /// signals and the body torque ports follow, in the order they were
    /// asked for.
    #[must_use]
    pub fn connectors(&self) -> usize {
        self.links.len()
    }

    /// The effort each connector holds during a step: zero without an
    /// actuator, and what the actuator's two ports received (force, then
    /// couple) with one.
    fn efforts(&self, inputs: &[f64]) -> Result<Vec<HeldInternalEffort>, RelaxationNetworkRefusal> {
        let mut efforts = vec![HeldInternalEffort::MaterialFirst([0.0; 6]); self.links.len()];
        for &Drive {
            connector,
            frame,
            input,
            ..
        } in &self.drives
        {
            let held: [f64; 6] = core::array::from_fn(|c| inputs[input + c]);
            if held.iter().any(|value| !value.is_finite()) {
                return Err(RelaxationNetworkRefusal::Effort { connector });
            }
            efforts[connector] = match frame {
                DriveFrame::MaterialFirst => HeldInternalEffort::MaterialFirst(held),
                DriveFrame::World => HeldInternalEffort::World(held),
            };
        }
        Ok(efforts)
    }

    /// Writes to `outputs` the rate conjugate to each actuator's held
    /// effort, and returns the size [J] of the drive ports' part of the
    /// energy allowance.
    ///
    /// The owner applies an effort in two sources of half a step each, at
    /// the start pose and at the drifted pose, and reports for each the
    /// mean relative twist in the first body's material axes at that pose;
    /// its work is the effort (in those axes) times that twist times half a
    /// step, summed. For an effort held in material axes the conjugate rate
    /// is therefore the mean of the two twists; for one held in world axes,
    /// the mean of the two twists turned to world axes at their own poses.
    /// A source changes momenta and no pose, so the second source's pose is
    /// that of the state the step returns.
    ///
    /// The two half-step works of a connector cancel when the relative
    /// motion reverses between the halves, while the rounding of each
    /// stays; so the allowance is sized by each half and not by their sum,
    /// and for each port by its own pair (a force with the linear rate, a
    /// couple with the angular rate).
    fn drive_rates(
        &self,
        step: &RelaxationNetworkMotionStep,
        old: &[RigidBodyState],
        dt_s: f64,
        inputs: &[f64],
        outputs: &mut [f64],
    ) -> f64 {
        let after = step.states_after();
        let [initial, endpoint] = &step.half_sources;
        let mut allowance_j = 0.0_f64;
        for &Drive {
            connector,
            frame,
            input,
            output,
        } in &self.drives
        {
            let twists = [
                initial.average_relative_twist[connector],
                endpoint.average_relative_twist[connector],
            ];
            let rate: [f64; 6] = match frame {
                DriveFrame::MaterialFirst => {
                    core::array::from_fn(|c| 0.5 * (twists[0][c] + twists[1][c]))
                }
                DriveFrame::World => {
                    let first = self.links[connector].viscous.bodies[0];
                    let poses = [old[first], after[first]];
                    let turned: [[f64; 6]; 2] = core::array::from_fn(|source| {
                        let rotation = poses[source].pose().orientation();
                        let t = twists[source];
                        let linear = rotation.rotate_body_to_world(Vec3::new(t[0], t[1], t[2]));
                        let angular = rotation.rotate_body_to_world(Vec3::new(t[3], t[4], t[5]));
                        [linear.x, linear.y, linear.z, angular.x, angular.y, angular.z]
                    });
                    core::array::from_fn(|c| 0.5 * (turned[0][c] + turned[1][c]))
                }
            };
            // The drive port, then the couple port.
            for port in 0..2 {
                let (mut effort_sq, mut twist_sq) = (0.0_f64, [0.0_f64; 2]);
                for c in 3 * port..3 * port + 3 {
                    outputs[output + c] = rate[c];
                    effort_sq += inputs[input + c] * inputs[input + c];
                    twist_sq[0] += twists[0][c] * twists[0][c];
                    twist_sq[1] += twists[1][c] * twists[1][c];
                }
                // Turning a twist to world axes does not change the length
                // of its linear part or of its angular part.
                allowance_j +=
                    0.5 * dt_s * effort_sq.sqrt() * (twist_sq[0].sqrt() + twist_sq[1].sqrt());
            }
            allowance_j += initial.connector_actuator_work_j[connector].abs()
                + endpoint.connector_actuator_work_j[connector].abs();
        }
        allowance_j
    }

    /// The temperature `connector`'s signal carries.
    fn signal(&self, connector: usize, inputs: &[f64]) -> Result<f64, RelaxationNetworkRefusal> {
        let kelvin = inputs[self.temperature_input + connector];
        if kelvin.is_finite() && kelvin > 0.0 {
            Ok(kelvin)
        } else {
            Err(RelaxationNetworkRefusal::DampingTemperature { connector, kelvin })
        }
    }

    /// The owner's step from `old` under `wrenches` and the efforts and
    /// temperatures in `inputs`: with the built connectors, with Kelvin
    /// forms from the damping law, or with kinetics from the shift law.
    fn step(
        &self,
        old: &RelaxationState,
        wrenches: &[Wrench],
        inputs: &[f64],
        dt_s: f64,
        cx: &Cx<'_>,
    ) -> Result<RelaxationNetworkMotionStep, RelaxationNetworkRefusal> {
        let efforts = self.efforts(inputs)?;
        let loads = HeldExternalLoads {
            wrenches,
            gravity: self.gravity,
        };
        let cancelled = || cx.is_cancel_requested();
        let refused = |connector: usize| {
            move |reason: String| RelaxationNetworkRefusal::Law { connector, reason }
        };
        let step = match &self.thermal {
            Thermal::Fixed => step_relaxation_network_with_loads(
                &old.bodies,
                &self.properties,
                &self.links,
                &old.memory,
                &efforts,
                loads,
                dt_s,
                self.controls,
                cancelled,
            )?,
            Thermal::Damping(law) => {
                // Always from the built form, never from a shifted one.
                let mut links = self.links.clone();
                for (connector, link) in links.iter_mut().enumerate() {
                    let kelvin = self.signal(connector, inputs)?;
                    link.viscous.damping = law
                        .at_temperature(connector, &self.links[connector].viscous.damping, kelvin)
                        .map_err(refused(connector))?;
                }
                step_relaxation_network_with_loads(
                    &old.bodies,
                    &self.properties,
                    &links,
                    &old.memory,
                    &efforts,
                    loads,
                    dt_s,
                    self.controls,
                    cancelled,
                )?
            }
            Thermal::Kinetics(law) => {
                let mut viscosity = Vec::with_capacity(self.links.len());
                let mut times = Vec::with_capacity(self.links.len());
                for (connector, link) in self.links.iter().enumerate() {
                    let kelvin = self.signal(connector, inputs)?;
                    let factor = law.factor(connector, kelvin).map_err(refused(connector))?;
                    // The owner admits a zero viscosity scale and refuses a
                    // zero time scale, so it would take a connector without
                    // branches at a factor of zero and refuse one with a
                    // branch. A factor is a ratio of times either way.
                    if !(factor.is_finite() && factor > 0.0) {
                        return Err(refused(connector)(format!(
                            "a shift factor of {factor} is not a positive ratio of times"
                        )));
                    }
                    viscosity.push(factor);
                    times.push(vec![factor; link.branches.len()]);
                }
                step_relaxation_network_with_loads_and_kinetics(
                    &old.bodies,
                    &self.properties,
                    &self.links,
                    &old.memory,
                    &efforts,
                    HeldRelaxationKinetics {
                        viscosity_scales: &viscosity,
                        relaxation_time_scales: &times,
                    },
                    loads,
                    dt_s,
                    self.controls,
                    cancelled,
                )?
            }
        };
        Ok(step)
    }

    /// The stored energy of one state and the sum of the sizes of its
    /// terms (kinetic terms, gravity potential terms, memory storage), for a
    /// roundoff allowance.
    fn energy(&self, state: &RelaxationState) -> Result<(f64, f64), RelaxationNetworkRefusal> {
        if state.bodies.len() != self.properties.len() {
            return Err(RelaxationNetworkRefusal::Bodies {
                expected: self.properties.len(),
                found: state.bodies.len(),
            });
        }
        if state.memory.len() != self.links.len() {
            return Err(RelaxationNetworkRefusal::Memory {
                connector: state.memory.len().min(self.links.len()),
            });
        }
        let integrator = RigidBodyIntegrator::new(self.gravity);
        let (mut total, mut terms) = (0.0_f64, 0.0_f64);
        for (body, properties) in state.bodies.iter().zip(&self.properties) {
            let diagnostics = integrator.diagnostics(*body, *properties)?;
            total += diagnostics.mechanical_energy;
            terms += diagnostics.translational_kinetic_energy
                + diagnostics.rotational_kinetic_energy
                + super::potential_terms_j(self.gravity, *properties, *body);
        }
        for (connector, (memory, link)) in state.memory.iter().zip(&self.links).enumerate() {
            if memory.len() != link.branches.len() {
                return Err(RelaxationNetworkRefusal::Memory { connector });
            }
            for amplitude in memory {
                let mut storage_j = 0.0_f64;
                for coordinate in amplitude {
                    if !coordinate.is_finite() {
                        return Err(RelaxationNetworkRefusal::Memory { connector });
                    }
                    storage_j += 0.5 * coordinate * coordinate;
                }
                total += storage_j;
                terms += storage_j;
            }
        }
        Ok((total, terms))
    }
}

impl Component for RelaxationNetwork {
    type State = RelaxationState;
    type Error = RelaxationNetworkRefusal;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn stored_energy_j(&self, state: &RelaxationState) -> Result<f64, RelaxationNetworkRefusal> {
        Ok(self.energy(state)?.0)
    }

    fn propose(
        &self,
        old: &RelaxationState,
        dt_s: f64,
        inputs: &[f64],
        new: &mut RelaxationState,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, RelaxationNetworkRefusal> {
        cx.checkpoint()
            .map_err(|_| RelaxationNetworkRefusal::Cancelled)?;
        let bodies = self.properties.len();
        let connectors = self.links.len();
        let (_, terms_before_j) = self.energy(old)?;
        let heat_inputs = &inputs[3 * bodies..3 * bodies + connectors];
        for (connector, &kelvin) in heat_inputs.iter().enumerate() {
            if !(kelvin.is_finite() && kelvin > 0.0) {
                return Err(RelaxationNetworkRefusal::Temperature { connector, kelvin });
            }
        }
        let wrenches = super::held_wrenches(bodies, &self.torques, inputs);
        let step = self.step(old, &wrenches, inputs, dt_s, cx)?;
        let after = step.states_after();

        // Work the kernel will book on the force ports, as a size.
        let mut paired_j = 0.0_f64;
        for body in 0..bodies {
            let start = old.bodies[body].pose().position_world();
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
        // A connector's heat is a sum of source integrals and is never
        // negative; a rounding residue below zero is not entropy destroyed.
        let (mut heat_terms_j, mut entropy_j_per_k) = (0.0_f64, 0.0_f64);
        for (connector, &heat_j) in step.connector_dissipation_j.iter().enumerate() {
            outputs[3 * bodies + connector] = -heat_j / dt_s;
            heat_terms_j += heat_j.abs();
            entropy_j_per_k += heat_j.max(0.0) / heat_inputs[connector];
        }
        let drives_j = self.drive_rates(&step, &old.bodies, dt_s, inputs, outputs);
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
        // Kinetic, gravity and memory storage together.
        let stored_change_j = step.mechanical_energy_change_j;
        let residual_j = step.work_heat_residual_j;
        let dissipated_j = step.dissipated_energy_j.max(0.0);

        new.bodies.clear();
        new.bodies.extend_from_slice(after);
        new.memory.clone_from(step.memory_after());
        let (_, terms_after_j) = self.energy(new)?;

        // The tolerance covers the kernel's f64 port products (each
        // half-step's, for a drive port), the roundoff of `stored_energy_j`
        // at both ends, and the owner's rounded sums.
        let scale_j = paired_j
            + heat_terms_j
            + terms_before_j
            + terms_after_j
            + stored_change_j.abs()
            + residual_j.abs()
            + drives_j
            + torques_j;
        Ok(Proposal::Step(Receipt {
            stored_change_j,
            dissipated_j,
            // The residual is `delta(K + V + memory) + heat - actuator work
            // - applied work`; a positive value is energy the splitting
            // created. The drive ports carry the actuator work in; a
            // connector without them holds a zero effort, whose work the
            // owner computes as exactly zero.
            numerical_loss_j: -residual_j,
            entropy_generated_j_per_k: entropy_j_per_k,
            energy_tolerance_j: 64.0 * f64::EPSILON * scale_j + f64::MIN_POSITIVE,
            residual: 0.0,
            tolerance: 0.0,
        }))
    }
}
