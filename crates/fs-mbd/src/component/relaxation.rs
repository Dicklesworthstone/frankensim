//! A rigid-body network with Maxwell memory as one composition-kernel
//! Component.
//!
//! [`RelaxationNetwork`] wraps [`step_relaxation_network_with_loads`]: rigid
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
//! Ports, in this order, as for [`super::network::ViscousNetwork`]:
//!
//! - `force<i>`, one per body, receiving effort: the world-frame force on
//!   the centre of mass [N], three components. It returns the mean
//!   centre-of-mass velocity of the step, `(x1 - x0) / dt` [m/s].
//! - `heat<e>`, one per connector, a heat port receiving effort: the
//!   temperature at which that connector's heat is released [K]. It returns
//!   minus the connector's mean dissipation rate, `-D_e / dt` [W], Kelvin
//!   and Maxwell loss together.
//!
//! Receipt: the stored change is kinetic plus uniform-gravity potential plus
//! elastic memory, as the owner reduces them; the dissipation is the owner's
//! total connector heat; the owner's work-heat residual is booked as
//! numerical loss; the entropy generated is `sum(D_e / T_e)`. Elastic
//! recovery can speed a body up while its connector still heats: only the
//! independently integrated loss is heat.
//!
//! No claim: the owner's held internal efforts (actuators on a connector)
//! are all zero here, because their work would enter with no port to carry
//! it; no torque port, as for the viscous network; the Kelvin viscosity and
//! the relaxation times are fixed (the owner's held-kinetics entry, which
//! scales those two, is not wrapped); the port map has no analytic tangent;
//! cancellation raised while the owner is inside a step is passed to it but
//! not exercised by a test; the owner's
//! accuracy statements (exact sources at fixed geometry, second order for
//! the moving geometry) are unchanged by the wrapping.

use core::num::NonZeroUsize;

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
    HeldInternalEffort, RelaxationNetworkMotionControls, RelaxationNetworkMotionError,
    step_relaxation_network_with_loads,
};
use crate::{
    DynamicsError, Gravity, MassProperties, RigidBodyIntegrator, RigidBodyState, Vec3, Wrench,
};

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
            Self::Identifier(error) => write!(f, "port identifier refused: {error:?}"),
            Self::Port(violation) => write!(f, "port refused: {violation}"),
        }
    }
}

impl std::error::Error for RelaxationNetworkBuildError {}

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
            ports,
        })
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
    /// their heat ports.
    #[must_use]
    pub fn connectors(&self) -> usize {
        self.links.len()
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
        let (_, terms_before_j) = self.energy(old)?;
        let heat_inputs = &inputs[3 * bodies..];
        for (connector, &kelvin) in heat_inputs.iter().enumerate() {
            if !(kelvin.is_finite() && kelvin > 0.0) {
                return Err(RelaxationNetworkRefusal::Temperature { connector, kelvin });
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
        // No actuator: every connector's held internal effort is zero.
        let idle = vec![HeldInternalEffort::MaterialFirst([0.0; 6]); self.links.len()];
        let step = step_relaxation_network_with_loads(
            &old.bodies,
            &self.properties,
            &self.links,
            &old.memory,
            &idle,
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
        // Kinetic, gravity and memory storage together.
        let stored_change_j = step.mechanical_energy_change_j;
        let residual_j = step.work_heat_residual_j;
        // Zero efforts do no work; whatever the owner rounds to is not
        // energy that entered through a port.
        let actuator_j = step.actuator_work_j;
        let dissipated_j = step.dissipated_energy_j.max(0.0);

        new.bodies.clear();
        new.bodies.extend_from_slice(after);
        new.memory.clone_from(step.memory_after());
        let (_, terms_after_j) = self.energy(new)?;

        // The tolerance covers the kernel's f64 port products, the roundoff
        // of `stored_energy_j` at both ends, and the owner's rounded sums.
        let scale_j = paired_j
            + heat_terms_j
            + terms_before_j
            + terms_after_j
            + stored_change_j.abs()
            + residual_j.abs()
            + actuator_j.abs();
        Ok(Proposal::Step(Receipt {
            stored_change_j,
            dissipated_j,
            // The residual is `delta(K + V + memory) + heat - actuator work
            // - applied work`; a positive value is energy the splitting
            // created.
            numerical_loss_j: -residual_j - actuator_j,
            entropy_generated_j_per_k: entropy_j_per_k,
            energy_tolerance_j: 64.0 * f64::EPSILON * scale_j + f64::MIN_POSITIVE,
            residual: 0.0,
            tolerance: 0.0,
        }))
    }
}
