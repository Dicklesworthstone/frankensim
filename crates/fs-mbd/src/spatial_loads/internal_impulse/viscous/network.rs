//! Sparse networks of full spatial viscous connectors on finite native bodies.
//!
//! Every edge delegates its exact fixed-geometry relaxation, wrench transport,
//! native momentum update and analytic heat to PreparedViscousConnector.
//! A palindromic composition resolves shared-body coupling to second order
//! under subdivision. It is not the simultaneous network exponential.
//! The supplied edge order is part of the finite-step method; cycles, parallel
//! connectors and isolated bodies are admitted without a dense global matrix.
//!
//! Geometry, inertia, damping and time interval are bound at preparation.
//! No clock, pose, thermal state, spring or constraint is advanced here.

/// Simultaneous held internal actuation and full spatial damping.
#[path = "network/driven.rs"]
pub mod driven;

/// Complete material-frame Maxwell histories on the same native graph.
pub use driven::relaxation;


use super::super::{
    InternalWrenchImpulse, InternalWrenchImpulseError, WideVec, add, cross, dd, finite,
    finite_vector,
    network::{InternalWrenchNetworkError, admit_limit, reserve},
    norm, poll, prepare, rotate_wide, rounded, subtract, validate_controls,
};
use super::{
    PassiveBushingMatrix, PreparedViscousConnector, ViscousConnectorControls, ViscousConnectorError,
};
use crate::{DynamicsError, MassProperties, Pose, RigidBodyState, RigidKineticEnergyChange, Vec3};
use fs_math::dd::Dd;

/// One original material connector in the supplied deterministic edge order.
#[derive(Clone, Debug, PartialEq)]
pub struct ViscousNetworkLink {
    /// Distinct original body indices. Shared bodies and cycles are permitted.
    pub bodies: [usize; 2],
    /// COM-relative attachment arms in each body's own material coordinates [m].
    pub arms_body_m: [Vec3; 2],
    /// Full original passive SI damping form, expressed in the first body's
    /// material axes and acting on the transported relative twist at the second.
    pub damping: PassiveBushingMatrix,
}

/// Hard resource caps checked before allocating the prepared source.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ViscousNetworkLimits {
    pub maximum_bodies: usize,
    pub maximum_connectors: usize,
    /// Actual scheduled pair calls: (2 * active_connectors - 1) * substeps.
    /// Zero duration and zero damping require no pair calls.
    pub maximum_pair_relaxations: usize,
    /// Conservative logical owned heap payload for this preparation plus one
    /// relaxation and its returned state/heat arrays [bytes]. Excludes borrowed
    /// inputs, stack storage, allocator metadata and capacity rounding.
    pub maximum_live_bytes: usize,
}

/// Accepted complete network source, in the original body and connector order.
#[derive(Clone, Debug, PartialEq)]
#[must_use]
pub struct ViscousNetworkReceipt {
    pub states_after: Vec<RigidBodyState>,
    /// Analytic heat from each original connector, including exact zero entries
    /// for inactive edges. No thermal destination or partition is inferred [J].
    pub connector_dissipation_j: Vec<f64>,
    pub source_duration_s: f64,
    pub substeps: usize,
    pub pair_relaxations: usize,
    pub angular_origin_world_m: Vec3,
    pub linear_momentum_change_world_ns: Vec3,
    pub angular_momentum_change_world_nms: Vec3,
    /// All original bodies' actual endpoint kinetic increments, reduced before
    /// rounding so opposing transfers cannot hide the physical loss [J].
    pub kinetic_energy_change_j: f64,
    pub dissipated_energy_j: f64,
    /// Sum of independently recovered native connector works [J].
    pub connector_work_j: f64,
    pub kinetic_heat_residual_j: f64,
    pub connector_heat_residual_j: f64,
    /// Maximum of every body's accumulated intended/actual transfer residual
    /// and the complete network's actual net momentum change.
    pub maximum_linear_transfer_residual_ns: f64,
    pub maximum_angular_transfer_residual_nms: f64,
    pub maximum_linear_speed_residual_m_s: f64,
    pub maximum_angular_speed_residual_rad_s: f64,
    pub planned_live_bytes: usize,
}

/// An error returns no partially relaxed network; all inputs are borrowed.
#[derive(Clone, Debug, PartialEq)]
pub enum ViscousNetworkError {
    Invalid(&'static str),
    Resource(InternalWrenchNetworkError),
    Impulse(InternalWrenchImpulseError),
    Connector {
        /// Original supplied connector index, including inactive entries.
        connector: usize,
        error: ViscousConnectorError,
    },
    Unresolved {
        linear_residual_ns: f64,
        angular_residual_nms: f64,
        kinetic_heat_residual_j: f64,
        connector_heat_residual_j: f64,
    },
    Cancelled,
}
impl From<InternalWrenchNetworkError> for ViscousNetworkError {
    fn from(value: InternalWrenchNetworkError) -> Self {
        match value {
            InternalWrenchNetworkError::Cancelled => Self::Cancelled,
            value => Self::Resource(value),
        }
    }
}
impl From<InternalWrenchImpulseError> for ViscousNetworkError {
    fn from(value: InternalWrenchImpulseError) -> Self {
        match value {
            InternalWrenchImpulseError::Cancelled => Self::Cancelled,
            value => Self::Impulse(value),
        }
    }
}
impl From<DynamicsError> for ViscousNetworkError {
    fn from(value: DynamicsError) -> Self {
        Self::Impulse(value.into())
    }
}
impl core::fmt::Display for ViscousNetworkError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "viscous connector network refused: {self:?}")
    }
}
impl std::error::Error for ViscousNetworkError {}

#[derive(Clone, Debug)]
struct PreparedLink {
    original_index: usize,
    bodies: [usize; 2],
    arms_body_m: [Vec3; 2],
    source: PreparedViscousConnector,
}

#[derive(Clone, Copy)]
struct IntendedImpulse {
    linear: WideVec,
    central_world: WideVec,
}
impl IntendedImpulse {
    const ZERO: Self = Self {
        linear: [Dd::ZERO; 3],
        central_world: [Dd::ZERO; 3],
    };
}

/// Prepared sparse passive source with exact native finite-body pair recoil.
///
/// One small source is prepared per active connector. All but the final active
/// edge use half a substep; the central edge uses a full substep. Forward and
/// reverse sweeps form a symmetric composition. Pair work is O(active_edges *
/// substeps), with O(bodies + supplied_edges) admission and endpoint accounting.
/// Storage is O(bodies + supplied_edges), without a global exponential or a
/// repeated eigensolve in the source loop.
///
/// Reuse is permitted only at these exact poses, inertias, damping forms and
/// interval. Momentum may change freely between calls. Composition error is
/// controlled by the explicitly supplied substep count, independently of the
/// physical momentum, work, heat and endpoint-rate admission tolerances.
#[derive(Clone, Debug)]
pub struct PreparedViscousNetwork {
    poses: Vec<Pose>,
    properties: Vec<MassProperties>,
    links: Vec<PreparedLink>,
    connector_count: usize,
    duration_s: f64,
    substeps: usize,
    pair_relaxations: usize,
    planned_live_bytes: usize,
}

impl PreparedViscousNetwork {
    /// Prepare the original source. Empty edge sets and zero source intervals
    /// are identities; isolated bodies retain their exact native state.
    /// Substeps must be positive even for an identity source. No coefficient is
    /// fitted or averaged and no positive damping direction is discarded.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        poses: &[Pose],
        properties: &[MassProperties],
        links: &[ViscousNetworkLink],
        duration_s: f64,
        substeps: usize,
        limits: ViscousNetworkLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ViscousNetworkError> {
        poll(&mut cancelled)?;
        if poses.len() != properties.len()
            || !duration_s.is_finite()
            || duration_s < 0.0
            || substeps == 0
        {
            return Err(ViscousNetworkError::Invalid(
                "matching bodies/inertias, finite nonnegative duration and positive substeps required",
            ));
        }
        admit_limit("bodies", poses.len() as u128, limits.maximum_bodies as u128)?;
        admit_limit(
            "connectors",
            links.len() as u128,
            limits.maximum_connectors as u128,
        )?;
        let mut active = 0;
        for link in links {
            poll(&mut cancelled)?;
            let [a, b] = link.bodies;
            if a == b || a >= poses.len() || b >= poses.len() {
                return Err(ViscousNetworkError::Invalid(
                    "connector requires two distinct existing bodies",
                ));
            }
            for arm in link.arms_body_m {
                finite_vector(arm)?;
            }
            if duration_s > 0.0 && nonzero(&link.damping) {
                active += 1;
            }
        }
        let calls = if active == 0 {
            0
        } else {
            (2 * active as u128 - 1)
                .checked_mul(substeps as u128)
                .ok_or(ViscousNetworkError::Invalid(
                    "pair-relaxation count overflows",
                ))?
        };
        admit_limit(
            "pair relaxations",
            calls,
            limits.maximum_pair_relaxations as u128,
        )?;
        let bytes = planned_live_bytes(poses.len(), links.len(), active);
        admit_limit("live bytes", bytes, limits.maximum_live_bytes as u128)?;
        let full_step = if active == 0 {
            0.0
        } else {
            duration_s / substeps as f64
        };
        if active != 0
            && (full_step == 0.0
                || (active > 1 && (full_step * 0.5 == 0.0 || 2.0 * (0.5 * full_step) != full_step)))
        {
            return Err(ViscousNetworkError::Invalid(
                "positive source subdivision is unrepresentable",
            ));
        }
        let mut bound_poses = reserve(poses.len())?;
        let mut bound_properties = reserve(properties.len())?;
        let mut prepared = reserve(active)?;
        for (&pose, &properties) in poses.iter().zip(properties) {
            poll(&mut cancelled)?;
            let state = RigidBodyState::new(pose, Vec3::ZERO, Vec3::ZERO)?;
            state.kinetic_energy_change_to(properties, state)?;
            bound_poses.push(pose);
            bound_properties.push(properties);
        }
        for (original_index, link) in links.iter().enumerate() {
            poll(&mut cancelled)?;
            if duration_s == 0.0 || !nonzero(&link.damping) {
                continue;
            }
            let step = if prepared.len() + 1 == active {
                full_step
            } else {
                full_step * 0.5
            };
            let source = PreparedViscousConnector::new(
                link.bodies.map(|i| poses[i]),
                link.bodies.map(|i| properties[i]),
                link.arms_body_m,
                link.damping.clone(),
                step,
                &mut cancelled,
            )
            .map_err(|error| connector_error(original_index, error))?;
            prepared.push(PreparedLink {
                original_index,
                bodies: link.bodies,
                arms_body_m: link.arms_body_m,
                source,
            });
        }
        poll(&mut cancelled)?;
        Ok(Self {
            poses: bound_poses,
            properties: bound_properties,
            links: prepared,
            connector_count: links.len(),
            duration_s,
            substeps,
            pair_relaxations: calls as usize,
            planned_live_bytes: bytes as usize,
        })
    }

    /// Original fixed-geometry source interval [s].
    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }
    /// Explicit number of symmetric source subdivisions.
    pub const fn substeps(&self) -> usize {
        self.substeps
    }
    /// Nonzero-law connectors whose original interval is nonzero.
    pub fn active_connectors(&self) -> usize {
        self.links.len()
    }
    /// Complete scheduled native pair calls after inactive-edge compaction.
    pub const fn pair_relaxations(&self) -> usize {
        self.pair_relaxations
    }
    /// Admitted logical peak payload for the prepared owner and one relaxation.
    pub const fn planned_live_bytes(&self) -> usize {
        self.planned_live_bytes
    }

    /// Stage all native body momenta, then admit the original pair laws and the
    /// whole-network balances before returning an accepted state bank.
    ///
    /// The supplied pair impulse tolerances also bound every body's accumulated
    /// command/readback defect, total P/L and both network heat/work residuals.
    /// Each pair's endpoint twist retains the original source tolerances.
    /// A failed later edge or cancelled sweep exposes no intermediate state.
    pub fn relax(
        &self,
        states: &[RigidBodyState],
        controls: ViscousConnectorControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ViscousNetworkReceipt, ViscousNetworkError> {
        poll(&mut cancelled)?;
        validate_controls(controls.impulse)?;
        if [
            controls.linear_speed_tolerance_m_s,
            controls.angular_speed_tolerance_rad_s,
        ]
        .into_iter()
        .any(|v| !v.is_finite() || v < 0.0)
        {
            return Err(ViscousNetworkError::Invalid(
                "finite nonnegative twist tolerances required",
            ));
        }
        if states.len() != self.poses.len()
            || states
                .iter()
                .zip(&self.poses)
                .any(|(s, pose)| s.pose() != *pose)
        {
            return Err(ViscousNetworkError::Invalid(
                "changed prepared body count or poses",
            ));
        }
        let mut staged = reserve(states.len())?;
        staged.extend_from_slice(states);
        let mut intended = reserve(states.len())?;
        intended.resize(states.len(), IntendedImpulse::ZERO);
        let mut edge_heat = reserve(self.connector_count)?;
        edge_heat.resize(self.connector_count, Dd::ZERO);
        let mut heat = Dd::ZERO;
        let mut work = Dd::ZERO;
        let mut linear_rate_error = 0.0_f64;
        let mut angular_rate_error = 0.0_f64;
        let mut exchange = |link: &PreparedLink,
                            cancelled: &mut dyn FnMut() -> bool|
         -> Result<(), ViscousNetworkError> {
            let pair = link.bodies.map(|i| staged[i]);
            let receipt = link
                .source
                .relax(pair, controls, cancelled)
                .map_err(|error| connector_error(link.original_index, error))?;
            // Recover the same unrounded intended moment transport used by the
            // pair owner. Summing body-frame rounded diagnostics alone would
            // discard frame/transport defects across repeated exchanges.
            let geometry = prepare(
                pair,
                link.bodies.map(|i| self.properties[i]),
                InternalWrenchImpulse {
                    arms_body_m: link.arms_body_m,
                    impulse_on_second_world_ns: receipt.impulse.bodies[1].linear_impulse_world,
                    couple_on_second_world_nms: receipt.impulse.attachment_couples_world_nms[1],
                },
            )?;
            let force = dd(receipt.impulse.bodies[1].linear_impulse_world);
            for (side, &body) in link.bodies.iter().enumerate() {
                intended[body].linear = if side == 0 {
                    subtract(intended[body].linear, force)
                } else {
                    add(intended[body].linear, force)
                };
                intended[body].central_world =
                    add(intended[body].central_world, geometry.central[side]);
                staged[body] = receipt.impulse.bodies[side].state_after;
            }
            let loss = Dd::from_f64(receipt.dissipated_energy_j);
            edge_heat[link.original_index] = edge_heat[link.original_index] + loss;
            heat = heat + loss;
            work = work + Dd::from_f64(receipt.impulse.connector_work_j);
            for value in receipt.relative_twist_residual[..3].iter() {
                linear_rate_error = linear_rate_error.max(value.abs());
            }
            for value in receipt.relative_twist_residual[3..].iter() {
                angular_rate_error = angular_rate_error.max(value.abs());
            }
            Ok(())
        };
        if !self.links.is_empty() {
            for _ in 0..self.substeps {
                for link in &self.links {
                    poll(&mut cancelled)?;
                    exchange(link, &mut cancelled)?;
                }
                for link in self.links[..self.links.len() - 1].iter().rev() {
                    poll(&mut cancelled)?;
                    exchange(link, &mut cancelled)?;
                }
            }
        }
        let mut energy = RigidKineticEnergyChange::new();
        let mut linear = [Dd::ZERO; 3];
        let mut angular = [Dd::ZERO; 3];
        let mut linear_error = 0.0_f64;
        let mut angular_error = 0.0_f64;
        for (i, (&before, &after)) in states.iter().zip(&staged).enumerate() {
            poll(&mut cancelled)?;
            energy.add_body(self.properties[i], before, after)?;
            let dp = subtract(
                dd(after.linear_momentum_world()),
                dd(before.linear_momentum_world()),
            );
            let dh = subtract(
                dd(after.angular_momentum_body()),
                dd(before.angular_momentum_body()),
            );
            let axes = [
                Vec3::new(1., 0., 0.),
                Vec3::new(0., 1., 0.),
                Vec3::new(0., 0., 1.),
            ]
            .map(|axis| before.pose().orientation().rotate_body_to_world(axis));
            let spin = rotate_wide(axes, dh);
            let arm = subtract(
                dd(before.pose().position_world()),
                dd(controls.impulse.angular_origin_world_m),
            );
            linear = add(linear, dp);
            angular = add(angular, add(cross(arm, dp), spin));
            linear_error = linear_error.max(norm(rounded(subtract(dp, intended[i].linear))?));
            angular_error =
                angular_error.max(norm(rounded(subtract(spin, intended[i].central_world))?));
        }
        let energy = energy.finish_retained()?;
        let linear_momentum_change_world_ns = rounded(linear)?;
        let angular_momentum_change_world_nms = rounded(angular)?;
        linear_error = linear_error.max(norm(linear_momentum_change_world_ns));
        angular_error = angular_error.max(norm(angular_momentum_change_world_nms));
        let kinetic_heat_residual_j = finite((energy + heat).to_f64())?;
        let connector_heat_residual_j = finite((work + heat).to_f64())?;
        if !linear_error.is_finite()
            || !angular_error.is_finite()
            || linear_error > controls.impulse.linear_tolerance_ns
            || angular_error > controls.impulse.angular_tolerance_nms
            || kinetic_heat_residual_j.abs() > controls.impulse.work_tolerance_j
            || connector_heat_residual_j.abs() > controls.impulse.work_tolerance_j
        {
            return Err(ViscousNetworkError::Unresolved {
                linear_residual_ns: linear_error,
                angular_residual_nms: angular_error,
                kinetic_heat_residual_j,
                connector_heat_residual_j,
            });
        }
        let mut connector_dissipation_j = reserve(self.connector_count)?;
        for loss in edge_heat {
            poll(&mut cancelled)?;
            connector_dissipation_j.push(finite(loss.to_f64())?);
        }
        let result = ViscousNetworkReceipt {
            states_after: staged,
            connector_dissipation_j,
            source_duration_s: self.duration_s,
            substeps: self.substeps,
            pair_relaxations: self.pair_relaxations,
            angular_origin_world_m: controls.impulse.angular_origin_world_m,
            linear_momentum_change_world_ns,
            angular_momentum_change_world_nms,
            kinetic_energy_change_j: finite(energy.to_f64())?,
            dissipated_energy_j: finite(heat.to_f64())?,
            connector_work_j: finite(work.to_f64())?,
            kinetic_heat_residual_j,
            connector_heat_residual_j,
            maximum_linear_transfer_residual_ns: linear_error,
            maximum_angular_transfer_residual_nms: angular_error,
            maximum_linear_speed_residual_m_s: linear_rate_error,
            maximum_angular_speed_residual_rad_s: angular_rate_error,
            planned_live_bytes: self.planned_live_bytes,
        };
        poll(&mut cancelled)?;
        Ok(result)
    }
}

// Shared with the finite-motion consumer so the complete trajectory admits
// its peak before allocating either original source preparation.
pub(in crate::spatial_loads) fn planned_live_bytes(
    bodies: usize,
    connectors: usize,
    active: usize,
) -> u128 {
    let retained = bodies as u128
        * (core::mem::size_of::<Pose>() + core::mem::size_of::<MassProperties>()) as u128
        + active as u128 * core::mem::size_of::<PreparedLink>() as u128;
    let staged = bodies as u128
        * (core::mem::size_of::<RigidBodyState>() + core::mem::size_of::<IntendedImpulse>())
            as u128
        + connectors as u128 * (core::mem::size_of::<Dd>() + core::mem::size_of::<f64>()) as u128;
    // Original six-coordinate Jacobi owner's conservative aggregate
    // 4*n*n+3*n envelope, including index storage at the wider scalar size.
    let spectral = if active == 0 {
        0
    } else {
        (4 * 6 * 6 + 3 * 6) * core::mem::size_of::<f64>().max(core::mem::size_of::<usize>()) as u128
    };
    retained + staged.max(spectral)
}

pub(in crate::spatial_loads) fn nonzero(damping: &PassiveBushingMatrix) -> bool {
    damping
        .native_coefficients()
        .iter()
        .flatten()
        .any(|v| *v != 0.0)
}
fn connector_error(connector: usize, error: ViscousConnectorError) -> ViscousNetworkError {
    match error {
        ViscousConnectorError::Impulse(InternalWrenchImpulseError::Cancelled) => {
            ViscousNetworkError::Cancelled
        }
        error => ViscousNetworkError::Connector { connector, error },
    }
}
