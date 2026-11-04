//! Simultaneous wrench transport and work on shared finite rigid bodies.

use super::{
    Dd, DynamicsError, InternalWrenchImpulse, InternalWrenchImpulseControls,
    InternalWrenchImpulseError, MassProperties, MomentumImpulseReceipt, RigidBodyState, Vec3,
    WideVec, add, dd, dot, finite, inspect_body, norm, poll, prepare, rounded, subtract,
    validate_controls,
};

/// One internal connector, indexed in the original native body slices.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InternalWrenchLink {
    /// Distinct first/second bodies. Cycles and shared bodies are permitted.
    pub bodies: [usize; 2],
    /// Original force/couple command on the second material attachment.
    pub impulse: InternalWrenchImpulse,
}

/// Whole-network resource bounds and original dimensional balance tolerances.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InternalWrenchNetworkControls {
    /// Per-body intended/actual increments and whole-network momentum/work.
    pub balance: InternalWrenchImpulseControls,
    pub maximum_bodies: usize,
    pub maximum_connectors: usize,
    /// Logical owned heap payload, excluding borrowed inputs, stack storage,
    /// allocator metadata and capacity rounding [bytes].
    pub maximum_live_bytes: usize,
}

/// One connector evaluated against the complete simultaneous body updates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct InternalConnectorReceipt {
    pub bodies: [usize; 2],
    pub attachment_separation_world_m: Vec3,
    pub attachment_couples_world_nms: [Vec3; 2],
    /// Signed work of this connector at the complete before/after midpoint
    /// velocities. Shared-body cross terms belong to every actual source;
    /// no connector is evaluated against a partially updated neighbor.
    pub connector_work_j: f64,
}

/// Atomic event: one original momentum update per body and every source work.
#[derive(Clone, Debug, PartialEq)]
pub struct InternalWrenchNetworkReceipt {
    /// Original body order, including unchanged isolated bodies.
    pub bodies: Vec<MomentumImpulseReceipt>,
    /// Original connector order. Physical results are order-independent up to
    /// the admitted finite arithmetic, without a sequential source selection.
    pub connectors: Vec<InternalConnectorReceipt>,
    pub linear_momentum_change_world_ns: Vec3,
    pub angular_momentum_change_world_nms: Vec3,
    pub kinetic_energy_change_j: f64,
    /// Compensated sum before rounding individual connector diagnostics [J].
    pub connector_work_j: f64,
    pub work_energy_residual_j: f64,
    pub planned_live_bytes: usize,
}

/// No partial network result is published on any refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum InternalWrenchNetworkError {
    Impulse(InternalWrenchImpulseError),
    Invalid(&'static str),
    Limit {
        resource: &'static str,
        required: u128,
        maximum: u128,
    },
    Allocation,
    Cancelled,
}
impl From<InternalWrenchImpulseError> for InternalWrenchNetworkError {
    fn from(error: InternalWrenchImpulseError) -> Self {
        match error {
            InternalWrenchImpulseError::Cancelled => Self::Cancelled,
            error => Self::Impulse(error),
        }
    }
}
impl From<DynamicsError> for InternalWrenchNetworkError {
    fn from(error: DynamicsError) -> Self {
        Self::Impulse(error.into())
    }
}
impl core::fmt::Display for InternalWrenchNetworkError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "internal wrench network refused: {self:?}")
    }
}
impl std::error::Error for InternalWrenchNetworkError {}

#[derive(Clone, Copy)]
struct BodyVectors {
    // Reused after native application for COM and world-spin midpoint velocity.
    linear: WideVec,
    angular: WideVec,
}
impl BodyVectors {
    const ZERO: Self = Self {
        linear: [Dd::ZERO; 3],
        angular: [Dd::ZERO; 3],
    };
}

/// Transfer all internal connector impulses on the original body poses at once.
///
/// The original pair geometry transports the opposite attachment couple for
/// every connector. Wide per-body reductions retain all applied force moments
/// and free couples; each body then calls its original native momentum owner
/// ONCE with that resultant. The same original per-body audit checks actual
/// increments against the unrounded command, so two swallowed opposite body
/// updates cannot masquerade as a successful transfer.
///
/// Connector work uses the complete native before/after midpoint velocities.
/// For example, two impulses on one shared body include their kinetic cross
/// term, with its correct contribution to each connector's work. Independently
/// summing isolated pair energies would lose that term. Net angular momentum
/// includes orbital change and transported couples at separated attachments.
///
/// Work and storage are O(bodies + connectors). One wide body bank is reused
/// for midpoint kinematics; a second bounded connector pass reuses the original
/// geometry instead of retaining another per-edge tensor/kinematic bank. The
/// complete peak payload is admitted before allocation. Inputs are borrowed,
/// and cancellation is polled between bounded native/body/connector operations.
///
/// This is an instantaneous prescribed internal drive, not a contact law,
/// constraint, elastic storage, finite-time force integration or heat source.
/// Positive and negative work require their authored actuator energy owner.
pub fn apply_internal_wrench_network(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    links: &[InternalWrenchLink],
    controls: InternalWrenchNetworkControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<InternalWrenchNetworkReceipt, InternalWrenchNetworkError> {
    poll(&mut cancelled)?;
    validate_controls(controls.balance)?;
    if states.len() != properties.len() {
        return Err(InternalWrenchNetworkError::Invalid(
            "body/property counts differ",
        ));
    }
    admit_limit(
        "bodies",
        states.len() as u128,
        controls.maximum_bodies as u128,
    )?;
    admit_limit(
        "connectors",
        links.len() as u128,
        controls.maximum_connectors as u128,
    )?;
    let bytes = planned_live_bytes(states.len(), links.len());
    admit_limit("live bytes", bytes, controls.maximum_live_bytes as u128)?;
    let mut bank = reserve(states.len())?;
    bank.resize(states.len(), BodyVectors::ZERO);
    let mut bodies = reserve(states.len())?;
    let mut connectors = reserve(links.len())?;
    for link in links {
        poll(&mut cancelled)?;
        let [a, b] = link.bodies;
        if a == b || a >= states.len() || b >= states.len() {
            return Err(InternalWrenchNetworkError::Invalid(
                "connector requires two distinct existing bodies",
            ));
        }
        let prepared = prepare(
            [states[a], states[b]],
            [properties[a], properties[b]],
            link.impulse,
        )?;
        let force = dd(link.impulse.impulse_on_second_world_ns);
        bank[a].linear = subtract(bank[a].linear, force);
        bank[b].linear = add(bank[b].linear, force);
        for (i, central) in link.bodies.into_iter().zip(prepared.central) {
            bank[i].angular = add(bank[i].angular, central);
        }
    }
    let mut linear = [Dd::ZERO; 3];
    let mut angular = [Dd::ZERO; 3];
    let mut energy = Dd::ZERO;
    let (mut linear_error, mut angular_error) = (0.0_f64, 0.0_f64);
    for i in 0..states.len() {
        poll(&mut cancelled)?;
        let expected = bank[i];
        let event = states[i].apply_momentum_impulse(
            properties[i],
            rounded(expected.linear)?,
            states[i]
                .pose()
                .orientation()
                .rotate_world_to_body(rounded(expected.angular)?),
        )?;
        let audit = inspect_body(
            properties[i],
            event,
            expected.linear,
            expected.angular,
            controls.balance.angular_origin_world_m,
        )?;
        linear = add(linear, audit.linear);
        angular = add(angular, audit.angular);
        energy = energy + audit.energy[0] + audit.energy[1];
        linear_error = linear_error.max(audit.linear_error);
        angular_error = angular_error.max(audit.angular_error);
        bank[i] = BodyVectors {
            linear: audit.velocity_mid,
            angular: audit.omega_world_mid,
        };
        bodies.push(event);
    }
    let mut work = Dd::ZERO;
    for link in links {
        poll(&mut cancelled)?;
        let [a, b] = link.bodies;
        let prepared = prepare(
            [states[a], states[b]],
            [properties[a], properties[b]],
            link.impulse,
        )?;
        // Cancel the common COM velocity before adding local spin work.
        let connector = dot(
            subtract(bank[b].linear, bank[a].linear),
            dd(link.impulse.impulse_on_second_world_ns),
        ) + dot(bank[a].angular, prepared.central[0])
            + dot(bank[b].angular, prepared.central[1]);
        work = work + connector;
        connectors.push(InternalConnectorReceipt {
            bodies: link.bodies,
            attachment_separation_world_m: rounded(prepared.separation)?,
            attachment_couples_world_nms: prepared.couples,
            connector_work_j: finite(connector.to_f64())?,
        });
    }
    let linear_momentum_change_world_ns = rounded(linear)?;
    let angular_momentum_change_world_nms = rounded(angular)?;
    linear_error = linear_error.max(norm(linear_momentum_change_world_ns));
    angular_error = angular_error.max(norm(angular_momentum_change_world_nms));
    let work_energy_residual_j = finite((energy - work).to_f64())?;
    if !linear_error.is_finite()
        || !angular_error.is_finite()
        || linear_error > controls.balance.linear_tolerance_ns
        || angular_error > controls.balance.angular_tolerance_nms
        || work_energy_residual_j.abs() > controls.balance.work_tolerance_j
    {
        return Err(InternalWrenchImpulseError::Unresolved {
            linear_residual_ns: linear_error,
            angular_residual_nms: angular_error,
            work_residual_j: work_energy_residual_j,
        }
        .into());
    }
    let result = InternalWrenchNetworkReceipt {
        bodies,
        connectors,
        linear_momentum_change_world_ns,
        angular_momentum_change_world_nms,
        kinetic_energy_change_j: finite(energy.to_f64())?,
        connector_work_j: finite(work.to_f64())?,
        work_energy_residual_j,
        planned_live_bytes: bytes as usize,
    };
    poll(&mut cancelled)?;
    Ok(result)
}

pub(in crate::spatial_loads) fn planned_live_bytes(bodies: usize, links: usize) -> u128 {
    bodies as u128
        * (core::mem::size_of::<MomentumImpulseReceipt>() + core::mem::size_of::<BodyVectors>())
            as u128
        + links as u128 * core::mem::size_of::<InternalConnectorReceipt>() as u128
}
pub(in crate::spatial_loads) fn retained_bytes(bodies: usize, links: usize) -> u128 {
    bodies as u128 * core::mem::size_of::<MomentumImpulseReceipt>() as u128
        + links as u128 * core::mem::size_of::<InternalConnectorReceipt>() as u128
}
pub(in crate::spatial_loads) fn admit_limit(
    resource: &'static str,
    required: u128,
    maximum: u128,
) -> Result<(), InternalWrenchNetworkError> {
    if required > maximum {
        Err(InternalWrenchNetworkError::Limit {
            resource,
            required,
            maximum,
        })
    } else {
        Ok(())
    }
}
pub(in crate::spatial_loads) fn reserve<T>(
    size: usize,
) -> Result<Vec<T>, InternalWrenchNetworkError> {
    let mut v = Vec::new();
    v.try_reserve_exact(size)
        .map_err(|_| InternalWrenchNetworkError::Allocation)?;
    Ok(v)
}
