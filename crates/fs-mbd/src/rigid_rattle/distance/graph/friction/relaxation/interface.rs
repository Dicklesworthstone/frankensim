//! Ordered dry-interface cards on the original native axial friction ports.
//!
//! A uniformly sliding rail uses mu*N. A cylindrical sleeve with one physical
//! radius uses mu*N*r and surface speed r*|omega|. That physical radius is
//! independent of the native hinge's numerical angular length. The complete
//! resolved speed range, including intermediate events, must stay inside the
//! original card's declared temperature/pressure/speed domain.

use super::{
    ArticulatedJoint, JointFrictionFlowControls, JointFrictionFlowError, JointFrictionFlowMotion,
    JointFrictionFlowPort, JointFrictionFlowReport, MAX_JOINT_FRICTION_PORTS, MassProperties,
    RigidBodyState, advance_joint_friction,
};
use fs_tribo::{DryFrictionState, DryInterfaceSystemCard, FrictionLaw, TriboError};

/// The actual source-authored slip geometry, with the native drive convention.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum JointDrySlipGeometry {
    /// Uniform axial slip of the carriage relative to its finite rail.
    Prismatic { drive_n: f64 },
    /// Full circumferential slip at a single physical sleeve radius. This is
    /// not a rolling-element bearing, distributed variable-radius contact,
    /// partial-slip patch or inferred bearing-pressure distribution.
    CylindricalSleeve { radius_m: f64, drive_nm: f64 },
}

/// One fixed-phase mechanical/constitutive binding. Card surfaces A/B follow
/// the joint's original authored body order. Reversing that order requires the
/// caller's corresponding ordered card; no material-pair lookup is inferred.
#[derive(Clone, Copy, Debug)]
pub struct JointDryInterface<'a> {
    pub joint_index: usize,
    pub card: &'a DryInterfaceSystemCard,
    pub geometry: JointDrySlipGeometry,
    /// Held total compressive force [N], supplied by its physical load owner.
    pub normal_force_n: f64,
    /// Positive nominal loaded area [m^2], used with that same force for N/A.
    pub nominal_area_m2: f64,
    /// Held interface temperature [K], not an inferred average of body states.
    pub temperature_k: f64,
}

/// Original card and actual admitted path, in canonical increasing joint order.
#[derive(Clone, Debug)]
pub struct JointDryInterfaceReceipt<'a> {
    pub binding: JointDryInterface<'a>,
    pub ordered_bodies: [usize; 2],
    pub nominal_pressure_pa: f64,
    pub source: JointFrictionFlowPort,
    /// Minimum/maximum actual surface speed, with the held temperature and
    /// nominal pressure. Every affine segment and measured native endpoint is
    /// enclosed by this computed range; it is not an interval-certified tube.
    pub admitted_state_range: [DryFrictionState; 2],
}

/// One atomic native candidate whose complete source-domain checks passed.
#[derive(Clone, Debug)]
#[must_use]
pub struct InterfaceJointFrictionReport<'a> {
    pub mechanical: JointFrictionFlowReport,
    pub interfaces: Vec<JointDryInterfaceReceipt<'a>>,
}

#[derive(Debug)]
pub enum InterfaceJointFrictionError {
    Invalid(&'static str),
    /// The original law needs a different continuous constitutive owner.
    UnsupportedLaw {
        joint_index: usize,
    },
    Constitutive {
        joint_index: usize,
        error: TriboError,
    },
    Mechanical(JointFrictionFlowError),
    Cancelled,
    Allocation,
}
impl core::fmt::Display for InterfaceJointFrictionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "native dry-interface friction: {self:?}")
    }
}
impl std::error::Error for InterfaceJointFrictionError {}

fn checkpoint(cancelled: &mut impl FnMut() -> bool) -> Result<(), InterfaceJointFrictionError> {
    if cancelled() {
        Err(InterfaceJointFrictionError::Cancelled)
    } else {
        Ok(())
    }
}
fn multiply(a: f64, b: f64) -> Result<f64, InterfaceJointFrictionError> {
    let result = a * b;
    if !result.is_finite() || (a != 0.0 && b != 0.0 && result == 0.0) {
        Err(InterfaceJointFrictionError::Invalid(
            "unrepresentable physical effort or speed",
        ))
    } else {
        Ok(result)
    }
}

/// Resolve constant Coulomb interface cards and advance their actual native flow.
///
/// The original law owns its static/kinetic coefficients and domain. Non-Coulomb
/// rungs refuse; their zero-speed values cannot replace their speed-dependent
/// laws. Inputs and source cards are borrowed, so a late speed-domain refusal,
/// ambiguous mechanical continuation, numerical failure or cancellation publishes
/// no partial momentum or source state.
///
/// Temperature, total normal load, nominal area and sleeve radius are held during
/// this phase. They are supplied by the caller's contact/material model, never
/// reconstructed from a net bilateral reaction. Distributed nonuniform pressure
/// requires its own applicability justification; N/A here is explicitly nominal.
/// These cards retain their original caller-declared/estimated/synthetic authority.
/// They are not substituted for fs-matdb source admission or empirical validation.
/// Poses, wear, heat partition and contact geometry evolution remain with their
/// existing owners. In particular, the hinge angular length is only conditioning.
#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
pub fn advance_interface_joint_friction<'a>(
    states: &[RigidBodyState],
    properties: &[MassProperties],
    joints: &[ArticulatedJoint],
    bindings: &[JointDryInterface<'a>],
    duration_s: f64,
    controls: JointFrictionFlowControls,
    mut cancelled: impl FnMut() -> bool,
) -> Result<InterfaceJointFrictionReport<'a>, InterfaceJointFrictionError> {
    checkpoint(&mut cancelled)?;
    if bindings.is_empty() || bindings.len() > MAX_JOINT_FRICTION_PORTS {
        return Err(InterfaceJointFrictionError::Invalid(
            "bounded nonempty dry joint interfaces required",
        ));
    }
    let mut ports = [JointFrictionFlowPort::Prismatic {
        joint_index: 0,
        capacity_n: 0.0,
        drive_n: 0.0,
    }; MAX_JOINT_FRICTION_PORTS];
    let mut pressures = [0.0; MAX_JOINT_FRICTION_PORTS];
    let mut body_order = [[0; 2]; MAX_JOINT_FRICTION_PORTS];
    for (i, binding) in bindings.iter().enumerate() {
        checkpoint(&mut cancelled)?;
        if bindings[..i]
            .iter()
            .any(|old| old.joint_index == binding.joint_index)
            || !binding.normal_force_n.is_finite()
            || binding.normal_force_n < 0.0
            || !binding.nominal_area_m2.is_finite()
            || binding.nominal_area_m2 <= 0.0
        {
            return Err(InterfaceJointFrictionError::Invalid(
                "distinct joints and physical normal load/area required",
            ));
        }
        let constitutive = |error| InterfaceJointFrictionError::Constitutive {
            joint_index: binding.joint_index,
            error,
        };
        let law = binding.card.friction_law();
        if !matches!(law, FrictionLaw::Coulomb { .. }) {
            return Err(InterfaceJointFrictionError::UnsupportedLaw {
                joint_index: binding.joint_index,
            });
        }
        let pressure = binding.normal_force_n / binding.nominal_area_m2;
        if !pressure.is_finite() || (binding.normal_force_n != 0.0 && pressure == 0.0) {
            return Err(InterfaceJointFrictionError::Invalid(
                "representable nominal pressure required",
            ));
        }
        // Preflight the held coordinates at the card's own valid speed bound.
        // This is no physical query receipt; actual speed extrema are checked below.
        binding
            .card
            .validate_state(
                DryFrictionState::new(
                    binding.temperature_k,
                    pressure,
                    binding.card.applicability().slip_speed_mps().minimum(),
                )
                .map_err(constitutive)?,
            )
            .map_err(constitutive)?;
        let static_force = multiply(
            law.static_coefficient().map_err(constitutive)?,
            binding.normal_force_n,
        )?;
        let kinetic_force = multiply(
            law.zero_speed_kinetic_coefficient().map_err(constitutive)?,
            binding.normal_force_n,
        )?;
        ports[i] = match (binding.geometry, joints.get(binding.joint_index)) {
            (
                JointDrySlipGeometry::Prismatic { drive_n },
                Some(ArticulatedJoint::Prismatic(joint)),
            ) => {
                body_order[i] = joint.bodies();
                JointFrictionFlowPort::PrismaticStiction {
                    joint_index: binding.joint_index,
                    static_capacity_n: static_force,
                    kinetic_capacity_n: kinetic_force,
                    drive_n,
                }
            }
            (
                JointDrySlipGeometry::CylindricalSleeve { radius_m, drive_nm },
                Some(ArticulatedJoint::Revolute(joint)),
            ) => {
                if !radius_m.is_finite() || radius_m <= 0.0 {
                    return Err(InterfaceJointFrictionError::Invalid(
                        "positive physical sleeve radius required",
                    ));
                }
                body_order[i] = joint.bodies();
                JointFrictionFlowPort::RevoluteStiction {
                    joint_index: binding.joint_index,
                    static_capacity_nm: multiply(static_force, radius_m)?,
                    kinetic_capacity_nm: multiply(kinetic_force, radius_m)?,
                    drive_nm,
                }
            }
            _ => {
                return Err(InterfaceJointFrictionError::Invalid(
                    "dry slip geometry must match the original joint",
                ));
            }
        };
        pressures[i] = pressure;
    }
    let mechanical = advance_joint_friction(
        states,
        properties,
        joints,
        &ports[..bindings.len()],
        duration_s,
        controls,
        &mut cancelled,
    )
    .map_err(InterfaceJointFrictionError::Mechanical)?;
    let mut interfaces = Vec::new();
    interfaces
        .try_reserve_exact(bindings.len())
        .map_err(|_| InterfaceJointFrictionError::Allocation)?;
    for port in &mechanical.ports {
        checkpoint(&mut cancelled)?;
        let i = bindings
            .iter()
            .position(|b| b.joint_index == port.joint_index)
            .ok_or(InterfaceJointFrictionError::Invalid(
                "missing original interface binding",
            ))?;
        let binding = bindings[i];
        let (length, radius, endpoint_rates) = match port.motion {
            JointFrictionFlowMotion::Revolute {
                angular_length_m,
                rate_before_rad_s,
                rate_after_rad_s,
                ..
            } => match binding.geometry {
                JointDrySlipGeometry::CylindricalSleeve { radius_m, .. } => (
                    angular_length_m,
                    radius_m,
                    [rate_before_rad_s.abs(), rate_after_rad_s.abs()],
                ),
                _ => {
                    return Err(InterfaceJointFrictionError::Invalid(
                        "changed original joint kind",
                    ));
                }
            },
            JointFrictionFlowMotion::Prismatic {
                velocity_before_m_s,
                velocity_after_m_s,
                ..
            } => (
                1.0,
                1.0,
                [velocity_before_m_s.abs(), velocity_after_m_s.abs()],
            ),
        };
        let mut range =
            [DryFrictionState::new(binding.temperature_k, pressures[i], 0.0).map_err(|error| {
                InterfaceJointFrictionError::Constitutive {
                    joint_index: binding.joint_index,
                    error,
                }
            })?; 2];
        for (index, (state, scaled)) in range
            .iter_mut()
            .zip([port.minimum_scaled_speed_m_s, port.maximum_scaled_speed_m_s])
            .enumerate()
        {
            let mut rate = scaled / length;
            if !rate.is_finite() || rate < 0.0 || (scaled != 0.0 && rate == 0.0) {
                return Err(InterfaceJointFrictionError::Invalid(
                    "representable original joint rate required",
                ));
            }
            // Keep unscaled native endpoints too: L multiplication/division
            // can round across a strict source-domain boundary.
            rate = if index == 0 {
                rate.min(endpoint_rates[0]).min(endpoint_rates[1])
            } else {
                rate.max(endpoint_rates[0]).max(endpoint_rates[1])
            };
            let error = |error| InterfaceJointFrictionError::Constitutive {
                joint_index: binding.joint_index,
                error,
            };
            *state =
                DryFrictionState::new(binding.temperature_k, pressures[i], multiply(rate, radius)?)
                    .map_err(error)?;
            binding.card.validate_state(*state).map_err(error)?;
        }
        interfaces.push(JointDryInterfaceReceipt {
            binding,
            ordered_bodies: body_order[i],
            nominal_pressure_pa: pressures[i],
            source: ports[i],
            admitted_state_range: range,
        });
    }
    checkpoint(&mut cancelled)?;
    Ok(InterfaceJointFrictionReport {
        mechanical,
        interfaces,
    })
}
