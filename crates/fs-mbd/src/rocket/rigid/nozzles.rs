//! Physical point-nozzle arrays for the native full-tensor rocket.
//!
//! Each off-center group has equal mass rates at opposite material arms. Central
//! nozzles are also admitted. This makes the exported mass first moment exactly
//! zero by construction, as required by the existing fixed-COM trajectory model.
//! Exhaust velocities are relative to the LOCAL rotating nozzle, not the COM.
//!
//! For total rate q, relative momentum pi, relative angular momentum j and
//! co-rotation tensor K, preparation derives
//! pi = sum(q_i c_i)/q, j = sum(q_i r_i cross c_i)/q,
//! K = sum(q_i (|r_i|^2 I - r_i r_i^T))/q.
//! The reaction is -q*pi, while ALL relative jet torque is exported as q*j.
//! The native thrust arm is therefore zero: adding a second arm torque would
//! double count the same particle momentum. Pure couples remain representable.
//!
//! The native owner retains its original mass/inertia law, burnout, temporal
//! resolution and engine-work assumptions. Balanced boundary flux does not
//! determine reservoir inertia, internal relative-flow momentum, slosh, nozzle
//! pressure thrust or combustion energy. Relative kinetic power below is a
//! stream diagnostic, not an extra heat debit or a complete energy balance.

use super::super::{RigidRocketError, RigidRocketSolve, ScalarSum, VectorSum, poll};
use super::{SpecificInertia, TensorRocketForcing, TensorRocketState, TensorRocketStep};
use crate::Vec3;
use fs_ga::Mat3;

/// Maximum number of central or opposed groups in one bounded preparation.
pub const MAX_NOZZLE_GROUPS: usize = 16_384;

/// A physical exhaust group with exactly zero mass-flow first moment.
///
/// Every velocity points in the actual particle travel direction relative to
/// its rotating material nozzle. All arms use the fixed material body axes.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum BalancedNozzleGroup {
    /// A single physical nozzle at the common center of mass.
    Central {
        /// Outgoing mass rate [kg/s], nonnegative.
        mass_flow_kg_s: f64,
        /// Physical exhaust velocity relative to the COM nozzle [m/s].
        exhaust_velocity_body_m_s: Vec3,
    },
    /// Two physical point nozzles at +arm and -arm with equal outgoing rates.
    OpposedPair {
        /// Rate of EACH nozzle [kg/s]; the group consumes twice this rate.
        mass_flow_per_nozzle_kg_s: f64,
        /// Material arm of the positive nozzle from the common COM [m].
        arm_body_m: Vec3,
        /// Physical relative exhaust velocity at +arm [m/s].
        positive_exhaust_velocity_body_m_s: Vec3,
        /// Physical relative exhaust velocity at -arm [m/s].
        negative_exhaust_velocity_body_m_s: Vec3,
    },
}

impl BalancedNozzleGroup {
    // A central group's second entry is an inactive point. No heap or retained
    // per-port storage is needed, and paired first moments are never estimated.
    fn ports(self) -> [(f64, Vec3, Vec3); 2] {
        match self {
            Self::Central {
                mass_flow_kg_s,
                exhaust_velocity_body_m_s,
            } => [
                (mass_flow_kg_s, Vec3::ZERO, exhaust_velocity_body_m_s),
                (0.0, Vec3::ZERO, Vec3::ZERO),
            ],
            Self::OpposedPair {
                mass_flow_per_nozzle_kg_s,
                arm_body_m,
                positive_exhaust_velocity_body_m_s,
                negative_exhaust_velocity_body_m_s,
            } => [
                (
                    mass_flow_per_nozzle_kg_s,
                    arm_body_m,
                    positive_exhaust_velocity_body_m_s,
                ),
                (
                    mass_flow_per_nozzle_kg_s,
                    arm_body_m.scale(-1.0),
                    negative_exhaust_velocity_body_m_s,
                ),
            ],
        }
    }
}

/// Constant-size physical stream moments, reusable throughout a held valve state.
///
/// Preparation traverses the groups twice with compensated reductions. Native
/// stepping and instantaneous flux queries then take constant storage/work,
// independent of nozzle count. Changing rates, gimbals or geometry requires a
/// new preparation; no state-dependent nozzle law is silently frozen.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PreparedNozzleArray {
    mass_flow_kg_s: f64,
    mean_exhaust_velocity_body_m_s: Vec3,
    relative_specific_angular_momentum_body_m2_s: Vec3,
    co_rotating_specific_inertia_m2: SpecificInertia,
    relative_kinetic_power_w: f64,
}

/// Instantaneous physical export at the current attitude and angular velocity.
///
/// These are transport rates, not interval impulses or an extra applied wrench.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct NozzleFlux {
    /// Total outgoing mass [kg/s], zero after native burnout.
    pub mass_flow_kg_s: f64,
    /// Outgoing world linear momentum [N], including COM velocity advection.
    pub linear_momentum_rate_world_n: Vec3,
    /// Outgoing angular momentum about the current COM [N m].
    /// Includes both orbital co-rotation and relative tangential jet transport.
    pub spin_angular_momentum_rate_world_n_m: Vec3,
    /// Outgoing angular momentum about the supplied world point [N m].
    /// Includes the COM orbital term once.
    pub angular_momentum_rate_about_origin_world_n_m: Vec3,
    /// Sum of one half mass rate times squared LOCAL relative exhaust speed [W].
    /// Excludes thermal/chemical energy and world-frame kinetic advection.
    pub relative_kinetic_power_w: f64,
}

impl NozzleFlux {
    /// Closed valves or complete native burnout.
    pub const ZERO: Self = Self {
        mass_flow_kg_s: 0.0,
        linear_momentum_rate_world_n: Vec3::ZERO,
        spin_angular_momentum_rate_world_n_m: Vec3::ZERO,
        angular_momentum_rate_about_origin_world_n_m: Vec3::ZERO,
        relative_kinetic_power_w: 0.0,
    };
}

impl PreparedNozzleArray {
    /// No mass transport, force, relative couple or co-rotation flux.
    pub const ZERO: Self = Self {
        mass_flow_kg_s: 0.0,
        mean_exhaust_velocity_body_m_s: Vec3::ZERO,
        relative_specific_angular_momentum_body_m2_s: Vec3::ZERO,
        co_rotating_specific_inertia_m2: SpecificInertia::ZERO,
        relative_kinetic_power_w: 0.0,
    };

    /// Prepare physical particle momentum from the same rates and port geometry.
    ///
    /// Empty/all-closed groups yield ZERO. Inactive inputs are still validated.
    /// Cancellation is polled initially, every 64 groups in BOTH passes, and
    /// before publishing. Non-finite or unrepresentable moments are refused.
    /// Rank-deficient point/line exit inertias are admitted without adding
    /// artificial inertia or projecting the tensor onto its diagonal.
    pub fn new(
        groups: &[BalancedNozzleGroup],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, RigidRocketError> {
        poll(&mut cancelled)?;
        if groups.len() > MAX_NOZZLE_GROUPS {
            return Err(RigidRocketError::Invalid("too many balanced nozzle groups"));
        }
        let mut mass = ScalarSum::default();
        for (index, group) in groups.iter().enumerate() {
            if index % 64 == 0 {
                poll(&mut cancelled)?;
            }
            for (rate, arm, velocity) in group.ports() {
                if !rate.is_finite() || rate < 0.0 || !arm.is_finite() || !velocity.is_finite() {
                    return Err(RigidRocketError::Invalid(
                        "finite nozzle geometry, velocity and nonnegative mass rate required",
                    ));
                }
                mass.add(rate);
            }
        }
        if !mass.value.is_finite() {
            return Err(RigidRocketError::Invalid("nozzle mass rate overflow"));
        }
        if mass.value == 0.0 {
            poll(&mut cancelled)?;
            return Ok(Self::ZERO);
        }

        let mut linear = VectorSum::default();
        let mut angular = VectorSum::default();
        let mut power_sum = ScalarSum::default();
        // xx, yy, zz, xy, xz, yz: symmetry is constructed, never repaired.
        let mut tensor: [ScalarSum; 6] = core::array::from_fn(|_| ScalarSum::default());
        for (index, group) in groups.iter().enumerate() {
            if index % 64 == 0 {
                poll(&mut cancelled)?;
            }
            for (rate, arm, velocity) in group.ports() {
                if rate == 0.0 {
                    continue;
                }
                let fraction = rate / mass.value;
                if fraction == 0.0 {
                    return Err(RigidRocketError::Invalid("nozzle rate fraction underflow"));
                }
                let root = fraction.sqrt();
                let r = arm.scale(root);
                let c = velocity.scale(root);
                linear.add(velocity.scale(fraction));
                angular.add(r.cross(c));
                let moments = [
                    r.y.mul_add(r.y, r.z * r.z),
                    r.x.mul_add(r.x, r.z * r.z),
                    r.x.mul_add(r.x, r.y * r.y),
                    -r.x * r.y,
                    -r.x * r.z,
                    -r.y * r.z,
                ];
                for (sum, value) in tensor.iter_mut().zip(moments) {
                    sum.add(value);
                }
                // Halve before squaring to avoid an unnecessary overflow at
                // the top of the representable kinetic-energy range.
                let energy_velocity = velocity.scale(rate.sqrt());
                power_sum.add((0.5 * energy_velocity.x).mul_add(
                    energy_velocity.x,
                    (0.5 * energy_velocity.y).mul_add(
                        energy_velocity.y,
                        (0.5 * energy_velocity.z) * energy_velocity.z,
                    ),
                ));
            }
        }
        let [xx, yy, zz, xy, xz, yz] = tensor.map(|sum| sum.value);
        let moment = Mat3 {
            m: [xx, xy, xz, xy, yy, yz, xz, yz, zz],
        };
        let power = power_sum.value;
        if moment.m.iter().any(|x| !x.is_finite())
            || !linear.value.is_finite()
            || !angular.value.is_finite()
            || !linear.value.scale(mass.value).is_finite()
            || !angular.value.scale(mass.value).is_finite()
            || !power.is_finite()
        {
            return Err(RigidRocketError::Invalid(
                "nozzle momentum or power overflow",
            ));
        }
        poll(&mut cancelled)?;
        Ok(Self {
            mass_flow_kg_s: mass.value,
            mean_exhaust_velocity_body_m_s: linear.value,
            relative_specific_angular_momentum_body_m2_s: angular.value,
            // Every summand is a nonnegative point-mass inertia. This
            // descendant of tensor retains its prepared-type invariant even
            // when the legitimate exit tensor has a zero principal moment.
            co_rotating_specific_inertia_m2: SpecificInertia {
                tensor_body_m2: moment,
            },
            relative_kinetic_power_w: power,
        })
    }

    /// Powered total outgoing mass rate [kg/s].
    #[must_use]
    pub const fn mass_flow_kg_s(&self) -> f64 {
        self.mass_flow_kg_s
    }

    /// Relative-particle reaction force held in body axes [N].
    #[must_use]
    pub fn reaction_force_body_n(&self) -> Vec3 {
        self.mean_exhaust_velocity_body_m_s
            .scale(-self.mass_flow_kg_s)
    }

    /// Relative-particle reaction couple held in body axes [N m].
    /// Co-rotating angular transport is additional and depends on angular velocity.
    #[must_use]
    pub fn reaction_couple_body_n_m(&self) -> Vec3 {
        self.relative_specific_angular_momentum_body_m2_s
            .scale(-self.mass_flow_kg_s)
    }

    /// Physical co-rotating outlet inertia per exported mass [m²].
    #[must_use]
    pub const fn co_rotating_specific_inertia_m2(&self) -> SpecificInertia {
        self.co_rotating_specific_inertia_m2
    }

    /// Powered LOCAL relative kinetic flux [W]; no thermal/chemical debit is made.
    #[must_use]
    pub const fn relative_kinetic_power_w(&self) -> f64 {
        self.relative_kinetic_power_w
    }

    /// Compile directly into the native full-tensor rocket law.
    ///
    /// Feed this to TensorRocketState::advance_spatial for additional prepared
    /// body/world forces and couples. Propulsion must not be applied a second
    /// time as an external load. Pressure thrust is NOT inferred from a particle
    /// velocity; a separately owned pressure resultant may be added explicitly.
    #[must_use]
    pub fn forcing(&self, gravity_world_m_s2: Vec3) -> TensorRocketForcing {
        TensorRocketForcing {
            mass_flow_kg_s: self.mass_flow_kg_s,
            effective_exhaust_velocity_body_m_s: self.mean_exhaust_velocity_body_m_s.scale(-1.0),
            thrust_application_arm_body_m: Vec3::ZERO,
            co_rotating_exhaust_specific_inertia_m2: self.co_rotating_specific_inertia_m2,
            exhaust_swirl_specific_angular_momentum_body_m2_s: self
                .relative_specific_angular_momentum_body_m2_s,
            gravity_world_m_s2,
            external_force_world_n: Vec3::ZERO,
            external_torque_body_n_m: Vec3::ZERO,
        }
    }

    /// Evaluate absolute linear and angular particle transport about a world point.
    ///
    /// The zero first moment cancels the COM-velocity spin term and net
    /// co-rotation linear transport. The remaining terms are q*(V+R*pi),
    /// R*q*(K*omega+j), and (X-origin) cross linear + spin. All valves close
    /// together at native propellant exhaustion; no dry mass is expelled.
    pub fn flux_at(
        &self,
        state: &TensorRocketState,
        origin_world_m: Vec3,
    ) -> Result<NozzleFlux, RigidRocketError> {
        if !origin_world_m.is_finite() {
            return Err(RigidRocketError::Invalid(
                "finite nozzle flux origin required",
            ));
        }
        if self.mass_flow_kg_s == 0.0 || state.translation().propellant_mass_kg == 0.0 {
            return Ok(NozzleFlux::ZERO);
        }
        let omega = state
            .mass_properties()?
            .angular_velocity_body_checked(state.angular_momentum_body_kg_m2_s())?;
        let k = self.co_rotating_specific_inertia_m2.tensor_body_m2().m;
        let corotation = Vec3::new(
            k[0].mul_add(omega.x, k[1].mul_add(omega.y, k[2] * omega.z)),
            k[3].mul_add(omega.x, k[4].mul_add(omega.y, k[5] * omega.z)),
            k[6].mul_add(omega.x, k[7].mul_add(omega.y, k[8] * omega.z)),
        );
        let linear = state
            .translation()
            .velocity_m_s
            .add(
                state
                    .orientation()
                    .rotate_body_to_world(self.mean_exhaust_velocity_body_m_s),
            )
            .scale(self.mass_flow_kg_s);
        let spin = state.orientation().rotate_body_to_world(
            corotation
                .add(self.relative_specific_angular_momentum_body_m2_s)
                .scale(self.mass_flow_kg_s),
        );
        let orbital = state
            .translation()
            .position_m
            .sub(origin_world_m)
            .cross(linear);
        let angular = orbital.add(spin);
        if !linear.is_finite() || !spin.is_finite() || !angular.is_finite() {
            return Err(RigidRocketError::Invalid(
                "absolute nozzle momentum flux overflow",
            ));
        }
        Ok(NozzleFlux {
            mass_flow_kg_s: self.mass_flow_kg_s,
            linear_momentum_rate_world_n: linear,
            spin_angular_momentum_rate_world_n_m: spin,
            angular_momentum_rate_about_origin_world_n_m: angular,
            relative_kinetic_power_w: self.relative_kinetic_power_w,
        })
    }
}

impl TensorRocketState {
    /// Advance simultaneous physical nozzles through the original bounded owner.
    ///
    /// Rates/geometry remain held during this interval. This uses native exact
    /// burnout segmentation and second-order attitude/outflow coupling, not
    /// sequential per-nozzle kicks. Use the prepared forcing with advance_spatial
    /// when additional loads are present.
    pub fn advance_nozzles(
        &self,
        nozzles: &PreparedNozzleArray,
        gravity_world_m_s2: Vec3,
        duration_s: f64,
        controls: RigidRocketSolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<TensorRocketStep, RigidRocketError> {
        self.advance(
            nozzles.forcing(gravity_world_m_s2),
            duration_s,
            controls,
            cancelled,
        )
    }
}

#[cfg(test)]
#[path = "nozzles/tests.rs"]
mod tests;
