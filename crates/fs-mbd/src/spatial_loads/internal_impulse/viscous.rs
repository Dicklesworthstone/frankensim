//! Exact passive spatial-wrench relaxation between two finite rigid bodies.
//!
//! The damping matrix acts on the relative twist at attachment B, expressed
//! in A's material axes. A's velocity is transported to B before subtraction:
//! s = [v_B - v_A - omega_A cross (x_B-x_A), omega_B-omega_A].
//! Thus one common rigid motion produces no viscous effort, even when the
//! attachment points differ. The opposite wrench includes the SAME native
//! moment transport; shear forces cannot create total angular momentum.
//!
//! Preparation uses the original passive six-coordinate form and native
//! Jacobi owner. At fixed geometry, G = B M^-1 B^T and s' = -G C s.
//! A symmetric eigendecomposition on the original damping range integrates this linear
//! source exactly, including stiff and semidefinite damping. Relaxation
//! applies one simultaneous native internal wrench, retaining original
//! momentum/work admission and independently integrated nonnegative heat.
//!
//! This is a composable source substep: neither pose moves, no constraint or
//! spring is imposed, and no complete changing-geometry trajectory is claimed.
//! Reprepare after either pose, attachment, inertia, coefficient or interval
//! changes. A symmetric composition with native free motion still requires
//! temporal refinement for changing geometry; no thermal destination is inferred.

/// Sparse symmetric compositions of the original finite-body connector source.
#[path = "viscous/network.rs"]
pub mod network;

/// Simultaneous held internal actuation and passive spatial relaxation.
#[path = "viscous/driven.rs"]
pub mod driven;

use super::{
    InternalWrenchImpulse, InternalWrenchImpulseControls, InternalWrenchImpulseError,
    InternalWrenchImpulseReceipt, apply_internal_wrench_impulse, dd, dot, finite, finite_vector,
    poll, prepare, rotate_wide, rounded, subtract, validate_controls,
};
pub use crate::articulated::motion::spatial_loads::bushing::coupled::{
    BushingCoordinateScales, BushingMatrixError, PassiveBushingMatrix,
};
use crate::{DynamicsError, MassProperties, Pose, RigidBodyState, Vec3};
use fs_la::eigen::jacobi_eigh;
use fs_math::{dd::Dd, det};

type Six = [Dd; 6];
type Matrix = [[f64; 6]; 6];
const SPECTRAL_TOLERANCE: f64 = 8192.0 * f64::EPSILON;

/// Native physical balance limits plus accuracy of the accepted relative twist.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ViscousConnectorControls {
    /// Original per-body and total momentum/work admission. Its work tolerance
    /// also bounds kinetic change plus heat and connector work plus heat.
    pub impulse: InternalWrenchImpulseControls,
    /// Maximum componentwise endpoint relative-translation-rate defect [m/s].
    pub linear_speed_tolerance_m_s: f64,
    /// Maximum componentwise endpoint relative-spin-rate defect [rad/s].
    pub angular_speed_tolerance_rad_s: f64,
}

/// Complete accepted momentum exchange and separate physical heat.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ViscousConnectorReceipt {
    /// Original native impulse receipt; both original poses remain unchanged.
    pub impulse: InternalWrenchImpulseReceipt,
    /// Transported relative linear velocity and angular velocity at B,
    /// in A's material axes, ordered [m/s; 3] then [rad/s; 3].
    pub relative_twist_before: [f64; 6],
    /// The same twist independently recovered from the actual native outputs.
    pub relative_twist_after: [f64; 6],
    /// Actual endpoint twist minus the analytic fixed-geometry relaxation.
    pub relative_twist_residual: [f64; 6],
    /// Integral of s^T C s over the source interval [J], computed from the
    /// original modes; never defined by a measured kinetic-energy loss.
    pub dissipated_energy_j: f64,
    /// Actual native kinetic change plus the independently integrated heat [J].
    pub kinetic_heat_residual_j: f64,
    /// Original native connector work plus the independently integrated heat [J].
    pub connector_heat_residual_j: f64,
}

/// Atomic refusal: neither input body, damping law nor preparation is changed.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum ViscousConnectorError {
    /// A duration, binding, control or numerical source is inadmissible.
    Invalid(&'static str),
    /// Original rigid-body impulse or physical-account refusal.
    Impulse(InternalWrenchImpulseError),
    /// Original passive-form evaluation refused.
    Damping(BushingMatrixError),
    /// Actual endpoint motion or heat/work differs from the declared source.
    Unresolved {
        /// Largest linear-rate component defect [m/s].
        linear_speed_residual_m_s: f64,
        /// Largest angular-rate component defect [rad/s].
        angular_speed_residual_rad_s: f64,
        /// Actual kinetic change plus heat [J].
        kinetic_heat_residual_j: f64,
        /// Original connector work plus heat [J].
        connector_heat_residual_j: f64,
    },
}
impl From<InternalWrenchImpulseError> for ViscousConnectorError {
    fn from(value: InternalWrenchImpulseError) -> Self {
        Self::Impulse(value)
    }
}
impl From<DynamicsError> for ViscousConnectorError {
    fn from(value: DynamicsError) -> Self {
        Self::Impulse(value.into())
    }
}
impl From<BushingMatrixError> for ViscousConnectorError {
    fn from(value: BushingMatrixError) -> Self {
        Self::Damping(value)
    }
}
impl core::fmt::Display for ViscousConnectorError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "spatial viscous connector refused: {self:?}")
    }
}
impl std::error::Error for ViscousConnectorError {}

/// Fixed-geometry six-coordinate passive relaxation on original native bodies.
///
/// C is a supplied symmetric positive-semidefinite SI damping matrix in
/// [translation, rotation] order, expressed in body A's material axes. Its
/// original admission and physical coordinate scales remain authoritative.
/// These coefficients act on spatial twist, not on a principal rotation-vector
/// derivative; the shared matrix type supplies numerical passivity admission.
/// All six coordinates relax simultaneously, including shear, spin, and
/// translation/rotation cross-coupling. Both finite masses and full inertias
/// contribute to mobility; neither body is an infinite prescribed bath.
///
/// Preparation has bounded six-coordinate work and one native small eigenproblem.
/// The retained object has fixed-size storage. Each subsequent relaxation is
/// allocation-free, performs no factorization/eigenproblem/exponential, and
/// calls the original native internal-impulse owner once.
#[derive(Clone, Debug)]
pub struct PreparedViscousConnector {
    poses: [Pose; 2],
    properties: [MassProperties; 2],
    arms_body_m: [Vec3; 2],
    axes_world: [Vec3; 3],
    force_columns_world: [Vec3; 6],
    angular_columns_body: [[Vec3; 2]; 6],
    damping: PassiveBushingMatrix,
    duration_s: f64,
    maximum_decay_rate_per_s: f64,
    relative_modes: Matrix,
    source_modes: Matrix,
    integrated_decay: [f64; 6],
    heat_root: [f64; 6],
}

impl PreparedViscousConnector {
    /// Bind original geometry, complete inertia, passive damping and interval.
    ///
    /// Zero duration or a zero matrix is an admitted identity source. Nonzero
    /// intervals must be positive. Coincident attachments are allowed: the
    /// constitutive matrix supplies the physical directions, not a guessed axis.
    /// Preparation refuses unresolved active mobility or spectral
    /// residuals; it never clips negative rates or adds a damping floor.
    pub fn new(
        poses: [Pose; 2],
        properties: [MassProperties; 2],
        arms_body_m: [Vec3; 2],
        damping: PassiveBushingMatrix,
        duration_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, ViscousConnectorError> {
        poll(&mut cancelled)?;
        if !duration_s.is_finite() || duration_s < 0.0 {
            return Err(ViscousConnectorError::Invalid(
                "finite nonnegative duration",
            ));
        }
        let states = [
            RigidBodyState::new(poses[0], Vec3::ZERO, Vec3::ZERO)?,
            RigidBodyState::new(poses[1], Vec3::ZERO, Vec3::ZERO)?,
        ];
        let axes_world = [
            Vec3::new(1.0, 0.0, 0.0),
            Vec3::new(0.0, 1.0, 0.0),
            Vec3::new(0.0, 0.0, 1.0),
        ]
        .map(|axis| poses[0].orientation().rotate_body_to_world(axis));
        let mut force_columns_world = [Vec3::ZERO; 6];
        let mut angular_columns_body = [[Vec3::ZERO; 2]; 6];
        for column in 0..6 {
            poll(&mut cancelled)?;
            let force = if column < 3 {
                axes_world[column]
            } else {
                Vec3::ZERO
            };
            let couple = if column < 3 {
                Vec3::ZERO
            } else {
                axes_world[column - 3]
            };
            let native = prepare(
                states,
                properties,
                InternalWrenchImpulse {
                    arms_body_m,
                    impulse_on_second_world_ns: force,
                    couple_on_second_world_nms: couple,
                },
            )?;
            force_columns_world[column] = force;
            for body in 0..2 {
                let moment = poses[body]
                    .orientation()
                    .rotate_world_to_body(rounded(native.central[body])?);
                finite_vector(moment)?;
                angular_columns_body[column][body] = moment;
            }
        }
        let mut owner = Self {
            poses,
            properties,
            arms_body_m,
            axes_world,
            force_columns_world,
            angular_columns_body,
            damping,
            duration_s,
            maximum_decay_rate_per_s: 0.0,
            relative_modes: [[0.0; 6]; 6],
            source_modes: [[0.0; 6]; 6],
            integrated_decay: [0.0; 6],
            heat_root: [0.0; 6],
        };
        if duration_s != 0.0
            && owner
                .damping
                .native_coefficients()
                .iter()
                .flatten()
                .any(|x| *x != 0.0)
        {
            owner.prepare_modes(&mut cancelled)?;
        }
        poll(&mut cancelled)?;
        Ok(owner)
    }

    /// Original interval [s], held during this fixed-geometry source.
    pub const fn duration_s(&self) -> f64 {
        self.duration_s
    }

    /// Largest active decay rate retained by the original prepared kinetic
    /// eigensystem [1/s]. This is zero when a zero duration or zero damping law
    /// skipped spectral preparation; zero duration does not certify a zero
    /// physical decay rate for a nonzero law.
    pub const fn maximum_decay_rate_per_s(&self) -> f64 {
        self.maximum_decay_rate_per_s
    }

    /// Original supplied passive form; no matrix projection or fitted damping.
    pub const fn damping(&self) -> &PassiveBushingMatrix {
        &self.damping
    }

    /// Relax arbitrary native momenta at the exact prepared poses.
    ///
    /// Actual source impulses, momentum transfer and work retain the original
    /// native gates. Relative-twist and analytic-heat checks additionally admit
    /// the accepted endpoint; opposite swallowed impulses cannot appear as
    /// damping. Invalid input, stale geometry or cancellation publishes nothing.
    pub fn relax(
        &self,
        states: [RigidBodyState; 2],
        controls: ViscousConnectorControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ViscousConnectorReceipt, ViscousConnectorError> {
        poll(&mut cancelled)?;
        validate_controls(controls.impulse)?;
        if [
            controls.linear_speed_tolerance_m_s,
            controls.angular_speed_tolerance_rad_s,
        ]
        .into_iter()
        .any(|x| !x.is_finite() || x < 0.0)
        {
            return Err(ViscousConnectorError::Invalid(
                "finite nonnegative twist tolerances",
            ));
        }
        if states.map(|state| state.pose()) != self.poses {
            return Err(ViscousConnectorError::Invalid(
                "changed prepared body poses",
            ));
        }
        let before = self.relative_twist(states)?;
        let mut after_analytic = before;
        let mut source = [Dd::ZERO; 6];
        let mut heat = Dd::ZERO;
        for mode in 0..6 {
            poll(&mut cancelled)?;
            let y = product(column(self.source_modes, mode), before);
            let change = -Dd::from_f64(self.integrated_decay[mode]) * y;
            for i in 0..6 {
                source[i] = source[i] + Dd::from_f64(self.source_modes[i][mode]) * change;
                after_analytic[i] =
                    after_analytic[i] + Dd::from_f64(self.relative_modes[i][mode]) * change;
            }
            let dissipated = Dd::from_f64(self.heat_root[mode]) * y;
            heat = heat + dissipated * dissipated;
        }
        let impulse = apply_internal_wrench_impulse(
            states,
            self.properties,
            InternalWrenchImpulse {
                arms_body_m: self.arms_body_m,
                impulse_on_second_world_ns: rounded(rotate_wide(
                    self.axes_world,
                    [source[0], source[1], source[2]],
                ))?,
                couple_on_second_world_nms: rounded(rotate_wide(
                    self.axes_world,
                    [source[3], source[4], source[5]],
                ))?,
            },
            controls.impulse,
            &mut cancelled,
        )?;
        poll(&mut cancelled)?;
        let after = self.relative_twist(impulse.bodies.map(|body| body.state_after))?;
        let residual = rounded_six(core::array::from_fn(|i| after[i] - after_analytic[i]))?;
        let linear_residual = residual[..3].iter().map(|x| x.abs()).fold(0.0, f64::max);
        let angular_residual = residual[3..].iter().map(|x| x.abs()).fold(0.0, f64::max);
        let dissipated_energy_j = finite(heat.to_f64())?;
        let kinetic_heat_residual_j =
            finite((Dd::from_f64(impulse.kinetic_energy_change_j) + heat).to_f64())?;
        let connector_heat_residual_j =
            finite((Dd::from_f64(impulse.connector_work_j) + heat).to_f64())?;
        if dissipated_energy_j < 0.0
            || linear_residual > controls.linear_speed_tolerance_m_s
            || angular_residual > controls.angular_speed_tolerance_rad_s
            || kinetic_heat_residual_j.abs() > controls.impulse.work_tolerance_j
            || connector_heat_residual_j.abs() > controls.impulse.work_tolerance_j
        {
            return Err(ViscousConnectorError::Unresolved {
                linear_speed_residual_m_s: linear_residual,
                angular_speed_residual_rad_s: angular_residual,
                kinetic_heat_residual_j,
                connector_heat_residual_j,
            });
        }
        let result = ViscousConnectorReceipt {
            impulse,
            relative_twist_before: rounded_six(before)?,
            relative_twist_after: rounded_six(after)?,
            relative_twist_residual: residual,
            dissipated_energy_j,
            kinetic_heat_residual_j,
            connector_heat_residual_j,
        };
        poll(&mut cancelled)?;
        Ok(result)
    }

    fn relative_twist(&self, states: [RigidBodyState; 2]) -> Result<Six, ViscousConnectorError> {
        let velocity = [
            states[0].center_of_mass_velocity_world(self.properties[0])?,
            states[1].center_of_mass_velocity_world(self.properties[1])?,
        ];
        let omega = [
            self.properties[0].angular_velocity_body_checked(states[0].angular_momentum_body())?,
            self.properties[1].angular_velocity_body_checked(states[1].angular_momentum_body())?,
        ];
        // Cancel bulk translation before adding either spin contribution.
        // The native wrench columns carry A's transport to B, so their dual
        // rate is objective even for separated, rotating attachments.
        let translation = subtract(dd(velocity[1]), dd(velocity[0]));
        Ok(core::array::from_fn(|i| {
            dot(dd(self.force_columns_world[i]), translation)
                + dot(dd(self.angular_columns_body[i][0]), dd(omega[0]))
                + dot(dd(self.angular_columns_body[i][1]), dd(omega[1]))
        }))
    }

    fn prepare_modes(
        &mut self,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(), ViscousConnectorError> {
        // C = W W^T comes from the ORIGINAL admitted passive form. Exact zero
        // factor columns do not enter the eigenproblem: a semidefinite free
        // direction cannot acquire a tiny numerical decay rate which a stiff
        // interval would amplify. No rank threshold or damping projection.
        let original = self.damping.gram_factor()?;
        let mut basis = [[0.0; 6]; 6];
        let mut count = 0;
        for j in 0..6 {
            if (0..6).any(|i| original[i][j] != 0.0) {
                for i in 0..6 {
                    basis[i][count] = original[i][j];
                }
                count += 1;
            }
        }
        if count == 0 {
            return Err(ViscousConnectorError::Invalid(
                "resolved nonzero damping factor",
            ));
        }

        // Twelve original mass-normalized body momentum coordinates factor
        // G. The small matrix W^T G W then has homogeneous units of 1/s.
        // All source resultants stay in range(W), including repeated rates
        // and cross-coupled semidefinite forms.
        let native = self.kinetic_columns()?;
        let mut kinetic = [[Dd::ZERO; 12]; 6];
        for mode in 0..count {
            for i in 0..6 {
                for coordinate in 0..12 {
                    kinetic[mode][coordinate] = kinetic[mode][coordinate]
                        + Dd::from_f64(basis[i][mode]) * native[i][coordinate];
                }
            }
        }
        let mut matrix = [0.0; 36];
        for i in 0..count {
            for j in i..count {
                let entry = finite(kinetic_product(kinetic[i], kinetic[j]).to_f64())?;
                matrix[count * i + j] = entry;
                matrix[count * j + i] = entry;
            }
        }
        let rate_scale = matrix[..count * count]
            .iter()
            .map(|x| x.abs())
            .fold(0.0_f64, f64::max);
        if rate_scale <= 0.0 {
            return Err(ViscousConnectorError::Invalid(
                "resolved nonzero damping spectrum",
            ));
        }
        for value in &mut matrix[..count * count] {
            *value /= rate_scale;
        }
        poll(cancelled)?;
        let (_, vectors) = jacobi_eigh(&matrix[..count * count], count);
        poll(cancelled)?;
        for mode in 0..count {
            let mut q = [0.0; 6];
            for i in 0..count {
                q[i] = vectors[count * i + mode];
            }
            let source = matrix_vector(basis, q)?;
            let mut response = [Dd::ZERO; 12];
            for i in 0..count {
                for coordinate in 0..12 {
                    response[coordinate] =
                        response[coordinate] + Dd::from_f64(q[i]) * kinetic[i][coordinate];
                }
            }
            // A sum of original squared kinetic factors supplies each rate.
            // Negative numerical eigenvalues are never clipped to passivity.
            let rate = finite(kinetic_product(response, response).to_f64())?;
            if rate <= 0.0 {
                return Err(ViscousConnectorError::Invalid(
                    "positive resolved active damping rate",
                ));
            }
            self.maximum_decay_rate_per_s = self.maximum_decay_rate_per_s.max(rate);
            let normalized_rate = rate / rate_scale;
            for i in 0..count {
                let mut applied = Dd::ZERO;
                let mut orthogonal = Dd::ZERO;
                for j in 0..count {
                    applied = applied + Dd::from_f64(matrix[count * i + j]) * Dd::from_f64(q[j]);
                    orthogonal =
                        orthogonal + Dd::from_f64(q[j]) * Dd::from_f64(vectors[count * j + i]);
                }
                if !normalized_rate.is_finite()
                    || (applied - Dd::from_f64(normalized_rate) * Dd::from_f64(q[i]))
                        .to_f64()
                        .abs()
                        > SPECTRAL_TOLERANCE
                {
                    return Err(ViscousConnectorError::Invalid("resolved viscous eigenpair"));
                }
                if (orthogonal.to_f64() - if mode == i { 1.0 } else { 0.0 }).abs()
                    > SPECTRAL_TOLERANCE
                {
                    return Err(ViscousConnectorError::Invalid("orthogonal viscous modes"));
                }
            }
            for i in 0..6 {
                self.source_modes[i][mode] = source[i];
                self.relative_modes[i][mode] =
                    finite(kinetic_product(native[i], response).to_f64())?;
            }
            let exposure = rate * self.duration_s;
            // Positive overflow means complete decay. In the weak limit use
            // h*phi rather than (1-exp(-lambda*h))/lambda: lambda*h can
            // underflow while the actual impulse or heat remains resolvable.
            let fraction = finite(-det::expm1(-exposure))?;
            let integrated = if exposure <= 1.0 {
                self.duration_s
                    * if exposure == 0.0 {
                        1.0
                    } else {
                        fraction / exposure
                    }
            } else {
                fraction / rate
            };
            self.integrated_decay[mode] = finite(integrated)?;
            self.heat_root[mode] = det::sqrt(integrated) * det::sqrt(1.0 - 0.5 * fraction);
            if !self.heat_root[mode].is_finite() || integrated <= 0.0 {
                return Err(ViscousConnectorError::Invalid("resolved viscous response"));
            }
        }
        Ok(())
    }

    // The same original finite-body factors define mobility for both passive
    // relaxation and simultaneous actuation. Wrench transport is prepared once
    // above; no consumer substitutes an infinite bath or diagonal inertia.
    fn kinetic_columns(&self) -> Result<[[Dd; 12]; 6], ViscousConnectorError> {
        let mut native = [[Dd::ZERO; 12]; 6];
        for i in 0..6 {
            for body in 0..2 {
                let linear = checked_vector_divide(
                    self.force_columns_world[i].scale(if body == 0 { -1.0 } else { 1.0 }),
                    det::sqrt(self.properties[body].mass()),
                )?;
                let angular = self.properties[body]
                    .angular_mobility_factor(self.angular_columns_body[i][body]);
                finite_vector(angular)?;
                let offset = 6 * body;
                native[i][offset..offset + 3].copy_from_slice(&dd(linear));
                native[i][offset + 3..offset + 6].copy_from_slice(&dd(angular));
            }
        }
        Ok(native)
    }
}

fn kinetic_product(a: [Dd; 12], b: [Dd; 12]) -> Dd {
    let mut result = Dd::ZERO;
    for i in 0..12 {
        result = result + a[i] * b[i];
    }
    result
}
fn column(matrix: Matrix, index: usize) -> [f64; 6] {
    core::array::from_fn(|i| matrix[i][index])
}
fn product(a: [f64; 6], b: Six) -> Dd {
    let mut result = Dd::ZERO;
    for i in 0..6 {
        result = result + Dd::from_f64(a[i]) * b[i];
    }
    result
}
fn matrix_vector(matrix: Matrix, vector: [f64; 6]) -> Result<[f64; 6], ViscousConnectorError> {
    rounded_six(core::array::from_fn(|i| {
        product(matrix[i], vector.map(Dd::from_f64))
    }))
}
fn rounded_six(value: Six) -> Result<[f64; 6], ViscousConnectorError> {
    let mut result = [0.0; 6];
    for i in 0..6 {
        result[i] = finite(value[i].to_f64())?;
    }
    Ok(result)
}
fn checked_vector_divide(value: Vec3, divisor: f64) -> Result<Vec3, ViscousConnectorError> {
    let result = Vec3::new(value.x / divisor, value.y / divisor, value.z / divisor);
    finite_vector(result)?;
    Ok(result)
}
